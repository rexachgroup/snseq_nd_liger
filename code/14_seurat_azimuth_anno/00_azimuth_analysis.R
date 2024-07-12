#!/usr/bin/env Rscript
library(tidyverse)
out_dir <- "../data/"
in_seurat <- "../../../analysis/pci_import/pci_seurat.rds"
in_seurat_meta <- "../../../analysis/seurat_lchen/liger_subcluster_metadata.rds"
excludes <- read_xlsx("../../../resources/subclusters_removed_byQC_final.xlsx")

plot_ts <- function() {
    cs <- system("md5sum 00_azimuth_analysis.R", intern = TRUE)
    str_glue("{cs} {base::date()}") 
}

# Ensure Seurat v4.0 or higher is installed
if (packageVersion(pkg = "Seurat") < package_version(x = "4.0.0")) {
  stop("Mapping datasets requires Seurat v4 or higher.", call. = FALSE)
}

# Ensure glmGamPoi is installed
if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    BiocManager::install("glmGamPoi")
  }
}

# Ensure Azimuth is installed
if (packageVersion(pkg = "Azimuth") < package_version(x = "0.3.1")) {
  stop("Please install azimuth - remotes::install_github('satijalab/azimuth')", call. = FALSE)
}

library(Seurat)
library(Azimuth)

# Download the Azimuth reference and extract the archive

# Load the reference
# Change the file path based on where the reference is located on your system.
reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_motorcortex")

# Load the query object for mapping
# Change the file path based on where the query file is located on your system.
#query <- LoadFileInput(path = in_seurat)
sobj <- readRDS(in_seurat)
sobj_meta <- readRDS(in_seurat_meta)
sobj_meta <- sobj_meta %>% filter(!ct_subcluster %in% excludes$ct_subcluster)

query <- subset(sobj, cells = sobj_meta$UMI)
stopifnot(all(rownames(query@meta.data) == sobj_meta$UMI))
query$ct_subcluster <- sobj_meta$ct_subcluster
query <- DietSeurat(query, assays = "RNA")
rm(sobj)
gc()

# Calculate nCount_RNA and nFeature_RNA if the query does not
# contain them already
if (!all(c("nCount_RNA", "nFeature_RNA") %in% c(colnames(x = query[[]])))) {
    calcn <- as.data.frame(x = Seurat:::CalcN(object = query))
    colnames(x = calcn) <- paste(
      colnames(x = calcn),
      "RNA",
      sep = '_'
    )
    query <- AddMetaData(
      object = query,
      metadata = calcn
    )
    rm(calcn)
}

# Calculate percent mitochondrial genes if the query contains genes
# matching the regular expression "^MT-"
if (any(grepl(pattern = '^MT-', x = rownames(x = query)))) {
  query <- PercentageFeatureSet(
    object = query,
    pattern = '^MT-',
    col.name = 'percent.mt',
    assay = "RNA"
  )
}

cells.use <- colnames(query)

# If the query contains mitochondrial genes, filter cells based on the
# thresholds for percent.mt you set in the app
if ("percent.mt" %in% c(colnames(x = query[[]]))) {
  cells.use <- cells.use & (query[["percent.mt", drop = TRUE]] <= 8 &
    query[["percent.mt", drop = TRUE]] >= 0)
}

# Filter cells based on the thresholds for nCount_RNA and nFeature_RNA
# you set in the app
cells.use <- query[["nCount_RNA", drop = TRUE]] <= 33185 &
  query[["nCount_RNA", drop = TRUE]] >= 212 &
  query[["nFeature_RNA", drop = TRUE]] <= 6486 &
  query[["nFeature_RNA", drop = TRUE]] >= 201 &
  query[["percent.mt", drop = TRUE]] <= 2

# Remove filtered cells from the query
query <- query[, cells.use]

# Preprocess with SCTransform
query <- SCTransform(
  object = query,
  assay = "RNA",
  new.assay.name = "refAssay",
  residual.features = rownames(x = reference$map),
  reference.SCT.model = reference$map[["refAssay"]]@SCTModel.list$refmodel,
  method = 'glmGamPoi',
  ncells = 2000,
  n_genes = 2000,
  do.correct.umi = FALSE,
  do.scale = FALSE,
  do.center = TRUE
)

gc()

# Find anchors between query and reference
anchors <- FindTransferAnchors(
  reference = reference$map,
  query = query,
  k.filter = NA,
  reference.neighbors = "refdr.annoy.neighbors",
  reference.assay = "refAssay",
  query.assay = "refAssay",
  reference.reduction = "refDR",
  normalization.method = "SCT",
  features = intersect(rownames(x = reference$map), VariableFeatures(object = query)),
  dims = 1:50,
  n.trees = 20,
  mapping.score.k = 100
)

gc()

# Transfer cell type labels and impute protein expression
#
# Transferred labels are in metadata columns named "predicted.*"
# The maximum prediction score is in a metadata column named "predicted.*.score"
# The prediction scores for each class are in an assay named "prediction.score.*"
# The imputed assay is named "impADT" if computed

refdata <- lapply(X = c("subclass", "cluster", "class"), function(x) {
  reference$map[[x, drop = TRUE]]
})
names(x = refdata) <- c("subclass", "cluster", "class")
if (FALSE) {
  refdata[["impADT"]] <- GetAssayData(
    object = reference$map[['ADT']],
    slot = 'data'
  )
}
query <- TransferData(
  reference = reference$map,
  query = query,
  dims = 1:50,
  anchorset = anchors,
  refdata = refdata,
  n.trees = 20,
  store.weights = TRUE
)

gc() 
# Calculate the embeddings of the query data on the reference SPCA
query <- IntegrateEmbeddings(
  anchorset = anchors,
  reference = reference$map,
  query = query,
  reductions = "pcaproject",
  reuse.weights.matrix = TRUE
)

gc() 
# Calculate the query neighbors in the reference
# with respect to the integrated embeddings
query[["query_ref.nn"]] <- FindNeighbors(
  object = Embeddings(reference$map[["refDR"]]),
  query = Embeddings(query[["integrated_dr"]]),
  return.neighbor = TRUE,
  l2.norm = TRUE
)

gc() 


# NNTransform isn't exported in Azimuth 0.4.3 -- define manually
NNTransform <- function(
  object,
  meta.data,
  neighbor.slot = "query_ref.nn",
  key = 'ori.index'
) {
  on.exit(expr = gc(verbose = FALSE))
  ind <- Indices(object[[neighbor.slot]])
  ori.index <- t(x = sapply(
    X = 1:nrow(x = ind),
    FUN = function(i) {
      return(meta.data[ind[i, ], key])
    }
  ))
  rownames(x = ori.index) <- rownames(x = ind)
  slot(object = object[[neighbor.slot]], name = "nn.idx") <- ori.index
  return(object)
}
# The reference used in the app is downsampled compared to the reference on which
# the UMAP model was computed. This step, using the helper function NNTransform,
# corrects the Neighbors to account for the downsampling.
query <- NNTransform(
  object = query,
  meta.data = reference$map[[]]
)

gc() 
# Project the query to the reference UMAP.
query[["proj.umap"]] <- RunUMAP(
  object = query[["query_ref.nn"]],
  reduction.model = reference$map[["refUMAP"]],
  reduction.key = 'UMAP_',
)


gc() 
# Calculate mapping score and add to metadata
query <- AddMetaData(
  object = query,
  metadata = MappingScore(anchors = anchors),
  col.name = "mapping.score"
)

gc() 
# VISUALIZATIONS
# 
# First predicted metadata field, change to visualize other predicted metadata
id <- c("subclass", "cluster", "class")[1]
predicted.id <- paste0("predicted.", id)

pdf(file.path(out_dir, "dimplots.pdf"), width = 30, height = 10)
# 
# DimPlot of the reference
#DimPlot(object = reference$plot, reduction = "refUMAP", group.by = id, label = TRUE) + labs(title = "ref UMAP ids")

# DimPlot of the query, colored by predicted cell type
DimPlot(object = query, reduction = "proj.umap", group.by = predicted.id, label = TRUE, pt.size = 0)

DimPlot(object = query, reduction = "proj.umap", group.by = "cluster_cell_type", label = TRUE, pt.size = 0, raster = TRUE)

DimPlot(object = query, reduction = "proj.umap", group.by = "predicted.cluster", label = TRUE, pt.size = 0) + NoLegend()
FeaturePlot(object = query, reduction = "proj.umap", features = "percent.mt", split.by = "region", raster = TRUE, label = TRUE, pt.size = 0) + NoLegend()
graphics.off()

pdf(file.path(out_dir, "dimplots_cluster.pdf"), width = 100, height = 10)
DimPlot(object = query, reduction = "proj.umap", group.by = "predicted.cluster", split.by = "predicted.subclass", label = FALSE, pt.size = 0, raster = TRUE)
graphics.off()

pdf(file.path(out_dir, "dimplots_meta.pdf"), width = 30, height = 10)
DimPlot(object = query, reduction = "proj.umap", group.by = "clinical_dx", split.by = "region", label = TRUE, pt.size = 0)
DimPlot(object = query, reduction = "proj.umap", group.by = "region", split.by = "clinical_dx", label = TRUE, pt.size = 0)

graphics.off()
# 
# # Plot the score for the predicted cell type of the query
# FeaturePlot(object = query, features = paste0(predicted.id, ".score"), reduction = "proj.umap")
# VlnPlot(object = query, features = paste0(predicted.id, ".score"), group.by = predicted.id) + NoLegend()
# 
# # Plot the mapping score
# FeaturePlot(object = query, features = "mapping.score", reduction = "proj.umap")
# VlnPlot(object = query, features = "mapping.score", group.by = predicted.id) + NoLegend()
# 
# # Plot the prediction score for the class CD16 Mono
# FeaturePlot(object = query, features = "CD16 Mono", reduction = "proj.umap")
# VlnPlot(object = query, features = "CD16 Mono", group.by = predicted.id) + NoLegend()
# 
# # Plot an RNA feature
# FeaturePlot(object = query, features = "PVALB", reduction = "proj.umap")
# VlnPlot(object = query, features = "PVALB", group.by = predicted.id) + NoLegend()
# 
# # Plot an imputed protein feature
# if (TRUE) {
#   FeaturePlot(object = query, features = "CD3-1", reduction = "proj.umap")
#   VlnPlot(object = query, features = "CD3-1", group.by = predicted.id) + NoLegend()
# }
# 
# graphics.off()

saveRDS(query, file.path(out_dir, "sobj.rds"), compress = FALSE)
write_csv(query@meta.data, file.path(out_dir, "meta.csv"))

