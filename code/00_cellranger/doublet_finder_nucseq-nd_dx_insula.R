# doublet_finder_nucseq-nd_dx_region.R

# 2020-07-22
# For nucseq-nd project
# Running for each of the 3 brain regions
# Have seurat objects already split by region
# Splitting by clinical_dx to run doublet finder proc


### options -----------------------------------------------------------------------------


### Requirements ------------------------------------------------------------------------
library(readxl)
library(cowplot)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(Seurat)
library(DoubletFinder)
library(reticulate)
# reticulate::use_python("/u/local/apps/python/3.7.2/bin/python3", required=TRUE)

### Functions ---------------------------------------------------------------------------
ggColorHue <- function(n) {
  hues = seq(15, 375, length = n + 2)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

createSubDirs <- function(out_dir, tar_dirs) {
  # create sub directories for for graphical output
  # Structure like so
  # tar_dir
  #    |
  #    |-> clustering
  #    |-> metadata
  #    |-> markers
  sub_dirs <- c("clustering/", "metadata/", "markers/")
  for (tar_dir in tar_dirs) {
    for (sub_dir in sub_dirs) {
      full_dir <- paste0(out_dir, tar_dir, sub_dir)
      if (!dir.exists(full_dir)) {
        dir.create(full_dir, recursive = TRUE)
      }
    } 
  }
}

createOuputDirs <- function(tar_dirs) {
  # Create master directory structure for graphical outputs
  out_dir <- paste0(format(Sys.Date(), "%Y_%m_%d"), "/")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive=FALSE)
  }

  tar_dirs <- paste(tar_dirs, "/", sep = "")
  createSubDirs(out_dir, tar_dirs)
  
  out_dir
}

printPlotToFile <- function(ggp, file, width = 10, height = 8) {
  png(file, width=width, height=height, res=200, units="in")
  print(ggp)
  dev.off()
  message(file, "\t...written to file")
} 

doSeuratClustering <- function(
  seurat_object, 
  out_dir,
  dx_str,
  start_res=0.4,
  end_res=1.0
) {
  for(i in seq(from=start_res,to=end_res, by=0.1)){
    print("")
    print(paste0("clustering for resolution: ", i))
    seurat_object <- FindClusters(
      object=seurat_object, 
      resolution = i, 
      n.start = 100
    )
    
    clustering_name <- paste0("RNA_snn_res.", i)
    print(unique(seurat_object@meta.data[ ,clustering_name]))
    
    png(
      paste0(out_dir,"clustering/", dx_str, "_clusters_res_",i,"_umap.png" ), 
      width=10, 
      height=8, 
      res=300, 
      units="in"
    )
    print(
      DimPlot(
        seurat_object, 
        reduction="umap", 
        label=TRUE, 
        pt.size=0.05, 
        group.by=clustering_name
      )
    )
    dev.off()
    
  }
  
  seurat_object 
}

expGraphGGP <- function(graph.df, graph.type, graph.title, tar.column = "exp") {
  # Function to graph coexpression of marker genes in dropseq dataset
  # THis is best for graphing expression
  ggp <- ggplot(graph.df, aes_string(x=paste0(graph.type, "_1"), 
                                     y=paste0(graph.type, "_2"))) +
    geom_point(aes_string(color=tar.column), shape=20,size=PT_SIZE) +
    # scale_color_distiller(name="expression\n", direction=1, type="div", 
    #                       palette="Reds",breaks=c(0,3), limits=c(0,3)) + 
    scale_color_gradientn(name=paste0(tar.column, "\n"), colours = c(
      "lightgrey", "#fee090", "#fdae61", "#f46d43", "#ca0020")) +
    ggtitle(graph.title) +
    theme_classic() +
    theme(legend.position="bottom",
          axis.ticks=element_blank(),
          axis.text=element_blank(),
          axis.title=element_text(size=12))
  if (graph.type=="tSNE") {
    ggp + xlab("tSNE_1") + ylab("tSNE_2") 
  } else {
    ggp + xlab("UMAP_1") + ylab("UMAP_2")
  }
  
  ggp
}

umapGraphGGP <- function(g.df, tar.col, col.legend, g.title) {
  # Similar to above but only generates a uamp plot ggplot object
  ggplot(g.df, 
         aes_string(x=paste0("UMAP_1"), 
                    y=paste0("UMAP_2"))) +
    geom_point(aes_string(color=tar.col), shape=20,size=PT_SIZE) +
    scale_color_distiller(name=col.legend, direction=-1, palette="Spectral") +
    ggtitle(g.title) +
    xlab("UMAP_1") + 
    ylab("UMAP_2") +
    theme_classic() +
    theme(axis.ticks=element_blank(),
          axis.text=element_blank(),
          axis.title=element_text(size=12))
}

doMetadataGraphs <- function(seurat_object, tar_metadata, target_str, out_dir) {
  # Output UMAP graphs of target meta data for seurat object
  cat("\n...Generating metadata plots\n")
  metadata_dir <- paste0(out_dir, "metadata/")
  
  for(md in tar_metadata){
    
    out_md_file <- paste0(metadata_dir, target_str, "_", md, "_umap.png")
    printPlotToFile(
      DimPlot(
        seurat_object, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = PT_SIZE, 
        group.by = md
      ),
      out_md_file
    )
  }
  
  umap_coords  <- Embeddings(seurat_object, reduction="umap")
  # UMI UMAP graph
  umi_graph_df <- as.data.frame(
    cbind(umap_coords, seurat_object@meta.data[ ,"nCount_RNA"])
  )
  colnames(umi_graph_df)[3] <- "umi"
  #  printPlotToFile(ggp, file) 
  umi_ggp <- umapGraphGGP(
    umi_graph_df, 
    tar.col="umi",
    col.legend="UMI",
    g.title="UMI counts")
  
  umi_file <- paste0(metadata_dir, target_str,"_umi_counts_umap.png")
  printPlotToFile(umi_ggp, umi_file)
  
  # Gene detection UMAP graph
  gd_umap_df <- as.data.frame(
    cbind(umap_coords, seurat_object@meta.data[ ,"nFeature_RNA"])
  )
  colnames(gd_umap_df)[3] <- "gene_det"
  gd_graph <- umapGraphGGP(
    gd_umap_df, 
    tar.col = "gene_det",
    col.legend = "gene detected",
    g.title = "Gene detection counts"
  )
  
  gd_file <- paste0(
    metadata_dir, target_str,"_gene_detection_counts_umap.png"
  )
  printPlotToFile(gd_graph, gd_file)
  
  # percent mitochondrial
  mito_umap_df <- as.data.frame(
    cbind(umap_coords, seurat_object@meta.data[ ,"percent.mt"])
  )
  colnames(mito_umap_df)[3] <- "per_mito"
  mito_graph <- umapGraphGGP(
    mito_umap_df, 
    tar.col = "per_mito",
    col.legend = "percent mito",
    g.title = "Percent Mitochondrial"
  )
  
  mito_file <- paste0(
    metadata_dir, target_str, "_percent_mito_umap.png"
  )
  printPlotToFile(mito_graph, mito_file)
  
}

doMarkersGraphs <- function(
  seurat_object, 
  marker_genes, 
  target_str, 
  is_mouse = FALSE,
  out_dir
) {
  # Generate and output gene marker graph expression on UMAP reductions
  cat("\n...Generating marker plots\n")
  markers_dir <- paste0(out_dir, "markers/")
  umap_coords <- Embeddings(seurat_object, reduction="umap")
  combined_markers_ggp_list <- list()
  
  if(is_mouse) {
    # Make sure genes are all lowercase for mouse
    marker_genes$gene_symbol <- sapply(marker_genes$gene_symbol, function(x) simpleCap(tolower(x)))
  }
  
  for(marker_type in unique(marker_genes$marker_for)) {
    cat("\t...", marker_type,"\n")
    subset_markers <- marker_genes[marker_genes$marker_for==marker_type, ]
    
    marker_exp <- FetchData(
      seurat_object,
      vars = subset_markers$gene_symbol, 
      slot = "data"
    )
    # output combined UMAP plots for each gene in marker class
    umap_marker_ggp_list <- lapply(
      colnames(marker_exp), 
      function(gene) {
        temp_df <- as.data.frame(cbind(umap_coords, marker_exp[ ,gene]))
        colnames(temp_df)[3] <- "exp"
        temp_df$exp[(temp_df$exp >= 3)] <- 3
        temp_df$exp[(temp_df$exp < 0)] <- 0
        
        expGraphGGP(temp_df, graph.type = "UMAP", graph.title =  gene)
      })
    # Dynamically generate graph height by # of umap graphs
    out_h <- ceiling(length(umap_marker_ggp_list)/3) * 533.3
    marker_genes_out_file <- paste0(
      markers_dir, target_str, "_", marker_type, "_umap.png"
    )
    
    png(
      marker_genes_out_file,
      units='px', 
      width=1600, 
      height=out_h
    )
    print(do.call("grid.arrange", c(umap_marker_ggp_list, ncol=3)))
    dev.off()
    cat("\n",marker_genes_out_file, "\t...written to file\n")
    
    # Save graph of mean marker expression for combined graphs of all markers
    umap_marker_exp_df <- as.data.frame(cbind(umap_coords, rowMeans(marker_exp)))
    colnames(umap_marker_exp_df)[3] <- "exp"
    umap_marker_exp_df$exp[(umap_marker_exp_df$exp >= 3)] <- 3
    umap_marker_exp_df$exp[(umap_marker_exp_df$exp < 0)] <- 0
    
    umap_marker_exp_ggp <- expGraphGGP(
      umap_marker_exp_df,
      graph.type = "UMAP",
      graph.title = marker_type
    )
    combined_markers_ggp_list[[marker_type]] <- umap_marker_exp_ggp
    
  }
  
  combined_markers_out_file <- paste0(
    markers_dir, target_str,"_combined_markers_mean_exp_umap.png"
  )
  
  out_h2 <- ceiling(length(combined_markers_ggp_list)/3) * 533.3
  
  png(combined_markers_out_file, units='px', width=1600, height=out_h2)
  print(do.call("grid.arrange", c(combined_markers_ggp_list, ncol=3)))
  dev.off()
  cat("\n\n",combined_markers_out_file, "\t...written to file\n")
  
}

runSeuratPreProcess <- function(
  s_obj, 
  target_str, 
  target_out_dir,
  pc_dims = 70, 
  do_clustering = FALSE
) {
  # Target string is the name of data of the seurat object, for example "AD"
  # Run through Seurat Pipeline to process data for doubletFinder
  # Takes subsetted object from main object and re-processes
  # Returns a seurat object
  s_obj <- FindVariableFeatures(
    s_obj, 
    selection.method = "vst", 
    nfeatures = 2000
  )
  s_obj <- ScaleData(s_obj, features = VariableFeatures(s_obj))
  s_obj <- RunPCA(
    s_obj, 
    features = VariableFeatures(object = s_obj), 
    npcs = pc_dims
  )
  
  elbow_file <- paste0(target_out_dir, target_str,"_pc_elbow_plot.png")
  png(
    elbow_file,  
    width=10, height=8, res=300, units="in"
  )
  print(ElbowPlot(s_obj, ndims = pc_dims))
  dev.off()
  message(elbow_file, "\t...written to file")
  
  s_obj <- FindNeighbors(
    object = s_obj, 
    reduction = "pca", 
    dims = 1:pc_dims, 
    nn.eps = 0, 
    k.param = 30
  )
  
  s_obj <- RunUMAP(object = s_obj, dims = 1:pc_dims)
  if (do_clustering) {
    s_obj <- doSeuratClustering(
      s_obj, 
      out_dir = target_out_dir, 
      dx_str = target_str,
      start_res = 0.3,
      end_res = 0.4
    )
  }
  s_obj
}

### Procedure ---------------------------------------------------------------------------
main <- function(so_list) {
  set.seed(0xABCDEF)
  marker_genes <- read.csv(MARKERS_FILE)
  # create output dirs
  dx <- sapply(region_list, function(r_obj) unique(r_obj@meta.data$clinical_dx))
  region <- sapply(region_list, function(r_obj) unique(r_obj@meta.data$region))
  # dx <- as.vector(unique(s_obj@meta.data$clinical_dx))
  # region <- unique(s_obj@meta.data$region)
  tar_dirs <- paste(region, dx, sep="_")
  out_dir <- createOuputDirs(tar_dirs)
  
  # i <- 1
  for (i in 1:length(so_list)) {
    tar_so <- region_list[[i]]
    # tar_so is the target seurat object, in this case the seurat object split by clinical dx
    dx_str <- tar_dirs[i]
    dx_dir <- paste0(out_dir, tar_dirs[i], "/")
    message("\n", dx_str)
    message("...begining process")
    
    # function(s_obj, target_str, target_out_dir, pc_dims = 70) 
    tar_so <- runSeuratPreProcess(
      tar_so, 
      target_str = dx_str,
      target_out_dir = dx_dir,
      pc_dims = NUM_PCs,
      do_clustering = TRUE
    )
    
    ## DoubletFinder Procedure ----------------------------------------------------------
    so_sweep_list <- paramSweep_v3(tar_so, PCs=1:NUM_PCs, sct=FALSE)
    save(
      so_sweep_list, 
      file = paste0(dx_dir, dx_str,"_parameter_sweep_out_",NUM_PCs,"_pcs.rdata")
    )
    
    sweep_stats <- summarizeSweep(so_sweep_list, GT=FALSE) 
    bcmvn_s <- find.pK(sweep_stats)
    # get pk  as.numeric(as.character(bcmvn_s[which.max(bcmvn_s$BCmetric), "pK"]))
    tar_pk <- as.numeric(as.character(bcmvn_s[which.max(bcmvn_s$BCmetric), "pK"]))
    
    bcmvn_s$BCmetric <- as.numeric(bcmvn_s$BCmetric)
    bcmvn_s$pK <- as.numeric(as.character(bcmvn_s$pK))
    # output graph
    printPlotToFile(
      ggplot(bcmvn_s, aes(x = pK, y = BCmetric)) +
        geom_point(color = "cyan3") +
        geom_line(color = "cyan3") +
        geom_vline(xintercept = tar_pk, linetype = "dashed", color = "red", size = 0.5) +
        xlab("pK") +
        ylab("BCmvn") + 
        theme_classic(),
      file = paste0(dx_dir, dx_str, "_pK_BCmvn_plot.png")
    )
    
    # Initial Run
    nExp_poi <- round(DOUBLET_RATE*ncol(tar_so))  
    tar_so <- doubletFinder_v3(
      tar_so, 
      PCs = 1:NUM_PCs, 
      pN = TAR_PN,
      pK = tar_pk, 
      nExp = nExp_poi, 
      reuse.pANN = FALSE,
      sct = FALSE
    )
    # With homotypic modeled
    annotations <- tar_so@meta.data$cluster_cell_type
    # annotations <- annotations[annotations != "mixed"]
    homotypic_prop <- modelHomotypic(annotations)
    nExp_poi.adj <- round(nExp_poi*(1-homotypic_prop))
    
    pANN_column <- paste0("pANN_",TAR_PN,"_", tar_pk, "_", nExp_poi)
    tar_so <- doubletFinder_v3(
      tar_so, 
      PCs = 1:NUM_PCs, 
      pN = TAR_PN,
      pK = tar_pk, 
      nExp = nExp_poi.adj, 
      reuse.pANN = pANN_column,
      sct = FALSE
    )
    
    # doublet_col <- grep("^DF", colnames(tar_so@meta.data), value = TRUE)
    doublet_col <- paste0("DF.classifications_", TAR_PN, "_", tar_pk, "_", nExp_poi)
    doublet_col_hc <- paste0("DF.classifications_", TAR_PN,"_", tar_pk, "_", nExp_poi.adj)
    
    tar_so@meta.data$doublets <- sapply(1:nrow(tar_so@meta.data), function(r) {
      if (tar_so@meta.data[r, doublet_col] == "Singlet" &
          tar_so@meta.data[r, doublet_col_hc] == "Singlet") {
       "Singlet"
      } else if (tar_so@meta.data[r, doublet_col_hc]  == "Doublet") {
        "Doublet_hc"
    } else {
        "Doublet_lc"
      }
    })
    
    tar_so@meta.data$ct_w_doublet <- sapply(1:nrow(tar_so@meta.data), function(r) {
      if (tar_so@meta.data[r, doublet_col] == "Singlet") {
        tar_so@meta.data[r, "cluster_cell_type"]
      } else {
        tar_so@meta.data[r, "doublets"]
      }
    })
    
    gg_df <- as.data.frame(Embeddings(tar_so, reduction="umap"))
    gg_df$doublets <- as.character(tar_so@meta.data$doublets)
    gg_df$doublets <- factor(gg_df$doublets, levels = c("Singlet", "Doublet_lc", "Doublet_hc"))
    
    printPlotToFile(
      ggplot(gg_df %>% arrange(doublets), aes(x=UMAP_1,y=UMAP_2)) +
        # ggplot(gg_df, aes(x=UMAP_1,y=UMAP_2)) +
        geom_point(aes(color=doublets), shape=20,size=PT_SIZE) +
        scale_color_manual(
          breaks = c("Singlet", "Doublet_lc", "Doublet_hc"),
          labels = c("Singlet", "Doublet-Low Confidence", "Doublet-High Confidence"),
          values =  c("grey70", "orange1","red3")) + 
        xlab("UMAP_1") + 
        ylab("UMAP_2") +
        theme_classic() +
        theme(
          axis.ticks=element_blank(),
          axis.text=element_blank(),
          axis.title=element_text(size=12)) +
        guides(colour = guide_legend(override.aes = list(size=5))),
      file = paste0(dx_dir, dx_str, "_doublet_v_singlet_umap_out.png")
    )
    

    gg.cols <- ggColorHue(length(unique(tar_so@meta.data$cluster_cell_type)))
    gg_df$ct_w_doublet <- as.character(tar_so@meta.data$ct_w_doublet)
    celltypes <- sort(unique(tar_so@meta.data$cluster_cell_type))
    gg_df$ct_w_doublet <- factor(gg_df$ct_w_doublet, levels = rev(c("Doublet_hc", "Doublet_lc", celltypes)))
    
    printPlotToFile(
      ggplot(gg_df %>% arrange(ct_w_doublet), aes(x=UMAP_1,y=UMAP_2)) +
        geom_point(aes(color=ct_w_doublet), shape=20,size=PT_SIZE) +
        scale_color_manual(
          breaks = c(celltypes, "Doublet_lc", "Doublet_hc"),
          labels = c(celltypes, "Doublet-Low Confidence", "Doublet-High Confidence"),
          values = c(gg.cols, "orange1", "red3")) + 
        xlab("UMAP_1") + 
        ylab("UMAP_2") +
        theme_classic() +
        theme(
          axis.ticks=element_blank(),
          axis.text=element_blank(),
          axis.title=element_text(size=12)) +
        guides(colour = guide_legend(override.aes = list(size=5))),
      file = paste0(dx_dir, dx_str, "_doublet_v_cell_type_umap_out.png")
    )
    
    # save doublets to csv
    doublet_df <- as.data.frame(tar_so@meta.data[ , c(pANN_column, "doublets")])
  
    doublet_file <- paste0(
      dx_dir, dx_str,
      "_doublets_classifications_doubletFinder.csv"
    )
    write.csv(doublet_df, file = doublet_file)
    
    ## Graphical output on umaps --------------------------------------------------------
    # Metadata to graph on umap
    doMetadataGraphs(
      tar_so, 
      tar_metadata = MD_GROUPS, 
      target_str = dx_str,
      out_dir = dx_dir
    ) 

    
    doMarkersGraphs(
      tar_so, 
      marker_genes = marker_genes, 
      target_str = dx_str,
      out_dir = dx_dir
    )
    
    saveRDS(tar_so, file = paste0(dx_dir, dx_str, "_so.rds"))
    
    message("\n")
    rm(tar_so)
    gc()
    
  }
  message("Process Complete!")
}

# memory of 160G was enough

######## Global 
MARKERS_FILE <-
  "/u/project/geschwind/dpolioud/nucseq_nd/resources/20200703_cell_markers_refined.csv"
PT_SIZE <- 0.5
NUM_PCs <- 40
DOUBLET_RATE <- 0.076
TAR_PN <- 0.25
MD_GROUPS <- c(
  "library_id",
  "clinical_dx", 
  "cluster_cell_type",
  "finalsite",
  "age",
  "sex"
)

## Data & dirs --------------------------------------------------------------------------

# seurat object: s_obj
# tar_file <- "../calcarine_nd/2020_07_21/calcarine_nucseq-nd_no_v3_seurat_object.rds"
tar_file <- "../insula_nd/2020_07_20/insula_nucseq-nd_seurat_object.rds"
# tar_file <- "../precg_nd/2020_07_28/precg_nucseq-nd_seurat_object.rds"
s_obj <- readRDS(tar_file)

s_obj@meta.data[!duplicated(s_obj@meta.data$library_id), c("library_id", "region", "clinical_dx")]
# fix insula/midinsula
# s_obj@meta.data$region <- gsub("midInsula", "insula", s_obj@meta.data$region)

# split by region
region_list <- SplitObject(s_obj, split.by = "clinical_dx") 

# for(e in region_list) { print(unique(e@meta.data$clinical_dx))}
# # create output dirs
# dx <- as.vector(unique(s_obj@meta.data$clinical_dx))
# region <- unique(s_obj@meta.data$region)
# tar_dirs <- paste(region, dx, sep="_")
# out_dir <- createOuputDirs(tar_dirs)

rm(s_obj)
gc()

main(region_list)





#########################################################################################
# sapply(ls(), function(x) object.size(get(x)))


# Pre-process with Seurat piplines
# clinical_so <- FindVariableFeatures(
#   clinical_so, 
#   selection.method = "vst", 
#   nfeatures = 2000
# )
# clinical_so <- ScaleData(clinical_so)
# clinical_so <- RunPCA(
#   clinical_so, 
#   features = VariableFeatures(object = clinical_so), 
#   npcs = 100
# )
# 
# elbow_file <- paste0(dx_dir, dx_str,"_pc_elbow_plot.png")
# png(
#   elbow_file,  
#   width=10, height=8, res=300, units="in"
# )
# print(ElbowPlot(clinical_so, ndims = 100))
# dev.off()
# message(elbow_file, "\t...written to file")
# 
# clinical_so <- FindNeighbors(
#   object = clinical_so, 
#   reduction = "pca", 
#   dims = 1:70, 
#   nn.eps = 0, 
#   k.param = 30
# )
# 
# clinical_so <- RunUMAP(object = clinical_so, dims = 1:70)
# clinical_so <- doSeuratClustering(
#   clinical_so, 
#   out_dir = dx_dir, 
#   dx_str = dx_str,
#   start_res = 0.3,
#   end_res = 0.5
# )

# for(so in clinical_so_list) {print(unique(so@meta.data$clinical_dx))}
# psp_so <- clinical_so_list[[1]]
# psp_chems_list <- SplitObject(psp_so, split.by = "library_chemistry")

# for(so in psp_chems_list) {print(unique(so@meta.data$library_chemistry)) }
