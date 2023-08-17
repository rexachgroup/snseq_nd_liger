set.seed(0)
options(deparse.max.lines = 5)
liblist <- c("Seurat", "tidyverse", "readxl", "batchtools", "lme4", "lmerTest", "future.apply", "broom.mixed")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

OUT_DIR <- "../../analysis/seurat_lchen/cellbender/expr_corr/"
batchtools <- file.path(OUT_DIR, "batchtools")
SEURAT_OBJ <- "../../analysis/pci_import/pci_seurat.rds"
in_seurat_meta <- "../../analysis/pci_import/pci_seurat_meta.rds"
CELLBENDER_TB <- "../../analysis/seurat_lchen/cellbender/merge/cellbender_region.rds"
SEURAT_DGE <- "../../analysis/seurat_lchen/seurat_cluster_lme/cluster_wk.rds"

RESOURCES <- list(
    ncpus = 1,
    memory = 64,
    walltime = 60 * 60 * 1,
    measure.memory = TRUE,
    nice = 100,
    chunks.as.arrayjobs = F,
    partition = "bigmem"
)
chunk_size <- 5

# 1. Read Seurat dge
# 2. Get dx-specific gene estimate + region
# 3. Read + split cellbender meta into cluster and celltype
# 4. In batchtools worker:
#       a. Read and subset seurat obj to target genes / clsuter
#       b. Read and subset cellbender obj to target genes / cluster
#       c. cor on remaining expression matrix
main <- function() {
    dir.create(OUT_DIR)
    if (dir.exists(batchtools)) {
        reg <- loadRegistry(batchtools, conf.file = "../batchtools.conf.R", writeable = TRUE)
        sweepRegistry()
    } else {
        reg <- makeRegistry(batchtools, packages = liblist, conf.file = "../batchtools.conf.R")
    }
    reg$packages <- liblist

    seurat_dge_tb <- readRDS(SEURAT_DGE)   
    seurat_dge_tb <- seurat_dge_tb %>% 
        mutate(dge_filt = map(broom_join, filter_dx_specific_dge))

    seurat_filt_tb <- seurat_dge_tb %>%
        pluck("dge_filt") %>%
        bind_rows

    seurat_filt_tbs <- dge_extract_dx(seurat_filt_tb) %>%
        group_by(ct_cluster) %>%
        group_nest(.key = "dge_data")

    meta <- readRDS(in_seurat_meta)
    cellbender_region <- readRDS(CELLBENDER_TB) %>%
        select(region, merge_out)
    
    meta_subset <- meta %>%
        mutate( 
            barcode = str_extract(rownames, "\\w+"),
            merge_barcode = str_glue("{barcode}-1-{library_id}"),
            ct_cluster = paste(region, cluster_cell_type, sep = "-"),
            clinical_dx = fct_relevel(clinical_dx, "Control"),
            log_number_umi = log(number_umi)
        ) %>%
        mutate(
            ct_cluster = factor(ct_cluster, levels = str_sort(unique(ct_cluster), numeric = TRUE))
        )
    
    # Split by cluster and dx contrast
    cluster_wk <- contr_treatment_vec(meta_subset$clinical_dx, meta_subset$rownames) %>%
        left_join(meta_subset, by = c("rownames")) %>%
        arrange(ct_cluster) %>%
        group_by(region, cluster_cell_type, ct_cluster, other_level) %>% 
        group_nest(.key = "cell_data") %>%
        mutate(seurat_file = SEURAT_OBJ) %>%
        inner_join(cellbender_region, by = c("region")) %>%
        inner_join(seurat_filt_tbs, by = c("ct_cluster"))

    print(cluster_wk, width = Inf)
    clearRegistry()
    #batchExport(mget("cluster_worker"))
    batchMap(
        cluster_worker,
        args = cluster_wk
    )
    ids <- findNotDone() %>%
        mutate(chunk = chunk(job.id, chunk.size = chunk_size))
    cluster_wk$job_id <- getJobTable()$job.id
    submitJobs(ids, resources = RESOURCES)
    saveRDS(cluster_wk, file.path(OUT_DIR, "cluster_wk.rds"))

    cluster_wk <- cluster_wk %>% mutate(cor_mat = map(job_id, function(x) {
        loadResult(x)
    }))

    cluster_wk <- cluster_wk %>% mutate(autocor_tb = map(cor_mat, function(mat) {
        return(tibble(gene = rownames(mat), val = diag(mat)))
    }))

    cluster_wk %>%
        mutate(dx_subset = paste0(other_level, "|Control")) %>%
        select(ct_cluster, dx_subset, autocor_tb) %>%
        unnest(autocor_tb) %>%
        write_csv(file.path(OUT_DIR, "cellbender_seurat_expr_cor.csv"))
}

contr_treatment_vec <- function(vec, rownames) {
    contr <- contrasts(vec)
    reference_level <- rownames(contr)[1]
    other_level <- colnames(contr)
    
    map(other_level, function(lvl) {
        idx = which(vec == reference_level | vec == lvl)
        return(data.frame(rownames = rownames[idx]))
    }) %>% 
    setNames(nm = other_level) %>%
    bind_rows(.id = "other_level") %>%
    mutate(ref_level = reference_level)
}

cluster_worker <- function(...) {
    cr <- list(...)
    browser()
    options(future.globals.maxSize = Inf)
    cellbender_so <- readRDS(file.path(unique(cr$merge_out))) 
    common_genes <- intersect(rownames(cellbender_so), cr$dge_data$gene)
    cellbender_barcodes <- intersect(colnames(cellbender_so), cr$cell_data$merge_barcode)
    seurat_barcodes <- cr$cell_data %>% filter(merge_barcode %in% cellbender_barcodes) %>% pluck("rownames")

    subcellbender_so <- subset(cellbender_so, cells = cellbender_barcodes, features = common_genes)
    rm(cellbender_so)
    gc()
    seurat_so <- readRDS(file.path(unique(cr$seurat_file)))  
    subseurat_so <- subset(seurat_so, cells = seurat_barcodes, features = common_genes)
    rm(seurat_so)
    gc()

    assay1 <- as.matrix(t(GetAssayData(subcellbender_so)))
    assay2 <- as.matrix(t(GetAssayData(subseurat_so)))
    gene_cor <- cor(assay1, assay2)
    return(gene_cor)
}

filter_dx_specific_dge <- function(dge) {
    dge_long <- dge_extract_dx(dge)

    # genes that are significant in exactly one dx.
    dge_signif_filter <- dge_long %>%
        filter(p.value.adj < 0.1 & abs(estimate) > 0.1) %>%
        group_by(gene) %>%
        summarize(n = n()) %>%
        filter(n == 1)

    # genes that are not significant in at least two dx.
    dge_pval_nonsignif <- dge_long %>%
        filter(p.value > 0.05) %>%
        group_by(gene) %>%
        summarize(n = n()) %>%
        filter(n == 2)


    genes_dx_specific <- inner_join(dge_signif_filter, dge_pval_nonsignif, by = "gene")
    dge_f <- dge %>% filter(gene %in% genes_dx_specific$gene)
    return(dge_f)
}

dge_extract_dx <- function(dge) {
    dge %>%
        pivot_longer(
            cols = -c("gene", "model", "region", "cell_type", "ct_cluster"), 
            names_pattern = "clinical_dx(.*)") %>%
        separate(
            col = c("name"),
            sep = "\\.",
            into = c("dx", "dge_type"),
            extra = "merge"
        ) %>%
        pivot_wider(
            names_from = c("dge_type"),
            values_from = c("value")
        )
}

if (!interactive()) main()
