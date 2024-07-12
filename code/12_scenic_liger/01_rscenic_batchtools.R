liblist <- c("tidyverse","Seurat", "batchtools", "argparse", "doParallel", "SCENIC", "reticulate")
l <- lapply(liblist, function(x) {library(x, character.only = T)})

REF_PATH <- "../../resources/cisTarget_databases/"
SCENIC_DIR <- "../../analysis/12_liger_scenic/"
SOBJ_FILE <- "../../analysis/seurat/20200816/pci_filtered_seurat_object.rdat"
CELLS_FILE <- "../../resources/scenic/cells_presampled.csv"
GENES_FILE <- "../../resources/scenic/genes_presampled.csv"
RSCRIPT <- "../../conda/bin/Rscript"
use_condaenv("../../conda/")

target_regions <- c("calcarine")
target_celltypes <- c("microglia")

RESOURCES <- list(
    ncpus = 1,
    memory = 120,
    walltime = 172800,
    measure.memory = TRUE,
    nice = 100,
    chunks.as.arrayjobs = TRUE,
    partition = "bigmem"
)

main <- function() {
    # prep target sobjs
    dir.create(SCENIC_DIR)

    if (!dir.exists(batchtools)) {
        reg <- makeRegistry(file.path(SCENIC_DIR, "batchtools"), config = "../batchtools.conf.R")
    } else {
        reg <- loadRegistry(file.path(SCENIC_DIR, "batchtools"), writeable = T)
    }
    reg$packages <- liblist
    cells_tb <- read_csv(CELLS_FILE)
    genes_tb <- read_csv(GENES_FILE)
    comb_tb <- expand_grid(region = target_regions, cluster_cell_type = target_celltypes)
    comb_tb <- comb_tb %>% mutate(
        out_prefix = str_glue("{SCENIC_DIR}/{region}/{cluster_cell_type}"),
        out_expr = str_glue("{out_prefix}/exp_mat_sub_{region}_{cluster_cell_type}.rds"),
        out_genes = str_glue("{out_prefix}/target_genes_{region}_{cluster_cell_type}.tsv"),
        out_meta = str_glue("{out_prefix}/cell_info_{region}_{cluster_cell_type}.rds")
    )

    write_expr_matrix_scenic(comb_tb, sobj_file)
    gc()
    
    runs_tb <- comb_tb %>% mutate(cmd = pmap(., function(...) {
        cr <- list(...)
        str_glue("{RSCRIPT} 01a_scenic_pipeline_param.R {cr$out_expr} {cr$out_meta} {cr$out_prefix} --cis_path {REF_PATH} --nCores 32 --title {cr$region}_{cr$cluster_cell_type} --genes_file {cr$out_genes}")
    }))

    clearRegistry()
    batchMap(function(x) {
        system(x, intern = T)
    }, runs_tb$cmd)
    submitJobs(getJobTable(), resources = RESOURCES)
    waitForJobs()
    
    runs_auc_tb <- runs_tb %>%
        select(region, cluster_cell_type, run_prefix) %>%
        mutate(rss_cmd = str_glue("{RSCRIPT} 01b_scenic_pipeline_rss.R {region} {cluster_cell_type} {run_prefix}"))
    
    clearRegistry()
    batchMap(function(x) {
        system(x, intern = T)
    }, runs_auc_tb$rss_cmd)
    submitJobs(getJobTable(), resources = RESOURCES)
    waitForJobs()
}

write_expr_matrix_scenic <- function(comb_tb, sobj_file) {
    load(sobj_file)
    sobj <- nd_so
    pwalk(comb_tb, function(...) {
        cr <- list(...)
        dir.create(cr$out_prefix, recursive = T, showWarnings = F)
        cells_s <- cells_tb %>% filter(region == cr$region, cluster_cell_type == cr$cluster_cell_type)
        genes_s <- genes_tb %>% filter(region == cr$region, cluster_cell_type == cr$cluster_cell_type)
        sobj_s <- subset(sobj, cells = cells_s$cells, features = genes_s$genes)
        cell_info <- sobj_s@meta.data %>% mutate(dx_celltype = str_glue("{clinical_dx}_{cr$cluster_cell_type}")) %>% select(dx_celltype) %>% as.data.frame
        writeLines(cr$out_expr)
        saveRDS(as.matrix(GetAssayData(sobj_s, assay = "RNA", slot = "counts")), cr$out_expr, compress = F)
        writeLines(cr$out_genes)
        write.table(genes_s$genes, cr$out_genes, sep = "\t")
        writeLines(cr$out_meta)
        saveRDS(cell_info, cr$out_meta, compress = F)
    }) 
}
