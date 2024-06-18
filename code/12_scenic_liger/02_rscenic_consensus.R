# Verify initial scenic run through consensus?
liblist <- c("tidyverse","Seurat", "batchtools", "argparse", "doParallel", "SCENIC", "reticulate")
l <- lapply(liblist, function(x) {library(x, character.only = T)})
#source("02a_scenic_pipeline_param.R")
old_rss = read_csv("../../analysis/12_liger_scenic/rscenic_old_old/microglia/calcarine/int/regulon_specificity_scoremicroglia_calcarine.csv")

REF_PATH <- "../../resources/cisTarget_databases/"
SCENIC_DIR <- "../../analysis/12_liger_scenic/rscenic_consensus/"
SOBJ_FILE <- "../../analysis/seurat/20200816/pci_filtered_seurat_object.rdat"
CELLS_FILE <- "../../resources/scenic/cells_presampled.csv"
GENES_FILE <- "../../resources/scenic/genes_presampled.csv"
RSCRIPT <- "../../conda/bin/Rscript"
use_condaenv("../../conda/")

target_regions <- c("calcarine")
target_celltypes <- c("microglia")
RUN_CT <- 10
batchtools <- file.path(SCENIC_DIR, "batchtools")

main <- function() {
    # prep target sobjs
    load(SOBJ_FILE)
    sobj <- nd_so
    dir.create(SCENIC_DIR)

    if (!dir.exists(batchtools)) {
        reg <- makeRegistry(file.path(SCENIC_DIR, "batchtools"))
    } else {
        reg <- loadRegistry(file.path(SCENIC_DIR, "batchtools"), writeable = T)
    }
    reg$cluster.functions <- makeClusterFunctionsMulticore(ncpus = 1)
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

    # subset
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

    # format run
    run_iter <- tibble(run = 1:RUN_CT)
    runs_tb <- cross_join(comb_tb, run_iter) %>%
        mutate(run_prefix = str_glue("{out_prefix}/run_{run}/"))

    runs_tb <- runs_tb %>% mutate(cmd = pmap(., function(...) {
        cr <- list(...)
        str_glue("{RSCRIPT} 02a_scenic_pipeline_param.R {cr$out_expr} {cr$out_meta} {cr$run_prefix} --cis_path {REF_PATH} --nCores 32 --title {cr$region}_{cr$cluster_cell_type} --genes_file {cr$out_genes}")
    }))
    saveRDS(runs_tb, file.path(SCENIC_DIR, "runs.rds"))
    runs_tb <- readRDS(file.path(SCENIC_DIR, "runs.rds"))

    # batchtools SCENIC script
    clearRegistry()
    batchMap(function(x) {
        system(x, intern = T)
    }, runs_tb$cmd)
    submitJobs(getJobTable())
    waitForJobs()

    # do rss processing
    runs_auc_tb <- runs_tb %>%
        select(region, cluster_cell_type, run_prefix) %>%
        mutate(rss_cmd = str_glue("{RSCRIPT} 02b_scenic_pipeline_rss.R {region} {cluster_cell_type} {run_prefix}"))

    clearRegistry()
    batchMap(function(x) {
        system(x, intern = T)
    }, runs_auc_tb$rss_cmd)
    submitJobs(getJobTable())
    waitForJobs()

    # load rss
    runs_rss_tb <- runs_tb %>%
        select(region, cluster_cell_type, run, run_prefix) %>%
        mutate(run = fct_relevel(as.character(run), str_sort(as.character(run), numeric = T))) %>%
        mutate(rss_file = str_glue("{run_prefix}/int/regulon_specificity_score_{cluster_cell_type}_{region}.csv")) %>%
        mutate(rss_rank = map(rss_file, read_csv))

    # extract tf order
    runs_rss_tb <- runs_rss_tb %>%
        mutate(rss_rank_tb = map(rss_rank, extract_tf_order))
    dx_ranking_cols <- str_subset(colnames(runs_rss_tb$rss_rank_tb[[1]]), "rank")

    # plot tf order comparison
    dx_ranking_plots <- map(dx_ranking_cols, function(dx_nm) {
        n = 10
        print(str_glue("dx_rank_comparision {dx_nm} {n}"))
        rank_tb <- rss_ranking_match(runs_rss_tb, old_rss, dx_nm, n)
        ggplot(rank_tb, aes(x = run, y = fct_rev(TF), fill = .data[[dx_nm]], label = .data[[dx_nm]])) +
            geom_tile() +
            scale_fill_viridis_c(direction = -1) +
            geom_text() +
            labs(
                title = str_glue("old {dx_nm} RSS rank vs. repeat RSS rank"),
                subtitle = str_glue("{runs_rss_tb$cluster_cell_type} {runs_rss_tb$region}")
            )
    })
    pdf(file.path(SCENIC_DIR, "rss_ranking_comparison.pdf"), width = 10, height = 8)
    print(dx_ranking_plots)
    graphics.off()

}

extract_tf_order <- function(x, ct) {
    x %>% 
        select(...1, contains("rank")) %>%
        rename(rownames = ...1) %>%
        mutate(TF = str_extract(rownames, "\\w+")) %>%
        mutate(TF = str_remove(TF, "_extended"))
}

rss_ranking_match <- function(runs_rss_tb, old_rss, col, n = 10) { 
    # get the top RSS ranking for the old rss table
    old_tf_head <- old_rss %>% extract_tf_order %>%
        arrange(.data[[col]]) %>%
        slice_head(n = n) %>%
        pluck("TF")
    # generate plot label mapping
    old_tf_label <- paste0("old_rank", 1:10, " | ", old_tf_head) %>%
        as.factor()
    levels(old_tf_label) <- str_sort(old_tf_label, numeric = T)
    names(old_tf_label) <- old_tf_head

    # match new rss table TF names and fetch rank number
    # break multiple extended regulons by min rank?
    tf_rank_tb <- runs_rss_tb %>%
        mutate(order_f_tb = map(rss_rank_tb, function(x) {
            xo <- x %>% filter(TF %in% old_tf_head)
            xo %>% select(TF, any_of(col)) %>%
                group_by(TF) %>%
                summarize(.data[[col]] = min(.data[[col]]))
        })) %>%
        select(run, order_f_tb) %>%
        unnest(order_f_tb) %>%
        filter(!is.na(TF)) %>%
        mutate(
            TF = old_tf_label[TF]
        )

    return(tf_rank_tb)
}
