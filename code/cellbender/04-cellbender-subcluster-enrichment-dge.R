set.seed(0)
options(deparse.max.lines = 5)
liblist <- c("Seurat", "tidyverse", "readxl", "batchtools", "lme4", "lmerTest", "future.apply", "broom.mixed")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

cellbender_tb <- "../../analysis/seurat_lchen/cellbender/merge/cellbender_region.rds"
SUBCLUSTER_CB_F <- "../../analysis/seurat_lchen/subcluster_excluded_dge/seurat_filter/02_subcluster_cellbender_filtered_meta.csv"
OUT_DIR <- "../../analysis/seurat_lchen/cellbender/subcluster_enrichment_dge/"
batchtools <- file.path(OUT_DIR, "batchtools")

RESOURCES <- list(
    ncpus = 1,
    memory = 32,
    walltime = 60 * 60 * 4,
    measure.memory = TRUE,
    nice = 100,
    chunks.as.arrayjobs = F,
    partition = "bigmem"
)
prop_detected_filter <- 0.1
chunk_size <- 5
model_designs <- c("expression ~ clinical_dx + pmi + age + sex + number_umi + percent_mito + (1 | library_id)")

main <- function() { 
    dir.create(OUT_DIR)
    if (dir.exists(batchtools)) {
        reg <- loadRegistry(batchtools, conf.file = "../batchtools.conf.R", writeable = TRUE)
        sweepRegistry()
        reg$packages <- liblist
    } else {
        reg <- makeRegistry(batchtools, packages = liblist, conf.file = "../batchtools.conf.R")
    }
    meta <- read_csv(SUBCLUSTER_CB_F)
    cellbender_region <- readRDS(cellbender_tb) %>%
        select(region, merge_out)
    
    meta_subset <- meta %>%
        left_join(cellbender_region, by = "region") %>%
        mutate(filepath = merge_out) %>%
        mutate( 
            barcode = str_extract(UMI, "\\w+"),
            merge_barcode = str_glue("{barcode}-1-{library_id}"),
            ct_subcluster = paste(region, cluster_cell_type, liger_clusters, sep = "-"),
            clinical_dx = fct_relevel(clinical_dx, "Control"),
            log_number_umi = log(number_umi)
        ) %>%
        mutate(
            ct_subcluster = factor(ct_subcluster, levels = str_sort(unique(ct_subcluster), numeric = TRUE))
        ) %>%
        filter(
            cluster_cell_type == cell_type,
            cluster_cell_type == "astrocyte",
        )
    
    model_designs <- c("expression ~ enrichment_cluster + pmi + age + sex + number_umi + percent_mito + (1 | library_id)")

    cluster_tbs <- meta_subset %>%
        inner_join(tibble(model_design = model_designs), by = character()) %>%
        arrange(ct_subcluster, model_design) %>%
        group_by(region, cluster_cell_type, model_design) %>%
        group_nest(keep = TRUE)

    subcluster_wk <- meta_subset %>%
        group_by(region, cluster_cell_type, ct_subcluster, liger_clusters) %>%
        group_nest(keep = TRUE) %>%
        inner_join(cluster_tbs, by = c("region", "cluster_cell_type")) %>%
        select(-data.x) %>% rename(data = data.y)
    
    clearRegistry()
    batchExport(mget(ls()))
    
    ids <- batchMap(subcluster_worker,
        args = list(test_cluster = subcluster_wk$liger_clusters, meta = subcluster_wk$data),
        more.args = list(prop_detected_filter = prop_detected_filter)
    )
    ids <- findNotDone() %>%
        mutate(chunk = chunk(job.id, chunk.size = chunk_size))
    subcluster_wk$job.id <- getJobTable()$job.id
    subcluster_wk$prop_detected_filter <- prop_detected_filter
    setJobNames(subcluster_wk$job.id, as.character(subcluster_wk$ct_subcluster))
    submitJobs(ids, resources = RESOURCES)
    saveRDS(subcluster_wk, file.path(OUT_DIR, "subcluster_wk.rds"), compress = FALSE)
    waitForJobs()
    
    subcluster_wk <- subcluster_wk %>%
        mutate(broom_long = pmap(., function(...) {
            cr <- list(...)
            tryCatch({
                broom_list <- loadResult(cr$job.id)
                writeLines(paste(cr$job.id))
                format_lm_long_output(broom_list, cr$data, cr$model_design, "enrichment_cluster")
            }, 
            error = function(x) { print(x); return(NA) })
        }))

    subcluster_wk <- subcluster_wk %>%
        mutate(broom_join = pmap(list(data, model_design, job.id), function(data, model_design, job.id) {
                tryCatch({
                    broom_list <- loadResult(job.id)
                    writeLines(paste(job.id))
                    format_lm_output(broom_list, data, model_design, "enrichment_cluster")
                }, 
                error = function(x) { print(x); return(NA) })
            }))
    gc()

    subcluster_wk %>% 
        select(-data) %>%
        saveRDS(file.path(OUT_DIR, "subcluster_lme_long.rds"), compress = FALSE)
    saveRDS(subcluster_wk, file.path(OUT_DIR, "subcluster_wk.rds"), compress = FALSE)

    # Merge all subcluster results per region / celltype and write out as separate csv files. 
    subcluster_wk %>%
        filter(!is.na(broom_join)) %>%
        group_by(region, cluster_cell_type) %>%
        group_split %>%
        walk(function(x) {
            join_cols <- c("gene", "cell_type", "model", "region")
            filepath <- file.path(
                OUT_DIR,
                "enrichment_tables",
                str_glue("{unique(x$region)}-{unique(x$cluster_cell_type)}.csv")
            )
            broom_tbl <- reduce(x$broom_join, full_join, by = join_cols)
            writeLines(as.character(filepath))
            dir.create(dirname(filepath), recursive = TRUE, showWarnings = FALSE)
            write_csv(broom_tbl, filepath)
        })

    subcluster_wk <- subcluster_wk %>%
        mutate(broom_join_noterm = pmap(., function(...) {
                cr <- list(...)
                gc()
                tryCatch({
                    broom_list <- loadResult(cr$job.id)
                    writeLines(paste(cr$job.id))
                    fmt_lm_output_noterm(broom_list, cr$data, cr$model_design, "enrichment_cluster")
                }, 
                error = function(x) { print(x); return(NA) })
            }))

    
}

subcluster_worker <- function(test_cluster, meta, prop_detected_filter) {
    # Load in subcluster file.
    options(future.globals.maxSize = Inf)
    sobj <- readRDS(unique(meta$filepath))
    
    # Subset.
    subcluster_so <- subset(sobj, cells = meta$merge_barcode)
    gc()

    print(dim(subcluster_so))
    # Calculate proportion detected.
    prop_detected_tbl <- prop_detected_generate(meta, subcluster_so)
    detected <- prop_detected_tbl %>% filter(prop_detected_all > prop_detected_filter)

    subcluster_so <- subset(subcluster_so, cells = meta$merge_barcode, features = detected$gene)
    print(dim(subcluster_so))

    # Run lmer_enrichment_de
    return(run_lmer_enrichment_de(subcluster_so, meta, unique(meta$model_design), test_cluster, down_sample_cells = 10000))
}

run_lmer_enrichment_de <- function(
    seurat_obj,
    meta,
    model_design,
    test_cluster,
    down_sample_cells = NULL,
    cores = NULL){

    # Downsample if down_sample_cells defined.
    if (!is.null(down_sample_cells) && ncol(seurat_obj) > down_sample_cells) {
        writeLines(str_glue("downsampling to {down_sample_cells}"))
        set.seed(0)
        seurat_obj <- subset(seurat_obj, cells = sample(colnames(seurat_obj), down_sample_cells))
    }

    expr_m <- GetAssayData(seurat_obj, slot = "data")

    # Convert model to formula.
    model <- as.formula(model_design)

    # Relabel clustering. 
    test_vars <- meta %>%
        filter(merge_barcode %in% colnames(seurat_obj)) %>%
        mutate(
            enrichment_cluster = fct_other(as.character(liger_clusters), keep = as.character(test_cluster), other_level = "all_other_cells") %>% fct_relevel("all_other_cells")
        )

    run_lmer(expr_m, test_vars, model, cores)
}

run_lmer <- function(expr_m, test_vars, model, cores = NULL) {
    if (!is.null(cores)) {
        plan(multicore, workers = cores)
    }
    lm_out_obj_l <- future_apply(X = expr_m, MARGIN = 1, FUN = function(expr_r) {
        dat <- data.frame(expression = as.vector(expr_r), test_vars)
        # use tryCatch to return NA when model can't be fit for a gene
        tryCatch({
            broom::tidy(lmerTest::lmer(model, data = dat))
        },
        error = function(x) { print(x); return(NA) })
    })
    return(lm_out_obj_l)
}

# proportion detected across all genes.
prop_detected_generate <- function(meta, subcluster_so) {
    writeLines(str_glue("prop_detected: {paste0(unique(meta$ct_subcluster), collapse = ' ')}, {paste0(dim(subcluster_so), collapse = ' ')}"))

    expr_m <- GetAssayData(subcluster_so, slot = "data")

    expr_cells_gz <- rowSums(expr_m > 0)
    frac_gz <- as.double(expr_cells_gz) / ncol(expr_m)
    prop_detected_all <- tibble(gene = rownames(expr_m), prop_detected_all = frac_gz)

    return(prop_detected_all)
}

format_lm_long_output <- function(
    lm_out_obj_l,
    liger_meta,
    model_design,
    beta_regex){

    lm_tb <- lm_out_obj_l[!is.na(lm_out_obj_l)] %>%
        bind_rows(.id = "gene") %>%
        filter(grepl(beta_regex, term) | is.na(term))
    lm_tb <- lm_tb %>%
        mutate(
            p.value.adj = p.adjust(p.value, method = "BH"),
            model = model_design,
            region = paste0(unique(liger_meta$region)),
            cell_type = paste0(unique(liger_meta$cluster_cell_type))
        )
    return(lm_tb)
}

format_lm_output <- function(
    lm_out_obj_l,
    liger_meta,
    model_design,
    beta_regex){

    # format lm output into tibble
    lm_tb <- lm_out_obj_l[!is.na(lm_out_obj_l)] %>%
        bind_rows(.id = "gene") %>%
        filter(grepl(beta_regex, term) | is.na(term))
    
    if (ncol(lm_tb) != 9) {
        browser()
    }

    if (nrow(lm_tb) == 0) {
        stop("lm output is of length zero. Check lm_out_obj_l and beta_regex arguments")
    }

    # pivot dx output.
    lm_wider_spec <- build_wider_spec(lm_tb,
        names_from = "term",
        values_from = c("estimate", "p.value", "std.error", "statistic", "df"),
        names_glue = "{term}.{.value}"
    ) %>% arrange(term)
    lm_wd_tb <- pivot_wider_spec(lm_tb, id_cols = "gene", lm_wider_spec)

    # fdr_pvalue correct
    lm_filter_fdr <- lm_wd_tb %>%
        mutate(across(contains("p.value"), ~p.adjust(.x, method = "BH"), .names = "{.col}.adj"))

    # print pval counts after correction
    lm_filter_fdr %>%
        summarize(across(contains("p.value"), ~sum(.x < 0.05))) %>%
        print(width = Inf)

    # formating
    lm_filter_out <- lm_filter_fdr %>%
        dplyr::select(
            gene,
            contains("estimate"),
            contains("p.value"),
            contains("statistic"),
            contains("std.error"),
            contains("p.value.fdr")
        ) %>%
        mutate(
            model = model_design,
            region = paste0(unique(liger_meta$region)),
            cell_type = paste0(unique(liger_meta$cluster_cell_type))
        )

    return(lm_filter_out)
}

fmt_lm_output_noterm <- function(
    lm_out_obj_l,
    liger_meta,
    model_design,
    beta_regex
) { 
    lm_tb <- lm_out_obj_l[!is.na(lm_out_obj_l)] %>%
        bind_rows(.id = "gene") %>%
        filter(grepl(beta_regex, term) | is.na(term))
    
    lm_wider_spec <- build_wider_spec(lm_tb,
        names_from = "term",
        values_from = c("estimate", "p.value", "std.error", "statistic", "df"),
        names_glue = "{.value}"
    ) %>% arrange(term)
    
    lm_wd_tb <- pivot_wider_spec(lm_tb, id_cols = "gene", lm_wider_spec)

    lm_filter_fdr <- lm_wd_tb %>%
        mutate(across(contains("p.value"), ~p.adjust(.x, method = "BH"), .names = "{.col}.adj"))
    
    lm_filter_out <- lm_filter_fdr %>%
        dplyr::select(
            gene,
            contains("estimate"),
            contains("p.value"),
            contains("statistic"),
            contains("std.error"),
            contains("p.value.fdr")
        ) %>%
        mutate(
            model = model_design,
            region = paste0(unique(liger_meta$region)),
            cell_type = paste0(unique(liger_meta$cluster_cell_type))
        )
    return(lm_filter_out)
}

summarize_gene_counts <- function(broom_join, liger_clusters) {
    estimate_col <- str_glue("enrichment_cluster{liger_clusters}.estimate")
    fdr_col <- str_glue("enrichment_cluster{liger_clusters}.p.value.adj")
    up_ct <- broom_join %>%
        filter(.data[[estimate_col]] > 0, .data[[fdr_col]] < 0.1) %>%
        summarize(n = n()) %>%
        pluck("n")
    dn_ct <- broom_join %>%
        filter(.data[[estimate_col]] < 0, .data[[fdr_col]] < 0.1) %>%
        summarize(n = n()) %>%
        pluck("n")
    return(list(up_genes = up_ct, down_genes = dn_ct))
}

if (!interactive()) main()
