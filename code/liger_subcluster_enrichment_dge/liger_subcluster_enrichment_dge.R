set.seed(0)
liblist <- c("Seurat", "tidyverse", "batchtools", "BiocParallel", "parallel", "future.apply", "lme4", "lmerTest", "broom.mixed")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)
options(future.globals.maxSize = Inf, deparse.max.lines = 5)

in_seurat_rds <-
  "../../analysis/pci_import/pci_seurat.rds"
in_seurat_liger <-
  "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"
in_subcluster_celltype_filter <- 
    "../../analysis/seurat_lchen/liger_subcluster_celltype_filter/celltype_filter.rds"

## Outputs
out_path_base <- "../../analysis/seurat_lchen/liger_subcluster_enrichment_dge/"
batchtools <- file.path(out_path_base, "batchtools")
dir.create(out_path_base, recursive = TRUE)

RESOURCES <- list(
    ncpus = 16,
    memory = 80,
    walltime = 172800,
    partition = "bigmem",
    measure.memory = TRUE
)
prop_detected_filter <- 0.1
chunk_size <- 30

main <- function() {
    if (dir.exists(batchtools)) {
        reg <- loadRegistry(batchtools, conf.file = "../batchtools.conf.R", writeable = TRUE)
        sweepRegistry()
        reg$packages <- liblist
    } else {
        reg <- makeRegistry(batchtools, packages = liblist, conf.file = "../batchtools.conf.R")
    }
    liger_meta <- readRDS(in_seurat_liger)
    subcluster_celltype_filter <- readRDS(in_subcluster_celltype_filter)

    # Filter metadata to get cell ids for matching region, cell_type.
    # Subset seurat object.
    liger_meta_subset <- liger_meta %>%
        mutate( 
            liger_clusters = fct_inseq(liger_clusters),
            ct_subcluster = paste(region, cluster_cell_type, liger_clusters, sep = "-"),
            log_number_umi = log(number_umi)
        ) %>%
        mutate(
            ct_subcluster = factor(ct_subcluster, levels = str_sort(unique(ct_subcluster), numeric = TRUE))
        ) %>%
        filter(cluster_cell_type == cell_type)

    # Add model design. MUST include enrichment_cluster as a fixed term.
    model_designs <- c("expression ~ enrichment_cluster + pmi + age + sex + number_umi + percent_mito + (1 | library_id)")

    cluster_tbs <- liger_meta_subset %>%
        inner_join(tibble(model_design = model_designs), by = character()) %>%
        arrange(ct_subcluster, model_design) %>%
        group_by(region, cluster_cell_type, model_design) %>%
        group_nest(keep = TRUE)

    subcluster_wk <- liger_meta_subset %>%
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
    saveRDS(subcluster_wk, file.path(out_path_base, "subcluster_wk.rds"), compress = FALSE)
    waitForJobs()

    subcluster_wk <- subcluster_wk %>%
        mutate(broom_join = pmap(list(data, model_design, job.id), function(data, model_design, job.id) {
                gc()
                tryCatch({
                    broom_list <- loadResult(job.id)
                    writeLines(paste(job.id))
                    format_lm_output(broom_list, data, model_design, "enrichment_cluster")
                }, 
                error = function(x) { print(x); return(NA) })
            }))

    
    saveRDS(subcluster_wk, file.path(out_path_base, "subcluster_wk.rds"))

    # Merge all subcluster results per region / celltype and write out as separate csv files. 
    subcluster_wk %>%
        filter(!is.na(broom_join)) %>%
        group_by(region, cluster_cell_type) %>%
        group_split %>%
        walk(function(x) {
            join_cols <- c("gene", "cell_type", "model", "region")
            filepath <- file.path(
                out_path_base,
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

    subcluster_wk <- subcluster_wk %>%
        filter(!is.na(broom_join_noterm)) %>%
        mutate(broom_join_fmt = pmap(., function(...) {
                cr <- list(...)
                cr$broom_join_noterm %>%
                    mutate(region = fmt_display_region(region), cell_type = fmt_display_celltype(cell_type), liger_clusters = cr$liger_clusters) %>%
                    mutate(is_signif = (p.value.adj < 0.05)) %>%
                    arrange(desc(is_signif), desc(abs(estimate)))
            }))

        
    subcluster_wk %>%
        pwalk(., function(...) {
            cr <- list(...)
            disp_ct_subcluster <- str_glue("{fmt_display_region(cr$region)}-{fmt_display_celltype(cr$cluster_cell_type)}-{cr$liger_clusters}")
            filepath <- file.path(
                out_path_base,
                "enrichment_tables_disp",
                str_glue("{disp_ct_subcluster}.csv")
            )
            writeLines(as.character(filepath))
            dir.create(dirname(filepath), showWarnings = F)
            write_csv(cr$broom_join_fmt, filepath)
        })
    
    bind_rows(subcluster_wk$broom_join_fmt) %>%
        saveRDS(file.path(out_path_base, "enrichment_tables_disp", "enrichment_dge.rds"))

    subcluster_gene_summary <- subcluster_wk %>%
        filter(!is.na(broom_join)) %>%
        mutate(gene_cts = pmap_df(., function(...) {
                cr <- list(...)
                summarize_gene_counts(cr$broom_join, cr$liger_clusters)
            })) %>%
        select(region, cluster_cell_type, ct_subcluster, liger_clusters, gene_cts) %>%
        unnest(gene_cts)
    write_csv(subcluster_gene_summary, file.path(out_path_base, "enrichment_dge_summary.csv"))
}

subcluster_worker <- function(test_cluster, meta, prop_detected_filter) {
    # Load in subcluster file.
    options(future.globals.maxSize = Inf)
    envdata <- new.env()
    load(file.path("..", unique(meta$filepath)), envdata)
    
    # Subset.
    subcluster_so <- subset(envdata$liger_so, cells = meta$cell_ids)
    rm(envdata)
    gc()

    print(dim(subcluster_so))
    # Calculate proportion detected.
    prop_detected_tbl <- prop_detected_generate(meta, subcluster_so)
    detected <- prop_detected_tbl %>% filter(prop_detected_all > prop_detected_filter)

    subcluster_so <- subset(subcluster_so, cells = meta$cell_ids, features = detected$gene)
    print(dim(subcluster_so))

    # Run lmer_enrichment_de
    return(run_lmer_enrichment_de(subcluster_so, unique(meta$model_design), test_cluster, down_sample_cells = 10000))
}

run_lmer_enrichment_de <- function(
    seurat_obj,
    model_design,
    test_cluster,
    down_sample_cells = NULL,
    cores = NULL){

    # Downsample if down_sample_cells defined.
    if (!is.null(down_sample_cells) && ncol(seurat_obj) > down_sample_cells) {
        writeLines(str_glue("downsampling to {down_sample_cells}"))
        set.seed(0)
        seurat_obj <- subset(seurat_obj, cells = sample(colnames(seurat_obj), down_sample_cells))
        seurat_obj@meta.data %>%
            select(region, clinical_dx, cluster_cell_type) %>% table %>% print
    }

    expr_m <- GetAssayData(seurat_obj, slot = "data")

    # Convert model to formula.
    model <- as.formula(model_design)
    test_vars <- seurat_obj@meta.data

    # Relabel clustering. 
    test_vars <- test_vars %>%
        mutate(
            enrichment_cluster = fct_other(liger_clusters, keep = test_cluster, other_level = "all_other_cells") %>% fct_relevel("all_other_cells")
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

# Rewrite region name to be consistent with display name.
# insula -> INS
# preCG -> BA4, Broadmann area 4 on precentral gyrus
# calcarine -> V1, visual cortex in calcarine fissure
fmt_display_region <- function(region) {
    str_replace_all(region, c("insula" = "INS", "preCG" = "BA4", "calcarine" = "V1"))
}

fmt_display_celltype <- function(ct) {
    str_replace_all(ct, c("opc" = "oligoprogenitor", "endothelia" = "endothelial"))
}

if (!interactive()) {
    main()
}
