set.seed(0)
options(deparse.max.lines = 5)
liblist <- c("Seurat", "tidyverse", "readxl", "batchtools", "lme4", "lmerTest", "future.apply", "broom.mixed")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

cellbender_tb <- "../../analysis/seurat_lchen/cellbender/merge/cellbender_region.rds"
in_seurat_meta <- "../../analysis/pci_import/pci_seurat_meta.rds"
OUT_DIR <- "../../analysis/seurat_lchen/cellbender/dge_no_logumi/"
batchtools <- file.path(OUT_DIR, "batchtools")

RESOURCES <- list(
    ncpus = 1,
    memory = 64,
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
    } else {
        reg <- makeRegistry(batchtools, packages = liblist, conf.file = "../batchtools.conf.R")
    }
    reg$packages <- liblist

    meta <- readRDS(in_seurat_meta)
    cellbender_region <- readRDS(cellbender_tb) %>%
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
    
    # Split by cluster and add model design and seurat obj.
    cluster_wk <- meta_subset %>%
        inner_join(tibble(model_design = model_designs), by = character()) %>%
        inner_join(cellbender_region, by = c("region")) %>%
        arrange(ct_cluster, model_design) %>%
        group_by(region, cluster_cell_type, ct_cluster, model_design) %>%
        group_nest(keep = TRUE)
    
    clearRegistry()
    batchExport(list(run_lmer_de = run_lmer_de, run_lmer = run_lmer,
                     prop_detected_generate = prop_detected_generate, prop_cells_detected_dx = prop_cells_detected_dx), reg = reg)
    ids <- batchMap(cluster_worker,
        args = list(meta = cluster_wk$data),
        more.args = list(prop_detected_filter = prop_detected_filter)
    )
    ids <- findNotDone() %>%
        mutate(chunk = chunk(job.id, chunk.size = chunk_size))
    cluster_wk$job_id <- getJobTable()$job.id
    cluster_wk$prop_detected_filter <- prop_detected_filter
    submitJobs(ids, resources = RESOURCES)
    saveRDS(cluster_wk, file.path(OUT_DIR, "cluster_wk.rds"))
    waitForJobs()

    cluster_tb <- cluster_wk %>%
        filter(job_id %in% findDone()$job.id) %>%
        mutate(broom_join = pmap(list(data, model_design, job_id), function(data, model_design, job_id) {
                gc()
                broom_list <- loadResult(job_id)
                writeLines(paste(job_id))
                if (all(is.na(broom_list))) {
                    return(NA)
                }
                format_lm_output(broom_list, data, model_design, "clinical_dx")
            }))

    saveRDS(cluster_tb, file.path(OUT_DIR, "cluster_wk.rds"), compress = FALSE)

    cluster_tbs <- cluster_tb %>%
        filter(map_lgl(broom_join, ~!any(is.na(.)))) %>%
        pluck("broom_join") %>%
        bind_rows

    saveRDS(cluster_tbs, file.path(OUT_DIR, "cluster_lme.rds"), compress = F)
    
    cluster_tb <- cluster_tb %>%
        mutate(filepath = file.path(
            OUT_DIR, 
            "lme_tables",
            paste0(ct_cluster, ".csv"))
        )

    pwalk(list(cluster_tb$filepath, cluster_tb$broom_join), function(filepath, tb) {
            dir.create(dirname(filepath), recursive = TRUE, showWarnings = FALSE)
            if (!all(is.na(tb))) {
                write_csv(tb, filepath)
            }
        })
    
    cluster_gene_summary <- cluster_tb %>%
        filter(!is.na(broom_join)) %>%
        mutate(gene_cts = pmap(., function(...) {
                cr <- list(...)
                summarize_gene_counts(cr$broom_join)
            })) %>%
        select(region, cluster_cell_type, ct_cluster, gene_cts)

    cluster_gene_summary %>%
        unnest(gene_cts, names_repair = "universal") %>%
        select(region, cluster_cell_type, ct_cluster = ct_cluster...3, dx, up, down) %>%
        write_csv(file.path(OUT_DIR, "dx_dge_summary.csv"))
}

cluster_worker <- function(meta, prop_detected_filter) {
    options(future.globals.maxSize = Inf)
    sobj <- readRDS(unique(meta$merge_out))
    common_bc <- intersect(colnames(sobj), meta$merge_barcode)
    cluster_so <- subset(sobj, cells = common_bc)
    meta_s <- meta[match(common_bc, meta$merge_barcode),]
    rm(sobj)
    gc()

    # Calculate proportion detected.
    prop_detected_tbl <- prop_detected_generate(meta_s, cluster_so)
    detected <- prop_detected_tbl %>% filter(prop_detected_AD > prop_detected_filter |
                      prop_detected_Control > prop_detected_filter |
                      prop_detected_bvFTD > prop_detected_filter |
                      `prop_detected_PSP-S` > prop_detected_filter)
    cluster_so <- subset(cluster_so, features = detected$gene)
    return(run_lmer_de(cluster_so, meta_s, down_sample_cells = 10000, cores = RESOURCES$ncores))
}

run_lmer_de <- function(
    seurat_obj,
    meta,
    down_sample_cells = NULL,
    cores = NULL){

    # Downsample if down_sample_cells defined.
    if (!is.null(down_sample_cells) && ncol(seurat_obj) > down_sample_cells) {
        writeLines(str_glue("downsampling to {down_sample_cells}"))
        use_cells <- sample(colnames(seurat_obj), down_sample_cells)
    } else {
        use_cells <- colnames(seurat_obj)
    }
    common_bc <- intersect(use_cells, meta$merge_barcode)
    meta_s <- meta[match(common_bc, meta$merge_barcode),]
    seurat_obj <- subset(seurat_obj, cells = common_bc)
    model <- as.formula(meta_s$model_design)

    expr_m <- GetAssayData(seurat_obj, slot = "data")

    run_lmer(expr_m, meta_s, model, cores)
}

run_lmer <- function(expr_m, test_vars, model, cores = NULL) {
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


prop_cells_detected_dx <- function(expr, meta, clinical_dx_val) {
    is_dx <- meta[["clinical_dx"]] == clinical_dx_val
    expr_genes_dx <- expr[,is_dx]
    expr_isexpressed <- rowSums(expr > 0)
    prop_detected <- expr_isexpressed / sum(is_dx)
    prop_detected_tb <- enframe(prop_detected, name = "gene", value = paste0("prop_detected_", clinical_dx_val))
    return(prop_detected_tb)
}

prop_detected_generate <- function(meta, subcluster_so) {
    #writeLines(str_glue("prop_detected: {unique(meta$ct_subcluster)}, {paste0(dim(subcluster_so), collapse = ' ')}"))

    expr <- GetAssayData(subcluster_so, slot = "data")

    # proportion of control cells expressing gene
    ctl_detect <- prop_cells_detected_dx(expr, meta, "Control")

    # proportion of dx cells expressing gene
    ad_detect <- prop_cells_detected_dx(expr, meta, "AD")
    ftd_detect <- prop_cells_detected_dx(expr, meta, "bvFTD")
    psp_detect <- prop_cells_detected_dx(expr, meta, "PSP-S")

    return(reduce(list(ad_detect, ftd_detect, psp_detect),
           inner_join,
           by = "gene",
           .init = ctl_detect))
}

format_lm_output <- function(
    lm_out_obj_l,
    meta,
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

    # formatting
    lm_filter_out <- lm_filter_fdr %>%
        dplyr::select(
            gene,
            contains("estimate"),
            contains("p.value"),
            contains("statistic"),
            contains("prop_detected"),
            contains("fdr.p.value")
        ) %>%
        mutate(
            model = unique(model_design),
            region = paste0(unique(meta$region)),
            cell_type = paste0(unique(meta$cluster_cell_type)),
            ct_cluster = paste0(unique(meta$ct_cluster))
        )

    return(lm_filter_out)
}

summarize_gene_counts <- function(broom_join) {
    dx <- c("AD", "bvFTD", "PSP-S")
    dx_map <- map(dx, function(clinical_dx) {
        estimate_col <- str_glue("clinical_dx{clinical_dx}.estimate")
        fdr_col <- str_glue("clinical_dx{clinical_dx}.p.value.adj")
        if (estimate_col %in% colnames(broom_join) && fdr_col %in% colnames(broom_join)) {
            up_ct <- broom_join %>%
                group_by(ct_cluster) %>%
                filter(.data[[estimate_col]] > 0, .data[[fdr_col]] < 0.1) %>%
                summarize(up = n())
            dn_ct <- broom_join %>%
                group_by(ct_cluster) %>%
                filter(.data[[estimate_col]] < 0, .data[[fdr_col]] < 0.1) %>%
                summarize(down = n())
            return(tibble(up_genes = up_ct, down_genes = dn_ct))
        }
    }) %>% setNames(nm = dx) %>%
    bind_rows(.id = "dx") %>%
    unnest(c(up_genes, down_genes), names_repair = "universal") %>%
    select(dx, ct_cluster = ct_cluster...2, up, down)
    return(dx_map)
}

if (!interactive()) main()
