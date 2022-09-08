set.seed(0)
liblist <- c("Seurat", "tidyverse", "variancePartition", "batchtools", "edgeR", "BiocParallel", "future.apply", "lme4", "lmerTest", "broom.mixed")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)
options(future.globals.maxSize = Inf)

in_seurat_rds <-
  "../../analysis/pci_import/pci_seurat.rds"
in_seurat_liger <-
  "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"
in_subcluster_celltype_filter <- 
    "../../analysis/seurat_lchen/liger_subcluster_celltype_filter/celltype_filter.rds"

## Outputs
out_path_base <- "../../analysis/seurat_lchen/liger_subcluster_lme/"
batchtools <- file.path(out_path_base, "batchtools")
dir.create(out_path_base, recursive = TRUE)

RESOURCES <- list(
    ncpus = 1,
    memory = 80,
    walltime = 172800,
    measure.memory = TRUE,
    nice = 100,
    chunks.as.arrayjobs = TRUE,
    partition = "bigmem"
)
prop_detected_filter <- 0.1
chunk_size <- 40

main <- function() {
    if (dir.exists(batchtools)) {
        reg <- loadRegistry(batchtools, conf.file = "../batchtools.conf.R", writeable = TRUE)
        sweepRegistry()
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
        )
 

    # Split by subcluster + add model design. 
    model_designs <- c("expression ~ clinical_dx + pmi + age + sex + number_umi + percent_mito + (1 | library_id)")

    subcluster_wk <- liger_meta_subset %>%
        inner_join(tibble(model_design = model_designs), by = character()) %>%
        arrange(ct_subcluster, model_design) %>%
        group_by(region, cluster_cell_type, ct_subcluster, model_design) %>%
        group_nest(keep = TRUE)

    clearRegistry()
    batchExport(list(run_lmer_de = run_lmer_de, run_lmer = run_lmer,
                     prop_detected_generate = prop_detected_generate, prop_cells_detected_dx = prop_cells_detected_dx), reg = reg)
    ids <- batchMap(subcluster_worker,
        args = list(meta = subcluster_wk$data),
        more.args = list(prop_detected_filter = prop_detected_filter)
    )
    ids <- findNotDone() %>%
        mutate(chunk = chunk(job.id, chunk.size = chunk_size))
    subcluster_wk$job_id <- getJobTable()$job.id
    subcluster_wk$prop_detected_filter <- prop_detected_filter
    submitJobs(ids, resources = RESOURCES)
    saveRDS(subcluster_wk, file.path(out_path_base, "subcluster_wk.rds"))
    waitForJobs()

    # Fetch results.
    subcluster_tb <- subcluster_wk %>%
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
    
    saveRDS(subcluster_tb, file.path(out_path_base, "subcluster_wk.rds"), compress = FALSE)

    # Merge all subcluster results into one large table and save.
    subcluster_tbs <- bind_rows(subcluster_tb$broom_join[map_lgl(subcluster_tb$broom_join, ~!any(is.na(.))) ])

    saveRDS(subcluster_tbs, file.path(out_path_base, "subcluster_lme.rds"), compress = FALSE)

    # For csv export, write each subcluster to a separate file.
    subcluster_tb <- subcluster_tb %>%
        mutate(filepath = file.path(
            out_path_base, 
            "lme_tables",
            paste(region, cluster_cell_type, sep = "-"),
            paste0(ct_subcluster, ".csv"))
        )

    pwalk(list(subcluster_tb$filepath, subcluster_tb$broom_join), function(filepath, tb) {
            dir.create(dirname(filepath), recursive = TRUE, showWarnings = FALSE)
            if (!all(is.na(tb))) {
                write_csv(tb, filepath)
            }
        })
    
    subcluster_gene_summary <- subcluster_tb %>%
        filter(!is.na(broom_join)) %>%
        mutate(gene_cts = pmap(., function(...) {
                cr <- list(...)
                summarize_gene_counts(cr$broom_join)
            })) %>%
        select(region, cluster_cell_type, ct_subcluster, gene_cts)

    subcluster_gene_summary %>%
        unnest(gene_cts, names_repair = "universal") %>%
        select(region, cluster_cell_type, ct_subcluster = ct_subcluster...3, dx, up, down) %>%
        write_csv(file.path(out_path_base, "dx_dge_summary.csv"))

}

subcluster_worker <- function(meta, prop_detected_filter) {
    options(future.globals.maxSize = Inf)
    envdata <- new.env()
    load(file.path("..", unique(meta$filepath)), envdata)
    
    # Subset.
    subcluster_so <- subset(envdata$liger_so, cells = meta$cell_ids)
    rm(envdata)
    gc()

    # Calculate proportion detected.
    prop_detected_tbl <- prop_detected_generate(meta, subcluster_so)
    detected <- prop_detected_tbl %>% filter(prop_detected_AD > prop_detected_filter |
                      prop_detected_Control > prop_detected_filter |
                      prop_detected_bvFTD > prop_detected_filter |
                      `prop_detected_PSP-S` > prop_detected_filter)
    subcluster_so <- subset(subcluster_so, cells = meta$cell_ids, features = detected$gene)

    return(run_lmer_de(subcluster_so, unique(meta$model_design), down_sample_cells = 10000, cores = RESOURCES$ncores))
}

run_lmer_de <- function(
    seurat_obj,
    model_design,
    down_sample_cells = NULL,
    cores = NULL){

    # Downsample if down_sample_cells defined.
    if (!is.null(down_sample_cells) && ncol(seurat_obj) > down_sample_cells) {
        writeLines(str_glue("downsampling to {down_sample_cells}"))
        seurat_obj <- subset(seurat_obj, cells = sample(colnames(seurat_obj), down_sample_cells))
        seurat_obj@meta.data %>%
            select(region, clinical_dx, cluster_cell_type) %>% table %>% print
    }

    expr_m <- GetAssayData(seurat_obj, slot = "data")

    # Convert model to formula.
    # If dx is present relevel seurat_obj so that control is the 1st factor level.
    model <- as.formula(model_design)
    test_vars <- seurat_obj@meta.data

    if ("clinical_dx" %in% colnames(test_vars)) {
        test_vars <- mutate(test_vars, clinical_dx = fct_relevel(clinical_dx, "Control"))
    }

    run_lmer(expr_m, test_vars, model, cores)
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
    prop_detected_dx <- apply(expr, 1, function(row) {
        row_dx <- row[is_dx]
        prop_detected <- (sum(row_dx > 0) / length(row_dx))
        prop_detected[!is.finite(prop_detected)] <- 0
        return(prop_detected)
    }) %>%
    round(3) %>%
    enframe(name = "gene", value = paste0("prop_detected_", clinical_dx_val))
    return(prop_detected_dx)
}

prop_detected_generate <- function(meta, subcluster_so) {
    writeLines(str_glue("prop_detected: {unique(meta$ct_subcluster)}, {paste0(dim(subcluster_so), collapse = ' ')}"))

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
            contains("prop_detected"),
            contains("fdr.p.value")
        ) %>%
        mutate(
            model = model_design,
            region = paste0(unique(liger_meta$region)),
            cell_type = paste0(unique(liger_meta$cluster_cell_type)),
            ct_subcluster = paste0(unique(liger_meta$ct_subcluster))
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
                group_by(ct_subcluster) %>%
                filter(.data[[estimate_col]] > 0, .data[[fdr_col]] < 0.1) %>%
                summarize(up = n())
            dn_ct <- broom_join %>%
                group_by(ct_subcluster) %>%
                filter(.data[[estimate_col]] < 0, .data[[fdr_col]] < 0.1) %>%
                summarize(down = n())
            return(tibble(up_genes = up_ct, down_genes = dn_ct))
        }
    }) %>% setNames(nm = dx) %>%
    bind_rows(.id = "dx") %>%
    unnest(c(up_genes, down_genes), names_repair = "universal") %>%
    select(dx, ct_subcluster = ct_subcluster...2, up, down)
    return(dx_map)
}

if (!interactive()) {
    main()
}
