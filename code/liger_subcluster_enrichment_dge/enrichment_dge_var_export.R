# export top 100 variable genes per subcluster.

set.seed(0)
liblist <- c("Seurat", "tidyverse", "batchtools", "readxl", "ComplexHeatmap")
l <- lapply(liblist, function(x) suppressPackageStartupMessages(require(x, character.only = TRUE, quietly = TRUE)))
in_cluster_wk <- "../../analysis/seurat_lchen/liger_subcluster_enrichment_dge/subcluster_wk.rds"
clusters_exclude_file <- "../../resources/subclusters_removed_byQC_final.xlsx"
batchtools <- "../../analysis/seurat_lchen/liger_subcluster_enrichment_dge/batchtools"
out_path_base <- "../../analysis/seurat_lchen/liger_subcluster_enrichment_dge/enrichment_tables_var"
out_path_base_f <- "../../analysis/seurat_lchen/liger_subcluster_enrichment_dge/enrichment_tables_var_filtered"

main <- function() {
    dir.create(out_path_base)
    dir.create(out_path_base_f)
    cluster_wk_in <- readRDS(in_cluster_wk)
    excludes <- read_xlsx(clusters_exclude_file)
    reg <- loadRegistry(batchtools, conf.file = "../batchtools.conf.R", writeable = T)
    subcluster_wk <- cluster_wk_in %>%
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
    
    cluster_beta_var <- subcluster_wk %>%
        select(-liger_clusters, -region) %>%
        unnest(broom_join_fmt) %>%
        group_by(cluster_cell_type) %>%
        group_nest %>%
        mutate(top_var_genes = map(data, top_var_genes)) %>%
        select(-data)

    subcluster_cluster_wk <- subcluster_wk %>%
        left_join(cluster_beta_var, by = "cluster_cell_type") %>%
        mutate(region = fmt_display_region(region), cell_type = fmt_display_celltype(cluster_cell_type)) %>%
        mutate(gene_list = pmap(., function(...) {
            cr <- list(...)
            gene_list <- cr$broom_join_fmt %>%
                filter(p.value.adj < 0.05 & estimate > 0.1) %>%
                slice_head(n = 100) %>%
                left_join(cr$top_var_genes, by = "gene") %>%
                pluck("gene")
            cr$broom_join_fmt %>%
                filter(gene %in% gene_list) %>%
                mutate(ct_subcluster = str_glue("{cr$region}-{cr$cluster_cell_type}-{liger_clusters}"))
        }))

    subcluster_cluster_wk %>%
        group_by(cluster_cell_type, region) %>%
        group_nest %>%
        pmap(., function(...) {
            cr <- list(...)
            out_file <- file.path(out_path_base, str_glue("{cr$region}-{cr$cluster_cell_type}.csv"))
            writeLines(out_file)
            cr$data$gene_list %>%
                bind_rows %>%
                mutate(ct_subcluster = str_glue("{cr$region}-{cr$cluster_cell_type}-{liger_clusters}")) %>%
                group_by(liger_clusters) %>%
                arrange(ct_subcluster, desc(estimate), desc(is_signif)) %>%
                write_csv(out_file)
        })
    
    subcluster_cluster_wk %>%
        filter(!ct_subcluster %in% excludes$ct_subcluster) %>%
        group_by(cluster_cell_type, region) %>%
        group_nest %>%
        pmap(., function(...) {
            cr <- list(...)
            out_file <- file.path(out_path_base_f, str_glue("{cr$region}-{cr$cluster_cell_type}.csv"))
            writeLines(out_file)
            cr$data$gene_list %>%
                bind_rows %>%
                mutate(ct_subcluster = str_glue("{cr$region}-{cr$cluster_cell_type}-{liger_clusters}")) %>%
                group_by(liger_clusters) %>%
                arrange(ct_subcluster, desc(estimate), desc(is_signif)) %>%
                write_csv(out_file)
        })

    
}

fmt_display_region <- function(region) {
    str_replace_all(region, c("insula" = "INS", "preCG" = "BA4", "calcarine" = "V1"))
}

fmt_display_celltype <- function(ct) {
    str_replace_all(ct, c("opc" = "oligoprogenitor", "endothelia" = "endothelial"))
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

top_var_genes <- function(enrichment_dge_list) {
    enrichment_dge_list %>%
        group_by(gene) %>%
        filter(length(unique(region)) == length(unique(enrichment_dge_list$region))) %>%
        summarize(beta_var = var(estimate)) %>%
        arrange(desc(beta_var)) %>%
        slice_head(n = 100) 
}

if (!interactive()) main()
