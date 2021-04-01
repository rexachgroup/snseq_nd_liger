set.seed(0)
liblist <- c("Seurat", "tidyverse", "eulerr", "grid", "gridExtra")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)
options(future.globals.maxSize = Inf, deparse.max.lines = 5)

in_cluster_wk <- "../../analysis/seurat_lchen/seurat_cluster_lme/cluster_wk.rds"
out_path_base <- "../../analysis/seurat_lchen/seurat_cluster_lme"

beta_threshold <- 0.1
fdr_threshold <- 0.5

main <- function() {
    writeLines("load meta")
    cluster_wk <- readRDS(in_cluster_wk) %>%
        select(region, ct_cluster, job_id, broom_join)
    
    writeLines("filter each broom result to cluster_up and cluster_down")
    dx <- c("AD", "bvFTD", "PSP-S")
    fill_colors = c("#F8766D", "#00BA38", "#619CFF")
    clinical_dx_beta_cols <- str_glue("clinical_dx{dx}.estimate")
    clinical_dx_fdr_cols <- str_glue("clinical_dx{dx}.p.value.adj")
    beta_direction <- c("up", "down")
    filter_args <- full_join(tibble(dx, fill_colors, clinical_dx_beta_cols, clinical_dx_fdr_cols), tibble(beta_direction), by = character())
    glimpse(filter_args)

    filtered_gene_tests <- full_join(cluster_wk, filter_args, by = character())
    glimpse(filtered_gene_tests)

    # Filter each broom_join according to the filter_args parameter for that row in filtered_gene_tests.
    # Has to be defined separately from the mutate/pmap below: https://github.com/tidyverse/dplyr/issues/4575#issuecomment-537347938
    filter_call <- function(...) {
        cr <- list(...)
        beta_col <- cr$clinical_dx_beta_cols
        fdr_col <- cr$clinical_dx_fdr_cols
        writeLines(str_glue("{cr$ct_cluster}:  {cr$clinical_dx_beta_cols} {cr$clinical_dx_fdr_cols} {cr$beta_direction}"))
        if (all(c(beta_col, fdr_col) %in% colnames(cr$broom_join))) {
            if (cr$beta_direction == "up") {
                filtered_genes <- cr$broom_join %>%
                    filter(
                        .data[[beta_col]] > beta_threshold,
                        .data[[fdr_col]] < fdr_threshold
                    ) %>%
                    pluck("gene")
            } else if (cr$beta_direction == "down") {
                filtered_genes <- cr$broom_join %>%
                    filter(
                        .data[[beta_col]] < -beta_threshold,
                        .data[[fdr_col]] < fdr_threshold
                    ) %>%
                    pluck("gene")
            }
            writeLines(str_glue("{length(filtered_genes)}"))
        } else {
            writeLines(str_glue("missing col: {beta_col}, {fdr_col}"))
            filtered_genes <- c()
        }
        return(filtered_genes)
    }

    filtered_gene_list <- filtered_gene_tests %>%
        mutate(lme_signif_genes = pmap(., filter_call))

    euler_diags <- filtered_gene_list %>%
        group_by(region, ct_cluster, beta_direction) %>%
        group_nest() %>%
        mutate(euler_objs = pmap(., function(...) {
            cr <- list(...)
            gene_lists <- setNames(cr$data$lme_signif_genes, cr$data$dx)
            return(euler(gene_lists))
        }))
    
    pdf(file.path(out_path_base, "celltype_dge_overlaps.pdf"), width = 4, height = 4)
    pwalk(euler_diags, function(...) {
        cr <- list(...)
        print(plot(
            cr$euler_objs, 
            fills = cr$data$fill_colors,
            quantities = list(fontsize = 16),
            labels = list(fontsize = 16), 
            main = list(str_glue("{cr$ct_cluster}_{cr$beta_direction}", fontsize = 8))
        ))
    })
    graphics.off()
}

if (!interactive()) {
    main()
}
