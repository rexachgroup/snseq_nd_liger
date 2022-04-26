# Plot hla locus / cwow genelist for microglia.
liblist <- c("tidyverse", "patchwork")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

IN_CLUSTER_DGE_WK <- "../../analysis/seurat_lchen/liger_subcluster_lme/subcluster_wk.rds"
IN_CLUSTER_HIER <- "../../analysis/seurat_lchen/liger_subcluster_hier/hclust/dge_enrichment_var_beta_hclust.rds"
GENE_LIST <- c("../../resources/HLAlocusGenesNotes_mic.csv", "../../resources/HLAlocusGenesNotes_mic2.csv")
OUT_DIR <- "../../analysis/seurat_lchen/liger_subcluster_hier/genelist/"
plot_ts <- str_glue("{system('md5sum liger_subcluster_hla.R', intern = TRUE)} {date()}")

main <- function() {
    cluster_wk_in <- readRDS(IN_CLUSTER_DGE_WK)
    cluster_hier <- readRDS(IN_CLUSTER_HIER)

    imap(GENE_LIST, function(x, i) {
        print(x)
        genelist <- read_csv(x)
        
        hclust <- filter(cluster_hier, cluster_cell_type == "microglia") %>% pluck("hclust", 1)
        cluster_wk <- filter(cluster_wk_in, 
            cluster_cell_type == "microglia", ct_subcluster %in% hclust$labels)
        hclust_order <- hclust$labels[hclust$order]
        
        lme_tb <- fmt_cluster_dx_dge(cluster_wk, genelist)
        lme_tb <- lme_tb %>%
            mutate(
                gene = fct_relevel(gene, rev(genelist$gene)), 
                stars_label = stars_label(z),
                ct_subcluster = fct_relevel(ct_subcluster, hclust_order)
            )

        yheight <- (length(unique(lme_tb$gene)) * 0.25) + 3

        pdf(file.path(OUT_DIR, str_glue("hlalocus_{i}.pdf")), height = yheight, width = 8)
        print(plot_split_dx(lme_tb))
        graphics.off()
        
    })

}

stars_label <- function(vec) {
    as.character(symnum(vec,
        cutpoints = c(-Inf, -2, 2, Inf),
        symbols = c("-", " ", "+")))
}

plot_split_dx <- function(lme_tb) {
   lme_tb %>%
    group_split(dx) %>%
    map(
        ~ggplot(., aes(x = ct_subcluster, y = gene, fill = beta)) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
            geom_tile() +
            geom_text(aes(label = stars_label), show.legend = TRUE) +
            scale_fill_gradient2(
                low = scales::muted("blue"),
                mid = "white",
                high = scales::muted("red"),
                limits = c(-0.25, 0.25),
                oob = scales::squish,
            ) +
        labs(
            title = str_glue("LME DGE in {unique(.$dx)} vs. Control"),
            subtitle = str_glue(
                plot_ts,
                "\n",
                "beta -> color, ",
                "- -> t_statistic < -2, ",
                "+ -> t_statistic > 2"
            )
        )
    )
}

fmt_cluster_dx_dge <- function(cluster_dx_wk, markers) {
    dx <- c("AD", "bvFTD", "PSP-S")
    clinical_dx_beta_cols <- str_glue("clinical_dx{dx}.estimate")
    clinical_dx_pval_cols <- str_glue("clinical_dx{dx}.p.value")
    clinical_dx_fdr_cols <- str_glue("clinical_dx{dx}.p.value.adj")
    clinical_dx_z_cols <- str_glue("clinical_dx{dx}.statistic")
    filter_args <- tibble(dx, clinical_dx_beta_cols, clinical_dx_pval_cols, clinical_dx_fdr_cols, clinical_dx_z_cols)
   
    filtered_gene_tests <- cluster_dx_wk %>%
        mutate(ct_subcluster = fct_drop(ct_subcluster)) %>%
        full_join(., filter_args, by = character())

    filter_call <- function(...) {
        cr <- list(...)
        beta_col <- cr$clinical_dx_beta_cols
        pval_col <- cr$clinical_dx_pval_cols
        fdr_col <- cr$clinical_dx_fdr_cols
        z_col <- cr$clinical_dx_z_cols
        writeLines(str_glue("{cr$ct_cluster}:  {beta_col} {pval_col} {fdr_col} {z_col}"))
        if (!is.na(cr$broom_join) && nrow(cr$broom_join) > 1 && all(c(beta_col, fdr_col, z_col) %in% colnames(cr$broom_join))) {
            broom_filter <- cr$broom_join %>%
                dplyr::filter(gene %in% markers$gene) %>%
                rename(beta = .data[[beta_col]], pval = .data[[pval_col]], fdr = .data[[fdr_col]], z = .data[[z_col]]) %>%
                select(gene, beta, pval, fdr, z)
        }
    }
    filtered_gene_tb <- filtered_gene_tests %>%
        mutate(lme_marker_estimates = pmap(., filter_call)) %>%
        select(region, dx, ct_subcluster, lme_marker_estimates) %>%
        unnest(lme_marker_estimates) %>%
        complete(ct_subcluster, dx)
    return(filtered_gene_tb) 
}

if (!interactive()) main()
