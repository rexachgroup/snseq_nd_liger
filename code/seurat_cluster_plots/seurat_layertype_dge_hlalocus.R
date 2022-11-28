# Plot per-gene DGE for cwow? grant genelist.
liblist <- c("Seurat", "tidyverse", "scales", "readxl", "ComplexHeatmap")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)
options(future.globals.maxSize = Inf, deparse.max.lines = 5)

in_cluster_wk <- "../../analysis/seurat_lchen/seurat_cluster_lme/cluster_wk.rds"
genelist <- "../../resources/HLAlocusGenesNotes.csv"
out_path_base <- "../../analysis/seurat_lchen/seurat_cluster_lme"

args <- tibble(
    genelist = file.path("../../resources/", c("HLAlocusGenesNotes.csv", "HLAlocusGenesNotes2.csv")),
    out_path = file.path(out_path_base, c("hlalocus_sub1.pdf", "hlalocus_sub2.pdf")),
    limits = list(c(-0.3, 0.3), c(-0.3, 0.3))
)

main <- function() {
    writeLines("load cluster wk")
    cluster_wk <- readRDS(in_cluster_wk)

    pwalk(args, function(...) {
        cr <- list(...)
        marker_gene_tb <- read_csv(cr$genelist)
        plot_worker(cluster_wk, marker_gene_tb, cr$out_path, cr$limits)
    })
}

plot_worker <- function(cluster_wk, marker_gene_tb, filepath, limits) {
    dx <- c("AD", "bvFTD", "PSP-S")
    clinical_dx_beta_cols <- str_glue("clinical_dx{dx}.estimate")
    clinical_dx_pval_cols <- str_glue("clinical_dx{dx}.p.value")
    clinical_dx_fdr_cols <- str_glue("clinical_dx{dx}.p.value.adj")
    clinical_dx_z_cols <- str_glue("clinical_dx{dx}.statistic")
    filter_args <- tibble(dx, clinical_dx_beta_cols, clinical_dx_pval_cols, clinical_dx_fdr_cols, clinical_dx_z_cols)

    filtered_gene_tests <- full_join(cluster_wk, filter_args, by = character())
    glimpse(filtered_gene_tests)
    filter_call <- function(...) {
        cr <- list(...)
        beta_col <- cr$clinical_dx_beta_cols
        pval_col <- cr$clinical_dx_pval_cols
        fdr_col <- cr$clinical_dx_fdr_cols
        z_col <- cr$clinical_dx_z_cols
        writeLines(str_glue("{cr$ct_cluster}:  {beta_col} {pval_col} {fdr_col} {z_col}"))
        if (!is.na(cr$broom_join) && nrow(cr$broom_join) > 1 && all(c(beta_col, fdr_col, z_col) %in% colnames(cr$broom_join))) {
            broom_filter <- cr$broom_join %>%
                dplyr::filter(gene %in% marker_gene_tb$gene) %>%
                rename(beta = .data[[beta_col]], pval = .data[[pval_col]], fdr = .data[[fdr_col]], z = .data[[z_col]]) %>%
                select(gene, beta, pval, fdr, z)
        }
    }

    writeLines("gene filter")
    filtered_gene_list <- filtered_gene_tests %>%
        filter(!cluster_cell_type %in% c("ependymal", "pericyte")) %>%
        mutate(cluster_cell_type = fct_recode(cluster_cell_type, lymphocytes = "t_cell"), ct_cluster = paste(region, cluster_cell_type, sep = "-")) %>%
        mutate(lme_marker_estimates = pmap(., filter_call)) %>%
        glimpse
    
    plot_gene_list <- filtered_gene_list %>%
        select(region, dx, ct_cluster, lme_marker_estimates) %>%
        unnest(lme_marker_estimates) %>%
        mutate(
            gene = fct_rev(fct_relevel(gene, marker_gene_tb$gene)),
            pval_stars = pval_symnum(pval),
            stat_stars = statistic_symnum(z),
            ct_cluster = fct_drop(ct_cluster)
        ) %>%
        complete(ct_cluster, gene, dx) %>%
        filter(!is.na(region)) %>%
        glimpse

    browser()

    writeLines("plot")
    xwidth = (0.25 * length(unique(plot_gene_list$region)) * 3) + 10
    yheight = (0.25 * length(unique(plot_gene_list$gene))) + 2
    pdf(filepath, width = xwidth, height = yheight)
    plot_gene_list %>%
        group_split(region) %>%
        map(., fmt_ggplot, limits) %>%
        print
    graphics.off()
    
}

pval_symnum <- function(pval) {
    unclass(symnum(pval, 
        corr = FALSE, 
        na = FALSE,
        cutpoints = c(0, 0.0001, 0.005, 0.01, 0.05, 0.1, 1),
        symbols = c("****", "***", "**", "*", ".", " ")
    ))
}

statistic_symnum <- function(statistic) {
    unclass(symnum(statistic,
        cutpoints = c(-Inf, -2, 2, Inf),
        symbols = c("-", " ", "+")
    ))
}

fmt_ggplot <- function(plot_gene_list, limits) {
    plot_title <- str_glue("lme for snseq data")
    plot_ts <- str_glue("{system('md5sum seurat_layertype_dge_hlalocus.R', intern = TRUE)} {date()}")
    if (is.null(limits)) { 
        limits <- c(quantile(plot_gene_list$beta, 0.05), quantile(plot_gene_list$beta, 0.95))
    }
    ggplot_list <- list(
        geom_tile(),
        geom_text(),
        facet_wrap("dx"), 
        scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"), na.value = "grey50", limits = limits, oob = scales::squish),
        theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(), panel.background = element_rect(fill = "grey75"))
    )
    return(list(
        ggplot(plot_gene_list, aes(x = ct_cluster, y = gene, fill = beta, label = stat_stars)) +
            ggplot_list +
            labs(title = str_glue(plot_title, " -- t-statistic"), subtitle = plot_ts),
        ggplot(plot_gene_list, aes(x = ct_cluster, y = gene, fill = beta, label = pval_stars)) +
            ggplot_list +
            labs(title = str_glue(plot_title, " -- adj.P.Val"), subtitle = plot_ts)
    ))
}

if (!interactive()) {
    main()
}
