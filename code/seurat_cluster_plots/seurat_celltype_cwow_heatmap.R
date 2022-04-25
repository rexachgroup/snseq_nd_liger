# Heatmap results from seurat_cluster_haplotype_dge.R
# Filter by cwow tau genelist.
liblist <- c("tidyverse", "patchwork")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

HAPLOTYPE_DGE_RDS <- "../../analysis/seurat_lchen/seurat_cluster_lme/seurat_haplotype_dge/subcluster_lme.rds"
GENE_LIST <- "../../resources/cwow_genelist_WITHNOTES.csv"
OUT_DIR <- "../../analysis/seurat_lchen/seurat_cluster_lme/seurat_haplotype_plot/"
plot_ts <- str_glue("{system('md5sum seurat_celltype_cwow_heatmap.R', intern = TRUE)} {date()}")

main <- function() {
    dir.create(OUT_DIR, recursive = TRUE)
    subcluster_tb <- readRDS(HAPLOTYPE_DGE_RDS)
    genelist <- read_csv(GENE_LIST)
    out_pdf <- file.path(OUT_DIR, "tau_genes.pdf")
    subcluster_tb <- subcluster_tb %>%
        filter(!cell_type %in% c("pericyte", "t_cell")) %>%
        inner_join(genelist, by = "gene") %>%
        mutate(gene = fct_relevel(gene, rev(genelist$gene)))

    plot_tb <- subcluster_tb %>%
        mutate(stars = map_chr(`Tau_HH1/H2.statistic`, stars_label)) %>%
        select(x = cell_type, y = gene, fill = `Tau_HH1/H2.estimate`, facet = region, stars)
    plot_gg <- ggplot(plot_tb, aes(x = x, y = y, fill = fill)) +
        geom_tile() +
        facet_wrap("facet") + 
        scale_fill_gradient2(
            low = scales::muted("blue"),
            mid = "white",
            high = scales::muted("red"),
            limits = c(-0.25, 0.25),
            oob = scales::squish,
        ) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        geom_text(aes(label = stars)) +
        labs(
            title = str_glue("LME DGE on Tau H1/H2 haplotype"),
            subtitle = str_glue(
                plot_ts,
                "\n",
                "beta -> color, ",
                "- -> t_statistic < -2, ",
                "+ -> t_statistic > 2"
            )
        )
    yheight <- (length(unique(plot_tb$y)) * 0.15) + 5
    
    pdf(out_pdf, height = yheight, width = 10)
    tryCatch(print(plot_gg), error = print)
    graphics.off()
}

stars_label <- function(vec) {
    as.character(symnum(vec,
        cutpoints = c(-Inf, -2, 2, Inf),
        symbols = c("-", " ", "+")))
}

if (!interactive()) main()
