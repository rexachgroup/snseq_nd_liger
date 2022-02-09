liblist <- c("tidyverse")
l <- lapply(liblist, function(x) suppressPackageStartupMessages(require(x, character.only = TRUE, quietly = TRUE)))
csv_dir <- "/home2/jrexach/BFgrant_final_report"
out_path_base <- file.path(csv_dir, "figures")

args_tb <- tibble(
    cell_type = c("astrocyte", "excitatory", "inhibitory", "opc"),
    mouse_csv = file.path(csv_dir, c("DEA_AST_mt_nRNA_Sex_regressed.csv", "DEA_EX_mt_nRNA_Sex_regressed.csv", "DEA_INH_mt_nRNA_Sex_regressed.csv", "DEA_ODC_mt_nRNA_Sex_regressed.csv")),
    calcarine_csv = file.path(csv_dir, c("calcarine-astrocyte.csv", "calcarine-excitatory_trajectory.csv", "calcarine-inhibitory.csv", "calcarine-oligodendrocyte.csv")),
    insula_csv = file.path(csv_dir, c("insula-astrocyte.csv", "insula-excitatory.csv", "insula-inhibitory.csv", "insula-oligodendrocyte.csv")),
    precg_csv = file.path(csv_dir, c("preCG-astrocyte.csv", "preCG-excitatory.csv", "preCG-inhibitory.csv", "preCG-oligodendrocyte.csv"))
)

modules <- file.path(csv_dir, "geneInfo_TPR50_UnregressedData_Added.HumanConvert.csv")

main <- function() {
    module_tb <- read_csv(modules)
    turquoise_tb <- module_tb %>%
        select(ensembl = Human.Ensembl.Gene.ID, symbol = GeneSymbol, module = Initially.Assigned.Module.Color, kme = consensus.kMEturquoise)
    top_genes <- turquoise_tb %>%
        filter(module == "turquoise") %>%
        arrange(desc(kme)) %>%
        slice_max(kme, n = 10) %>%
        mutate(symbol = toupper(symbol))

    gene_tb <- pmap(args_tb, function(...) {
        cr <- list(...)
        mouse <- read_csv(cr$mouse_csv, 
            col_names = c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"), skip = 1) %>%
            mutate(gene = toupper(gene), cell_type = cr$cell_type)

        m_filter <- mouse %>%
            filter(gene %in% top_genes$symbol)
        glimpse(m_filter)
        return(m_filter)
    }) %>% bind_rows

    gg <- ggplot(gene_tb, aes(x = cell_type, y = avg_log2FC, label = signif(p_val_adj, 3), fill = cell_type)) + 
        geom_col() +
        geom_text(nudge_y = 0.001) +
        facet_wrap("gene", scales = "free_x", ncol = 4) +
        labs(title = "Turquoise module hits in mouse; bar height = log2FC, text = adjusted p value") +
        theme(legend.position = "none")

    pdf(file.path(out_path_base, "turquoise_expr.pdf"), width = 5 * 4, height = 5)
    print(gg)
    graphics.off()
}

if (!interactive()) {
    main()
}
