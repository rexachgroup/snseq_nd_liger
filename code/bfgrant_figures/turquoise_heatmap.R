liblist <- c("tidyverse", "ComplexHeatmap", "gtools", "circlize", "scales", "RColorBrewer")
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
        slice_max(kme, n = 50) %>%
        mutate(symbol = toupper(symbol))

    mouse_ast <- read_csv(file.path(csv_dir, "DEA_AST_mt_nRNA_Sex_regressed.csv"),
            col_names = c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"), skip = 1) %>%
            mutate(gene = toupper(gene), col = "mouse astrocyte") %>%
            select(gene = gene, estimate = avg_log2FC, pval = p_val_adj, col)
    insula_ast <- read_csv(file.path(csv_dir, "insula-astrocyte.csv")) %>%
            mutate(gene = toupper(gene), col = "insula ftd astrocyte") %>%
            select(gene = gene, estimate = clinical_dxbvFTD.p.value, pval = clinical_dxbvFTD.p.value.adj, col)
    mouse_odc <- read_csv(file.path(csv_dir, "DEA_ODC_mt_nRNA_Sex_regressed.csv"),
            col_names = c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"), skip = 1) %>%
            mutate(gene = toupper(gene), col = "mouse oligodendrocyte") %>%
            select(gene = gene, estimate = avg_log2FC, pval = p_val_adj, col)
    insula_odc <- read_csv(file.path(csv_dir, "insula-oligodendrocyte.csv")) %>%
            mutate(gene = toupper(gene), col = "insula ftd oligodendrocyte") %>%
            select(gene = gene, estimate = clinical_dxbvFTD.p.value, pval = clinical_dxbvFTD.p.value.adj, col)
    insula_mic <- read_csv(file.path(csv_dir, "insula-microglia.csv")) %>%
            mutate(gene = toupper(gene), col = "insula ftd microglia") %>%
            select(gene = gene, estimate = clinical_dxbvFTD.p.value, pval = clinical_dxbvFTD.p.value.adj, col)
    
    dge_tb <- bind_rows(mouse_ast, insula_ast, mouse_odc, insula_odc, insula_mic) %>%
        filter(gene %in% top_genes$symbol)

    estimate_mat <- pivot_matrix(dge_tb, "col", "estimate", "gene")
    pval_mat <- pivot_matrix(dge_tb, "col", "pval", "gene")
    stars_mat <- gtools::stars.pval(pval_mat)
    dimnames(stars_mat) <- dimnames(pval_mat)

    pdf(file.path(out_path_base, "turquoise_heatmap.pdf"), width = ncol(estimate_mat) * 0.75, height = nrow(estimate_mat) * 0.375)
    tryCatch(print(heatmap(estimate_mat, stars_mat)), error = print)
    graphics.off()

}

pivot_matrix <- function(tb, cols_from, values_from, rows_from) {
    tb_pivot <- tb %>%
        select(all_of(c(cols_from, values_from, rows_from))) %>%
        pivot_wider(names_from = all_of(cols_from), values_from = all_of(values_from)) %>%
        ungroup()

    tb_matrix <- select(tb_pivot, -all_of(rows_from)) %>%
        as.matrix()

    rownames(tb_matrix) <- tb_pivot %>% 
        rowwise() %>%
        summarize(join_rowname = paste(c_across(rows_from), collapse = "|"), .groups = "drop") %>% 
        pluck("join_rowname")

    return(tb_matrix)
}

heatmap <- function(heatmap_val_matrix, heatmap_label_matrix) { 
    plot_text_label <- function(j, i, x, y, w, h, col) {
        label <- heatmap_label_matrix[i, j]
        if (!is.na(label)) {
            grid.text(label, x, y)
        }
    }

    colormap_min <- quantile(heatmap_val_matrix, 0.05, na.rm = TRUE)
    colormap_max <- quantile(heatmap_val_matrix, 0.95, na.rm = TRUE)
    colormap <- colorRamp2(
        breaks = c(colormap_min, 0, colormap_max),
        colors = c(muted("blue"), "white", muted("red"))
    )

    heatmap_obj <- Heatmap(
        estimate_mat,
        col = colormap,
        na_col = "grey75",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        cell_fun = plot_text_label,
        column_names_rot = 90,
        name = " ",
    )
}
