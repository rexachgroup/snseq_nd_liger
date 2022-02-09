liblist <- c("tidyverse", "ggcorrplot", "GGally", "ggrepel")
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
pval_threshold <- 0.05

main <- function() {
    pwalk(args_tb, function(...) {
        cr <- list(...)
        mouse <- read_csv(cr$mouse_csv, 
            col_names = c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"), skip = 1) %>%
            mutate(gene = toupper(gene))
        calcarine <- read_csv(cr$calcarine_csv) %>%
            mutate(gene = toupper(gene))
        insula <- read_csv(cr$insula_csv) %>%
            mutate(gene = toupper(gene))
        precg <- read_csv(cr$precg_csv) %>%
            mutate(gene = toupper(gene))

        m_filter <- mouse %>%
            filter(p_val_adj < pval_threshold)
        c_filter <- calcarine %>%
            filter(clinical_dxbvFTD.p.value.adj < pval_threshold)
        i_filter <- insula %>%
            filter(clinical_dxbvFTD.p.value.adj < pval_threshold)
        p_filter <- precg %>%
            filter(clinical_dxbvFTD.p.value.adj < pval_threshold)

        mc_join <- inner_join(m_filter, c_filter, by = "gene") %>%
            select(gene = gene, mouse_estimate = avg_log2FC, calcarine_estimate = clinical_dxbvFTD.estimate)
        mi_join <- inner_join(m_filter, i_filter, by = "gene") %>%
            select(gene = gene, mouse_estimate = avg_log2FC, insula_estimate = clinical_dxbvFTD.estimate)
        mp_join <- inner_join(m_filter, p_filter, by = "gene") %>%
            select(gene = gene, mouse_estimate = avg_log2FC, precg_estimate = clinical_dxbvFTD.estimate)

        mc_cortest <- tryCatch(cor.test(mc_join$mouse_estimate, mc_join$calcarine_estimate), error = function(x) { return(list(estimate = NA, p.value = NA)) })
        mi_cortest <- tryCatch(cor.test(mi_join$mouse_estimate, mi_join$insula_estimate), error = function(x) { return(list(estimate = NA, p.value = NA)) })
        mp_cortest <- tryCatch(cor.test(mp_join$mouse_estimate, mp_join$precg_estimate), error = function(x) { return(list(estimate = NA, p.value = NA)) })

        pdf(file.path(out_path_base, str_glue("scatter_{cr$cell_type}.pdf")))
        mc_plot <- ggplot(mc_join, aes(x = mouse_estimate, y = calcarine_estimate, label = gene)) +
            geom_point() +
            geom_text_repel() +
            ggtitle(str_glue("{cr$cell_type} calcarine Pearson's correlation={signif(mc_cortest$estimate, 3)}, pval={signif(mc_cortest$p.value, 3)}"))
        mi_plot <- ggplot(mi_join, aes(x = mouse_estimate, y = insula_estimate, label = gene)) + 
            geom_point() +
            geom_text_repel() +
            ggtitle(str_glue("{cr$cell_type} insula Pearson's correlation={signif(mi_cortest$estimate, 3)}, pval={signif(mi_cortest$p.value, 3)}"))
        mp_plot <- ggplot(mp_join, aes(x = mouse_estimate, y = precg_estimate, label = gene)) + 
            geom_point() +
            geom_text_repel() +
            ggtitle(str_glue("{cr$cell_type} precg Pearson's correlation={signif(mp_cortest$estimate, 3)}, pval={signif(mp_cortest$p.value, 3)}"))
        print(mc_plot)
        print(mi_plot)
        print(mp_plot)

        dev.off()
        return()
    })
}

function() {
        m_filter <- mouse %>%
            filter(p_val_adj < pval_threshold)
        p_filter <- precg %>%
            filter(clinical_dxbvFTD.p.value.adj < pval_threshold)
        mp_join <- inner_join(m_filter, p_filter, by = "gene") %>%
            select(gene = gene, mouse_estimate = avg_log2FC, precg_estimate = clinical_dxbvFTD.estimate)
        mp_cortest <- tryCatch(cor.test(mp_join$mouse_estimate, mp_join$precg_estimate), error = function(x) { return(list(estimate = NA, p.value = NA)) }) 
        mp_plot <- ggplot(mp_join, aes(x = mouse_estimate, y = precg_estimate, label = gene)) + 
            geom_point() +
            geom_text_repel() +
            ggtitle(str_glue("{cr$cell_type} precg Pearson's correlation={signif(mp_cortest$estimate, 3)}, pval={signif(mp_cortest$p.value, 3)}"))
}

if (!interactive()) {
    main()
}
