liblist <- c("here", "tidyverse", "WGCNA", "BisqueRNA", "Biobase", "broom", "patchwork", "ggpubr")
l <- suppressPackageStartupMessages(
    lapply(liblist, library, character.only = TRUE, quietly = TRUE)
)

LMHCL_PATH <- here("Analysis_lchen_consensus/")
LMHCL_EXPR_LIST <- file.path(LMHCL_PATH, "Step01_PCA_OutlierRemoval/datExpr_regions_rmsampl.rds")
LMHCL_CONSENSUS_TREE <- file.path(LMHCL_PATH, "Step06_WGCNA/consensusTOM/consensusTree.rds")
GENE_NAMES <- file.path(LMHCL_PATH, "Step01_PCA_OutlierRemoval/geneNames.rds")
LAKE_FC <- file.path(LMHCL_PATH, "Step00_InputFiles/Lake2018_FrontalCortex_forEWCE.RData")
LAKE_VC <- file.path(LMHCL_PATH, "Step00_InputFiles/Lake2018_VisualCortex_forEWCE.RData")
ND_DATA <- file.path(LMHCL_PATH, "Step01_PCA_OutlierRemoval/individual_nd_filtered.rds")
source(here("Scripts_lchen/Analysis/COMMON_DEFS.R"))
OUT_DIR <- file.path(LMHCL_PATH, "Step09_Bisque/")

plot_ts <- function() {
    nm <- "bisque_lake_psp.R"
    str_glue("{nm} {date()} {tools::md5sum(nm)}")
}

#pal_stallion <- c("Control"="#D51F26","AD"="#272E6A","PiD"="#208A42","PSP"="#89288F")

main <- function() {
    dir.create(OUT_DIR)
    expr_list <- readRDS(LMHCL_EXPR_LIST)
    consensus_tree <- readRDS(LMHCL_CONSENSUS_TREE)
    consensus_meta <- rbindlist_mtd(consensus_tree$meta, "sample")
    individual_nd_filtered <- readRDS(ND_DATA) %>%
        select(Autopsy.ID, region, type, score) %>%
        group_by(Autopsy.ID, region) %>%
        group_nest()
    gene_names <- readRDS(GENE_NAMES)
    meta_nd <- left_join(consensus_meta, individual_nd_filtered, by = c("Autopsy.ID", "Region" = "region")) %>%
        mutate(Primary.Neuropath.Dx = fct_recode(Primary.Neuropath.Dx, c("FTD" = "PiD")) %>% fct_relevel("Control"))

    ba10_env <- new.env()
    ba17_env <- new.env()
    load(LAKE_FC, ba10_env)
    load(LAKE_VC, ba17_env)
    lake_fc <- make_es(ba10_env$ba10_umi, ba10_env$ba10_datMeta)
    lake_vc <- make_es(ba17_env$ba17_umi, ba17_env$ba17_datMeta)

    expr_region_es <- map(expr_list, make_data_es, gene_names)

    bisque_lake_fc <- map(expr_region_es, lake_bisque_wrap, lake_fc)
    bisque_lake_vc <- map(expr_region_es, lake_bisque_wrap, lake_vc)

    saveRDS(bisque_lake_fc, file.path(OUT_DIR, "bisqueRNA_expression_fc.rds"))
    saveRDS(bisque_lake_vc, file.path(OUT_DIR, "bisqueRNA_expression_vc.rds"))
    bisque_lake_fc <- readRDS(file.path(OUT_DIR, "bisqueRNA_expression_fc.rds"))
    bisque_lake_vc <- readRDS(file.path(OUT_DIR, "bisqueRNA_expression_vc.rds"))
    dt_fc <- bisque_tb_multi(bisque_lake_fc) %>%
        left_join(meta_nd, by = "sample")
    dt_vc <- bisque_tb_multi(bisque_lake_vc) %>%
        left_join(meta_nd, by = "sample")
    write_csv(dt_fc, file.path(OUT_DIR, "bisqueRNA_expession_fc.csv"))
    write_csv(dt_vc, file.path(OUT_DIR, "bisqueRNA_expression_vc.csv"))

    fc_test <- dt_fc %>%
        filter(subpopulation == "Ast") %>%
        bisque_pairwise("bulk_props", "Primary.Neuropath.Dx", c("Region", "subpopulation"))
    vc_test <- dt_vc %>%
        filter(!is.na(type), subpopulation == "Ast") %>%
        bisque_pairwise("bulk_props", "Primary.Neuropath.Dx", c("Region", "subpopulation"))

    # celltype bulk.props plot + dx t-test
    pdf(file.path(OUT_DIR, "ggplot_bisque_celltype_fc.pdf"), width = 30, height=7)
    print(ggplot_bisque_ttest(fc_test, "Lake BA10 bulk.props t.test, celltype ="))
    graphics.off() 

    pdf(file.path(OUT_DIR, "ggplot_bisque_celltype_vc.pdf"), width = 28, height=7)
    print(ggplot_bisque_ttest(vc_test, "Lake BA17 bulk.props t.test, celltype ="))
    graphics.off()
    
    saveRDS(fc_test, file.path(OUT_DIR, "bisque_celltype_fc_ttest.rds"))
    saveRDS(vc_test, file.path(OUT_DIR, "bisque_celltype_vc_ttest.rds"))
    fc_test %>% 
        select(Region, subpopulation, bulk_tidy, bulk_n, fdr_p_value) %>%
        unnest(c(bulk_tidy, fdr_p_value)) %>%
        write_csv(file.path(OUT_DIR, "bisque_celltype_fc_ttest.csv"))
    vc_test %>%
        select(Region, subpopulation, bulk_tidy, bulk_n, fdr_p_value) %>%
        unnest(c(bulk_tidy, fdr_p_value)) %>%
        write_csv(file.path(OUT_DIR, "bisque_celltype_vc_ttest.csv"))

    # celltype bulk.props / dx wilcoxon
    fc_wc_test <- dt_fc %>%
        filter(subpopulation == "Ast") %>%
        group_test_pairwise("bulk_props", "Primary.Neuropath.Dx", c("Region", "subpopulation"))
    
    vc_wc_test <- dt_vc %>%
        filter(subpopulation == "Ast") %>%
        group_test_pairwise("bulk_props", "Primary.Neuropath.Dx", c("Region", "subpopulation"))
    saveRDS(fc_wc_test, file.path(OUT_DIR, "bisque_celltype_fc_wilcoxon.rds"))
    saveRDS(vc_wc_test, file.path(OUT_DIR, "bisque_celltype_vc_wilcoxon.rds"))
    write_csv(fc_wc_test, file.path(OUT_DIR, "bisque_celltype_fc_wilcoxon.csv"))
    write_csv(vc_wc_test, file.path(OUT_DIR, "bisque_celltype_vc_wilcoxon.csv"))
    
    pdf(file.path(OUT_DIR, "ggplot_bisque_celltype_fc_wilcoxon.pdf"), width = 28, height = 7)
    print(ggplot_bisque_wilcoxon(fc_wc_test, "Lake BA10 bulk.props wilcoxon + fdr_p_value, celltype ="))
    graphics.off()

    pdf(file.path(OUT_DIR, "ggplot_bisque_celltype_vc_wilcoxon.pdf"), width = 28, height = 7)
    print(ggplot_bisque_wilcoxon(vc_wc_test, "Lake BA17 bulk.props wilcoxon + fdr_p_value, celltype ="))
    graphics.off()

    # nd score vs. bulk proportion test
    # for samples with valid score, do wilcox of nd vs. bulk proportion within each (region, dx, celltype) group.
    # fdr correct across number of (region, dx, celltype) per nd type.
    fc_score_test <- dt_fc %>% unnest(data) %>%
        filter(subpopulation == "Ast", !is.na(type)) %>%
        group_by(type) %>% group_nest() %>%
        mutate(test_tb = map(data, ~group_test_wilcox(.x, "score", "bulk_props", c("Region", "Primary.Neuropath.Dx", "subpopulation")))) %>%
        select(-data) %>% unnest(test_tb)

    vc_score_test <- dt_vc %>% unnest(data) %>%
        filter(subpopulation == "Ast", !is.na(type)) %>%
        group_by(type) %>% group_nest() %>%
        mutate(test_tb = map(data, ~group_test_wilcox(.x, "score", "bulk_props", c("Region", "Primary.Neuropath.Dx", "subpopulation")))) %>%
        select(-data) %>% unnest(test_tb)

    write_csv(fc_score_test, file.path(OUT_DIR, "bisque_fc_nd_wilcox.csv"))
    write_csv(vc_score_test, file.path(OUT_DIR, "bisque_vc_nd_wilcox.csv"))


    fc_ast_neu_props <- ct_proportion_calc(dt_fc, "Ast", "Ex*|In*") %>%
        unnest(cor_sum)
    vc_ast_neu_props <- ct_proportion_calc(dt_vc, "Ast", "Ex*|In*") %>%
        unnest(cor_sum)
    write_csv(fc_ast_neu_props, file.path(OUT_DIR, "bisque_fc_ast_neu_props.csv"))
    write_csv(vc_ast_neu_props, file.path(OUT_DIR, "bisque_vc_ast_neu_props.csv"))

    pdf(file.path(OUT_DIR, "ggplot_fc_bisque_celltype_scatter.pdf"), width = 28, height = 7)
    print(ct_proportion_comp(dt_fc, "Ast", "Ex*|In*", "Bisque BA10 Ast v. summed neuron bulk.props"))
    graphics.off()
    
    pdf(file.path(OUT_DIR, "ggplot_vc_bisque_celltype_scatter.pdf"), width = 28, height = 7)
    print(ct_proportion_comp(dt_vc, "Ast", "Ex*|In*", "Bisque BA17 Ast v. summed neuron bulk.props"))
    graphics.off()
    

}

make_es <- function(umi, meta) {
    varMeta <- data.frame(labelDescription = colnames(meta), row.names = colnames(meta))
    meta <- as.data.frame(meta)
    rownames(meta) <- meta$SampleName2
    phenoData <- new("AnnotatedDataFrame", data = meta, varMetadata = varMeta)
    es <- ExpressionSet(assayData = as.matrix(umi), phenoData = phenoData)
    return(es)
}

make_data_es <- function(dat_expr, gene_names) {
    gene_info_subset <- gene_names %>%
        filter(gene_id %in% rownames(dat_expr)) %>%
        filter(!duplicated(gene_name))
    dat_expr <- dat_expr[gene_info_subset$gene_id, ]
    rownames(dat_expr) <- gene_info_subset$gene_name
    dat_expr_es <- Biobase::ExpressionSet(assayData = as.matrix(dat_expr))
    return(dat_expr_es)
}

lake_bisque_wrap <- function(bulk_es, lake_es) {
	BisqueRNA::ReferenceBasedDecomposition(
		bulk_es, 
		lake_es,
		use.overlap = FALSE,
		cell.types = "Identity",
		subject.names = "Patient UMB#"
	)
}

bisque_tb_multi <- function(bisque_list) {
    map(bisque_list, function(x) {
        x$bulk.props %>%
            as.data.frame %>%
            rownames_to_column("subpopulation") %>%
            pivot_longer(cols = -"subpopulation", names_to = "sample", values_to = "bulk_props")
    }) %>%
    bind_rows
}

rbindlist_mtd <- function(x, rownames = "rownames") {
    map(x, function(x1) {
        x1$data %>%
            as.data.frame %>%
            rownames_to_column(rownames)
    }) %>% bind_rows
}

group_test_wilcox <- function(tb, var1, var2, group_cols) {
   tb_test <- tb %>%
        group_by(across(all_of(group_cols))) %>% group_nest() %>%
        mutate(
            test_res = pmap(., function(...) {
                cr <- list(...)
                var1 <- as.numeric(cr$data[[var1]])
                var2 <- as.numeric(cr$data[[var2]])
                tryCatch({
                    res = tidy(wilcox.test(var1, var2, paired = T))
                    return(res)
                },
                error = function(x) {print(x); return(NA)})
            })
        ) %>%
        unnest(test_res) %>%
        mutate(
            n = n(),
            fdr_p_value = p.adjust(p.value, method = "fdr", n = n())
        )
}

group_test_pairwise <- function(tb, var, group, nest_cols) {
    tb %>%
        group_by(across(all_of(nest_cols))) %>% group_nest() %>%
        mutate(
            test_res = pmap(., function(...) {
                cr <- list(...)
                var <- as.numeric(cr$data[[var]])
                group <- cr$data[[group]]
                tidy(pairwise.wilcox.test(x = var, g = group, p.adjust.method = "none"))
            })
        ) %>%
        unnest(test_res) %>%
        filter(group2 == levels(tb[[group]])[1]) %>% # only control tests
        mutate(
            test_n = n(),
            fdr_p_value = p.adjust(p.value, method = "fdr", n = test_n)
        )
}

bisque_pairwise <- function(bisque_dt, bulk_col, test_group_col, nest_cols) {
    bisque_dt_wk <- bisque_dt %>%
        group_by_at(nest_cols) %>%
        group_nest() %>%
        mutate(
            bulk_test = pmap(., function(data = data, ...){
                bulk <- data[[ bulk_col ]]
                group <- data[[ test_group_col ]]
                pairwise.t.test(bulk, group, p.adjust.method = "none")
            }),
            bulk_avg = pmap_dbl(., function(data = data, ...) mean(data[[ bulk_col ]])),
        )
    bisque_dt_wk <- bisque_dt_wk %>%
        mutate(
            bulk_tidy = map(bulk_test, function(x) { 
                tryCatch(tidy(x), error = function(x) { return(NA) }) }),
            bulk_sum = map(bulk_test, summary),
        ) %>%
        filter(!is.na(bulk_tidy)) %>%
        mutate(bulk_n = n()) %>%
        mutate(
            fdr_p_value = pmap(., function(bulk_n = bulk_n, bulk_tidy = bulk_tidy, ...) {
                p.adjust(bulk_tidy$p.value, n = bulk_n, method = "fdr")
            })
        )
}

bisque_lm_nd_fit <- function(bisque_dt) {
   bisque_dt %>%
       filter(!is.na(type)) %>%
       group_by(subpopulation,Primary.Neuropath.Dx, type) %>%
       group_nest(keep = T) %>%
       mutate(lm_dat = pmap(. ,function(...) {
            cr <- list(...)
            form <- "bulk_props ~ Region + score"
            tryCatch(
                tidy(lm(as.formula(form), cr$data)) %>% mutate(formula = form), 
            error = function(x) {NA})
       }))
}

pair_ext <- function(dt) {
    dt %>%
        unnest(cols = c(bulk_tidy, fdr_p_value)) %>%
        filter(group2 == "Control") %>%
        select(-data, -bulk_sum, -bulk_test)
}

ggplot_bisque_ttest <- function(bisque_ttest_tb, title_prefix) {
    ts <- plot_ts()
    bracket_data <- pair_ext(bisque_ttest_tb)
    plots <- bisque_ttest_tb %>%
        group_by(subpopulation, Region) %>%
        group_nest() %>%
        mutate(plot = pmap(., function(...) {
            cr <- list(...)
            point_data <- cr$data$data[[1]]
            bracket_plot <- filter(bracket_data, subpopulation %in% cr$subpopulation, Region %in% cr$Region) %>%
                mutate(label = signif(p.value, 2))
            gg <- ggplot(point_data, aes(x = Primary.Neuropath.Dx, y = bulk_props, fill = Primary.Neuropath.Dx)) +
                geom_violin() +
                geom_point(position = position_jitter(width = 0.01)) +
                geom_bracket(
                    data = bracket_plot, 
                    aes(xmin = group2, xmax = group1, label = label, y.position = max(point_data$bulk_props)),
                    step.increase = 0.05,
                    inherit.aes = F
                ) +
                ggtitle(cr$Region)
            return(gg)
        }))
    
    plots %>%
        ungroup %>%
        group_by(subpopulation) %>% 
        group_map(function(.x, .y) {
            wrap_plots(.x$plot, ncol = length(unique(.x$Region))) +
                plot_layout(guides = "collect") +
                plot_annotation(title = str_glue("{title_prefix} {.y$subpopulation}"), subtitle = ts)
        })
}

ggplot_bisque_wilcoxon <- function(bisque_wilcox_tb, title_prefix) {
    ts <- plot_ts()
    plots <- bisque_wilcox_tb %>%
        group_by(subpopulation, Region) %>%
        group_nest() %>%
        mutate(plot = pmap(., function(...) {
            cr <- list(...)
            point_data <- cr$data$data[[1]]
            bracket_plot <- cr$data %>%
                mutate(label = signif(fdr_p_value, 2))
            gg <- ggplot(point_data, aes(x = Primary.Neuropath.Dx, y = bulk_props, fill = Primary.Neuropath.Dx)) +
                geom_violin() +
                geom_point(position = position_jitter(width = 0.01)) +
                geom_bracket(
                    data = bracket_plot, 
                    aes(xmin = group2, xmax = group1, label = label, y.position = max(point_data$bulk_props)),
                    step.increase = 0.05,
                    inherit.aes = F
                ) +
                ggtitle(cr$Region)
            return(gg)
        }))
    
    plots %>%
        ungroup %>%
        group_by(subpopulation) %>% 
        group_map(function(.x, .y) {
            wrap_plots(.x$plot, ncol = length(unique(.x$Region))) +
                plot_layout(guides = "collect") +
                plot_annotation(title = str_glue("{title_prefix} {.y$subpopulation}"), subtitle = ts)
        })
}

ct_proportion_calc <- function(bisque_meta, ct1, ct2) {
    bisque_meta %>%
        group_by(Region) %>%
        group_nest(keep = F) %>%
        mutate(cor_sum = pmap(., function(...) {
            cr <- list(...)
            ct1_rename <- setNames(c("bulk_props"), paste0(ct1, "_bulk_props"))
            ct2_rename <- setNames(c("bulk_props"), paste0(ct2, "_bulk_props"))
            meta <- cr$data %>% select(sample, Source, Primary.Neuropath.Dx, Autopsy.ID) %>% filter(!duplicated(sample))
            ct1_tb <- cr$data %>%
                filter(str_detect(subpopulation, ct1)) %>%
                group_by(sample) %>%
                summarize(bulk_props = sum(bulk_props), .groups = "drop") %>%
                rename(all_of(ct1_rename))
            ct2_tb <- cr$data %>%
                filter(str_detect(subpopulation, ct2))  %>%
                group_by(sample) %>%
                summarize(bulk_props = sum(bulk_props), .groups = "drop") %>%
                rename(all_of(ct2_rename))
            join_tb <- inner_join(ct1_tb, ct2_tb, by = "sample") %>%
                left_join(meta, by = "sample")
                
            return(join_tb)
        }))

}

ct_proportion_comp <- function(bisque_meta, ct1, ct2, title_prefix) {
    ts <- plot_ts()
    plots <- bisque_meta %>%
        group_by(Region) %>%
        group_nest() %>%
        mutate(plot = pmap(., function(...) {
            cr <- list(...)
            ct1_tb <- cr$data %>%
                filter(str_detect(subpopulation, ct1)) %>%
                group_by(Autopsy.ID) %>%
                summarize(bulk_props = sum(bulk_props), .groups = "drop")
            ct2_tb <- cr$data %>%
                filter(str_detect(subpopulation, ct2))  %>%
                group_by(Autopsy.ID) %>%
                summarize(bulk_props = sum(bulk_props), .groups = "drop")
            plot_tb <- inner_join(ct1_tb, ct2_tb, by = "Autopsy.ID") %>%
                left_join(cr$data, by = "Autopsy.ID")
            test_res <- tidy(cor.test(ct1_tb$bulk_props, ct2_tb$bulk_props))
            ggplot(plot_tb, aes(x = `bulk_props.x`, y = `bulk_props.y`, color = Primary.Neuropath.Dx)) +
                geom_point() +
                labs(
                    title = str_glue("{cr$Region} Bisque Ast v. Neuron"),
                    subtitle = str_glue("(cor.test estimate = {signif(test_res$estimate, 2)}, p.vaue = {signif(test_res$p.value, 2)})"), 
                    x = str_glue("{ct1} summed bulk.props"),
                    y = str_glue("{ct2} summed bulk.props")
                )
        }))
    wrap_plots(plots$plot, ncol = 7) +
        plot_layout(guides = "collect") + 
        plot_annotation(title = str_glue("{title_prefix}"), subtitle = ts)
}

main()
