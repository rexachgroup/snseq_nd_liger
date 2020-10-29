set.seed(27)
require(Seurat)
require(cowplot)
require(viridis)
require(tidyverse)
require(RColorBrewer)
require(future.apply)
require(ggpubr)
require(patchwork)
require(lme4)
require(lmerTest)
require(broom.mixed)
require(ComplexHeatmap)
source("function_library.R")
source("ggplot_theme.R")
source("seurat_function_library.R")
options(future.globals.maxSize = Inf)

in_model <- "../analysis/clinical_dx_de_lme_microglia/20201026/tables/clinical_dx_de_lme_microglia_lm_tb.rds" 
in_seurat_rds <-
  "../analysis/pci_import/20201028/tables/pci_seurat.rds"

out_graph <- file.path(dirname(in_model), "../graphs/")
if (!dir.exists(out_graph)) dir.create(out_graph)
out_graph <- normalizePath(out_graph)

nd_so <- readRDS(in_seurat_rds)
cell_type <- "microglia"
region <- "preCG"

# Filter metadata to get cell ids for matching region, cluster_cell_type.
# Subset seurat object.
cell_id_subset <- nd_so@meta.data %>%
    filter(cell_type == {{cell_type}},
           region == {{region}}) %>%
    rownames
mc_so <- subset(nd_so, cells = cell_id_subset)
meta <- mc_so@meta.data
rm(nd_so)
gc()

lm_tb <- readRDS(in_model)

clinical_dx <- c("bvFTD", "PSP-S", "AD")
plot_vars <- c(
      "pmi",
      "age",
      "rin",
      "sex",
      "number_umi",
      "percent_mito",
      "prep",
      "finalsite")


up_plot_tbs <- lm_tb %>%
    group_split(term) %>%
    map(function(tb) {
        dx <- str_match(unique(tb$term), "clinical_dx(.*)")[1, 2]
        writeLines(str_glue(" ========== {dx} =========== "))
        up_tb <- tb %>% filter(estimate > 0.1, t_statistic > 2)

        meta_dx_subset <- mc_so@meta.data %>%
            as_tibble() %>%
            filter(clinical_dx == {{dx}} | clinical_dx == "Control") %>%
            select(clinical_dx, cell_ids, all_of(plot_vars))
            
        up_so <- subset(mc_so, features = up_tb$gene, cells = meta_dx_subset$cell_ids)

        up_mat <- GetAssayData(up_so) %>%
            as_tibble(rownames = "gene") %>%
            pivot_longer(cols = -gene,
                names_to = "cell_ids",
                values_to = "expr")

        plot_tb <- inner_join(meta_dx_subset, up_mat, by = "cell_ids")
    })


up_aveExpr_plots <- map(up_plot_tbs, function(up_tb) {
    dx <- unique(up_tb$clinical_dx)
    title_str <- str_glue("LME of {dx[[1]]} v. {dx[[2]]} in PreCG microglia")
    subtitle_str <- str_glue("average expression per dx of significant upregulated genes \n (beta > 0.1, t_statistic > 2)")
    ave_up <- up_tb %>%
        group_by(gene, clinical_dx) %>%
        summarize(ave_expr = mean(expr))
    gg <- ggplot(ave_up, aes_string(x = "clinical_dx", y = "ave_expr", fill = "clinical_dx")) +
        geom_boxplot() + 
        stat_compare_means(paired = FALSE) +
        labs(title = title_str, subtitle = subtitle_str)

    return(gg)
})

up_ggplots <- map(up_plot_tbs, function(up_tb) {
        dx <- unique(up_tb$clinical_dx)
        title_str <- str_glue("LME of {dx[[1]]} v. {dx[[2]]} in PreCG microglia")
        subtitle_str <- str_glue("expression of significant upregulated genes (beta > 0.1, t_statistic > 2) vs. metadata")
        plots <- map(plot_vars, function(plot_var) {
            plot_subset <- up_tb %>%
                select("clinical_dx", "gene", "expr", {{plot_var}})
            gg <- ggplot(plot_subset, aes_string(x = {{plot_var}}, y = "expr"))
            if (is.numeric(plot_subset[[plot_var]]))
                gg <- gg + geom_point(aes_string(color = "clinical_dx")) +
                    geom_smooth(aes_string(group = "clinical_dx", fill = "clinical_dx"), method = "lm", se = FALSE)
            else
                gg <- gg + geom_boxplot(aes_string(fill = "clinical_dx")) +
                    stat_compare_means(paired = FALSE)
        })
        wrap_plots(plots, ncol = 2) +
            plot_annotation(title = title_str, subtitle = subtitle_str)
    })

system.time({
    png(file.path(out_graph, "upreg-dx-%d.png"), units = "in", width = 14, height = length(plot_vars) * 2, pointsize = 20, res = 300)
    walk(up_ggplots, plot)
    dev.off()
})

system.time({
    png(file.path(out_graph, "aveexpr-dx-%d.png"), units = "in", width = 7, height = 7, pointsize = 20, res = 300)
    walk(up_aveExpr_plots, plot)
    dev.off()
})
