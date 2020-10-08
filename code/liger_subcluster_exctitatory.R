liblist <- c("tidyverse", "Seurat", "ggpubr", "broom")
lapply(liblist, require, character.only = TRUE, quiet = TRUE)

in_seurat_metadata <- readRDS("../analysis/seurat_lchen/liger_subcluster_subset.rds")

plot_celltype_by_dxregion_boxplot <- function(seurat_meta){
  dat <- seurat_meta %>% 
    mutate(clinical_dx = fct_relevel(clinical_dx, "Control")) %>%
    group_by(library_id, clinical_dx, region, cluster_ids) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(clinical_dx, region, cluster_ids) %>%
    mutate(percent = n / sum(n))

  dat %>% group_split %>%
      pluck(1) %>% glimpse()
  
  my_comparisons <- list(c(1, 2), c(1, 3), c(1, 4))
  plt <- dat %>%
    ggplot(aes(x = clinical_dx, y = percent, fill = clinical_dx)) +
      geom_boxplot(outlier.shape = NA) +
      facet_grid("cluster_ids ~ region") +
      geom_jitter(color = "black", alpha = 0.5, position = position_jitter(seed = 1)) +
      stat_compare_means(aes(group = clinical_dx)
                         , comparisons = my_comparisons
                         , method = "t.test"
                         , p.adjust.method = "fdr"
                         , label = "p.format"
                         ) +
      coord_flip() +
      ggtitle("cell counts normalized per (dx, cluster) + comparison across dx")
  return(plt)
}

pdf(width = 15, height = 2 * 30)
tryCatch({in_seurat_metadata %>% 
  filter(cluster_cell_type == "excitatory") %>%
  plot_celltype_by_dxregion_boxplot},
  error = print,
  warning = print)
dev.off()


dat <- in_seurat_metadata %>%
    filter(cluster_cell_type == "excitatory") %>%
    mutate(clinical_dx = fct_relevel(clinical_dx, "Control")) %>%
    group_by(library_id, clinical_dx, region, cluster_ids) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(clinical_dx, region, cluster_ids) %>%
    mutate(percent = n / sum(n)) %>%
    arrange(region, clinical_dx, cluster_ids)

dat_ttest <- dat %>%
    group_by(region, cluster_ids) %>%
    group_nest() %>%
    mutate(pval = purrr::map(data, function(data) {
            if (nrow(data) > 1) {
                broom::tidy(pairwise.t.test(data$percent, data$clinical_dx, p.adjust.method = "fdr"))
            }
        })
    ) %>%
    unnest(pval) %>%
    select(-data) %>%
    filter(group2 == "Control")
