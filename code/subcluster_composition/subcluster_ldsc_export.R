# export significant genes from subclusters with significant dx-specific cell composition changes
liblist <- c("biomaRt", "tidyverse")
lapply(liblist, require, character.only = TRUE)
dge_dir <- "../../analysis/clinical_dx_lme_and_lm_liger/20200930_2/tables/"
dge_files <- list.files(dge_dir, recursive = TRUE, pattern = "*.csv", full.names = TRUE)

#clinical_dx_de_lme_microglia <- readRDS(
limma_sc_dx_pmi <- read_rds("../../analysis/seurat_lchen/subcluster_composition/limma/subcluster_composition_dx_pmi.rds")
META_FILE <- "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"
out_dir <- "../../analysis/seurat_lchen/tables/"

SUBCLUSTER_FILTER_FILE <- "../../analysis/seurat_lchen/liger_subcluster_filtered_props.rds"
if (!dir.exists(out_dir)) dir.create(out_dir)

limma_sc_unnest <- limma_sc_dx_pmi %>%
    unnest(limma_pval) %>%
    glimpse

limma_sc_ad <- limma_sc_unnest %>%
    filter(`dxAD - (dxControl + dxbvFTD + dxPSP.S) / 3` < 0.05) %>%
    pluck("ct_subcluster")
limma_sc_ftd <- limma_sc_unnest %>%
    filter(`dxbvFTD - (dxControl + dxAD + dxPSP.S) / 3` < 0.05) %>%
    pluck("ct_subcluster")
limma_sc_psp <- limma_sc_unnest %>%
    filter(`dxPSP.S - (dxControl + dxAD + dxbvFTD) / 3` < 0.05) %>%
    pluck("ct_subcluster")

mart <- useMart('ENSEMBL_MART_ENSEMBL', 'hsapiens_gene_ensembl')
mart_map <- getBM(attributes = c('entrezgene_id', 'ensembl_gene_id'),
						 mart = mart)
mart_map <- as_tibble(mart_map)
dge_tables <- map(dge_files, read_csv)
dge_tb <- bind_rows(dge_tables) %>%
    mutate(ct_subcluster = paste(region, cell_type, liger_cluster, sep = "-")) %>%
    left_join(mart_map, by = c("ensembl" = "ensembl_gene_id")) %>%
    drop_na(ensembl)

good_modules <- c("insula-microglia-7",
                  "preCG-microglia-1",
                  "preCG-microglia-0",
                  "preCG-microglia-3")

t_statistic_slice <- function(tb, n, mode = "up") {
    if (mode == "down") {
        tb %>% group_by(ct_subcluster) %>%
            slice_min(t_statistic_ad_lme, n = n)
    } else if (mode == "up") {
        tb %>% group_by(ct_subcluster) %>%
            slice_max(t_statistic_ad_lme, n = n) 
    }
}

dge_tb_list <- list(
    ad_up = dge_tb %>% t_statistic_slice(1000, "up")  %>% filter(ct_subcluster %in% good_modules),
    ad_down = dge_tb %>% t_statistic_slice(1000, "down") %>% filter(ct_subcluster %in% good_modules)
)
# dge_tb_list <- list(
#     ad_up = dge_tb %>% filter(t_statistic_ad_lme > 2),
#     ad_down = dge_tb %>% filter(t_statistic_ad_lme < -2),
#     ftd_up = dge_tb %>% filter(t_statistic_ftd_lme > 2),
#     ftd_down = dge_tb %>% filter(t_statistic_ftd_lme < -2),
#     psp_up = dge_tb %>% filter(t_statistic_psp_lme > 2),
#     psp_down = dge_tb %>% filter(t_statistic_psp_lme < -2)
# 

dge_genelist <- imap(dge_tb_list, function(x, nm) {
    x %>%
        mutate(subcluster_dx = paste(region, cell_type, liger_cluster, nm, sep = "_")) %>%
        dplyr::select(ensembl, entrezgene_id, gene, ct_subcluster, subcluster_dx, contains("t_statistic"))
})

dge_all <- bind_rows(dge_genelist, .id = "dge_type") %>% ungroup


# jaccardBinary <- function(vec1, vec2) {
#     return( length(intersect(vec1, vec2)) / length(union(vec1, vec2)) )
# }
# 
# subclusters <- unique(dge_all$ct_subcluster)
# cat_combn <- t(combn(subclusters, m = 2)) %>%
#     as_tibble()
# colnames(cat_combn) <- c("cat1", "cat2")
# dge_jaccard <- cat_combn %>%
#     mutate(jacc = future_pmap_dbl(list(cat1, cat2), function(cat1, cat2) {
#         genes1 <- dge_all %>% filter(ct_subcluster == cat1)
#         genes2 <- dge_all %>% filter(ct_subcluster == cat2)
#         jaccardBinary(genes1$gene, genes2$gene)
#     }))
# 
# sparse_index_i <- match(dge_jaccard$cat1, subclusters)
# sparse_index_j <- match(dge_jaccard$cat2, subclusters)
# jacc_matrix <- sparseMatrix(i = sparse_index_i, j = sparse_index_j, x = dge_jaccard$jacc)

write_csv(dge_all, file.path(out_dir, "liger_genesets.csv"))
write_tsv(dge_all %>% select(ensembl, subcluster_dx), file.path(out_dir, "liger_modules.tsv"), col_names = FALSE)
write_tsv(dge_all %>% filter(dge_type %in% c("ad_up", "ad_down")) %>% select(ensembl, subcluster_dx),
          file.path(out_dir, "ad_modules.tsv"),
          col_names = FALSE)
#writeMM(jacc_matrix, file.path(out_dir, "clinical_dx_lme_and_lm_genesets_jaccard.mtx"))

dge_counts <- dge_all %>% group_by(subcluster_dx) %>% summarize(n = n(), ad = mean(t_statistic_ad_lme))
write_tsv(dge_counts, file.path(out_dir, "gene_counts.tsv"))
