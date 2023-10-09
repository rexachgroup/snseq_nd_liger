liblist <- c("Seurat", "tidyverse", "WGCNA", "broom.mixed")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

META <- "../../analysis/seurat_lchen/"
SUBCLUSTER_UNFILTERED <- "../../analysis/seurat_lchen/subcluster_excluded_dge/seurat_filter/00_unfiltered_meta.csv"
SUBCLUSTER_FILTERED <- "../../analysis/seurat_lchen/subcluster_excluded_dge/seurat_filter/01_subcluster_filtered_meta.csv"
SUBCLUSTER_CB_F <- "../../analysis/seurat_lchen/subcluster_excluded_dge/seurat_filter/02_subcluster_cellbender_filtered_meta.csv"
OUT_DIR <- "../../analysis/seurat_lchen/subcluster_excluded_dge/cluster_composition_test/"

main <- function() {
    dir.create(OUT_DIR)
    s_un <- read_csv(SUBCLUSTER_UNFILTERED)
    s_un <- s_un %>%
        mutate(clinical_dx = fct_relevel(clinical_dx, "Control"))
    library_match <- s_un %>%
        group_by(cluster_cell_type, library_id) %>%
        summarize(frac_matching = sum(cell_type == cluster_cell_type) / n()) %>%
        ungroup() %>%
        write_csv(file.path(OUT_DIR, "library_composition_matching.csv"))

    dx_tb <- s_un %>%
        select(library_id, clinical_dx, region) %>%
        group_by(library_id) %>% slice_head(n = 1)

    library_dx <- left_join(library_match, dx_tb, by = c("library_id"))

    test_contrasts <- expand_grid(
        region = unique(s_un$region),
        cluster_cell_type = unique(s_un$cluster_cell_type)
    )

    lm_formula <- "frac_matching ~ clinical_dx"
    lm_tests <- test_contrasts %>%
        mutate(
            comp_results = pmap(., function(...) {
                cr <- list(...)
                matching_f <- library_dx %>% 
                    filter(cluster_cell_type == cr$cluster_cell_type, region == cr$region)
                tryCatch(broom::tidy(lm(as.formula(lm_formula), data = matching_f)), error = print)
            })
        )

    lm_res <- lm_tests %>%
        filter(unlist(map(comp_results, is.data.frame))) %>%
        unnest(comp_results) %>%
        mutate(formula = lm_formula)%>%
        write_csv(file.path(OUT_DIR, "library_composition_lm.csv"))
    
    test_contrasts <- expand_grid(
        region = unique(s_un$region),
        dx = c("AD", "PSP-S", "bvFTD"),
        cluster_cell_type = unique(s_un$cluster_cell_type)
    )


    cor_tests <- test_contrasts %>%
        mutate(
            comp_results = pmap(., function(...) {
                cr <- list(...)
                #matching_f <- library_dx %>% 
                matching_ctl <- library_dx %>% 
                    filter(clinical_dx == "Control", cluster_cell_type == cr$cluster_cell_type, region == cr$region)
                matching_dx <- library_dx %>% 
                    filter(clinical_dx == cr$dx, cluster_cell_type == cr$cluster_cell_type, region == cr$region)
                matching_f <- bind_rows(matching_ctl, matching_dx)
                tryCatch(
                    broom::tidy(cor.test(as.numeric(matching_f$clinical_dx == cr$dx), matching_f$frac_matching)), 
                error = print)
            })
        )

    cor_res <- cor_tests %>%
        filter(unlist(map(comp_results, ~!"error" %in% class(.x)))) %>%
        #mutate(comp_results = map(comp_results, as.data.frame)) %>%
        unnest(comp_results) %>%
        print(width = Inf, n = Inf) %>%
        write_csv(file.path(OUT_DIR, "library_composition_cor.csv"))
    
}

if (!interactive()) main()
