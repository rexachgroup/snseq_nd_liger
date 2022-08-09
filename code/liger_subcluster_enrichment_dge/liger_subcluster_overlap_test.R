# Fishers' exact test with external geneset.

set.seed(0)
liblist <- c("Seurat", "tidyverse", "broom")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

OUT_DIR <- "../../analysis/seurat_lchen/liger_subcluster_enrichment_dge/overlap/"
DGE_RDS <- "../../analysis/seurat_lchen/liger_subcluster_enrichment_dge/subcluster_wk.rds"
GENELIST <- "../../resources/Cravatt_broadTarget_list.csv"

main <- function() {
    dir.create(OUT_DIR)
    subcluster_wk <- readRDS(DGE_RDS)
    genelist <- read_csv(GENELIST)

    genelist_grp <- genelist %>%
        rename(reference_group = group) %>%
        group_by(reference_group) %>%
        group_nest(.key = "ref_genes")
    
    subcluster_wk <- subcluster_wk %>%
        filter(!is.na(broom_join))
    subcluster_wk <- mutate(subcluster_wk, broom_gene_list = map(broom_join, fmt_broom_tb, filter = T))

    cluster_gene_bg <- subcluster_wk %>%
        group_by(region, cluster_cell_type) %>%
        group_nest %>%
        mutate(bg_gene_list = map(data, function(dat) {
            bg_gene_lists <- map(dat$broom_join, fmt_broom_tb, filter = F) %>%
                map(pluck("gene"))
            Reduce(union, bg_gene_lists)
        })) %>%
        select(region, cluster_cell_type, bg_gene_list)
    
    subcluster_test <- subcluster_wk %>%
        left_join(cluster_gene_bg, by = c("region", "cluster_cell_type")) %>%
        left_join(genelist_grp, by = character())

    subcluster_test <- subcluster_test %>%
        mutate(fisher_dat = pmap(., function(...) {
            cr <- list(...)
            writeLines(str_glue("{cr$subcluster_ct}"))
            ora(cr$broom_gene_list$gene, cr$bg_gene_list, cr$ref_genes$gene)
        }))

    subcluster_padj <- subcluster_test %>%
        group_by(region, cluster_cell_type, reference_group) %>%
        group_nest() %>%
        mutate(pval_adj = map(data, function(dat) {
            dat %>%
                select(ct_subcluster, fisher_dat) %>%
                unnest(fisher_dat) %>%
                mutate(padj = p.adjust(p.value, method = "BH"))
        }))
        
    subcluster_padj %>% 
        select(-data) %>%
        saveRDS(file.path(OUT_DIR, "subcluster_overlap.rds"))

    subcluster_padj %>%
        select(-data) %>%
        unnest(pval_adj) %>%
        print(width = Inf, n = Inf) %>%
        write_csv(file.path(OUT_DIR, "subcluster_overlap.csv"))
}

fmt_broom_tb <- function(tb, filter) {
    tb <- tb %>%
        select(gene, estimate = tidyselect::contains("estimate"), padj = tidyselect::contains("p.value.adj"))
    if (filter) {
        tb <- tb %>%
            filter(!is.na(estimate) & !is.na(padj)) %>%
            filter(abs(estimate) > 0.1, padj < 0.05)
    }
    return(tb)
}

ora <- function(test, bg, reference) { 
    q <- length(intersect(test,reference))
    k <- length(intersect(reference,bg))
    m <- length(intersect(test,bg))
    t <- length(bg)
    writeLines(str_glue("{q} {k} {m} {t}"))
    fisher_obj <- fisher.test(matrix(c(q, k-q, m-q, t-m-k+q), 2, 2),conf.int=TRUE)
    return(broom::tidy(fisher_obj))
}

if (!interactive()) main()
