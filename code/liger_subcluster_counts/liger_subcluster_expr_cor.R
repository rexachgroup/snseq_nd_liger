# Calculate correlation between gene expression values across a subcluster contrast.
liblist <- c("Seurat", "tidyverse", "broom")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

in_seurat_rds <- "../../analysis/pci_import/pci_seurat.rds"
in_liger_meta <- "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"

out_path_base <- "../../analysis/seurat_lchen/liger_subcluster_expr_cor"
batchtools <- file.path(out_path_base, "batchtools")

test_statements <- tribble(
    ~name,                  ~gene1,     ~gene2,     ~subset_quo,
    "insula-exc-2-ftd",     "RORB",     "NPTX2",    quo(ct_subcluster == "insula-excitatory-2" & clinical_dx == "bvFTD"),
    "insula-exc-25-ftd",    "RORB",     "NPTX2",    quo((ct_subcluster %in% c("insula-excitatory-2", "insula-excitatory-5")) & clinical_dx == "bvFTD"),
    "insula-exc-5-ftd",     "RORB",     "NPTX2",    quo(ct_subcluster == "insula-excitatory-5" & clinical_dx == "bvFTD"),
    "precg-exc-4-ad",       "RORB",     "NPTX2",    quo(ct_subcluster == "preCG-excitatory-4" & clinical_dx == "AD"),
    "precg-exc-47-ad",      "RORB",     "NPTX2",    quo((ct_subcluster %in% c("preCG-excitatory-4", "preCG-excitatory-7")) & clinical_dx == "AD"),
    "insula-exc-2-ftd",     "RORB",     "ELK4",     quo(ct_subcluster == "insula-excitatory-2" & clinical_dx == "bvFTD"),
    "insula-exc-25-ftd",    "RORB",     "ELK4",     quo((ct_subcluster %in% c("insula-excitatory-2", "insula-excitatory-5")) & clinical_dx == "bvFTD"),
    "insula-exc-5-ftd",     "RORB",     "ELK4",     quo(ct_subcluster == "insula-excitatory-5" & clinical_dx == "bvFTD"),
    "precg-exc-4-ad",       "RORB",     "ELK4",     quo(ct_subcluster == "preCG-excitatory-4" & clinical_dx == "AD"),
    "precg-exc-47-ad",      "RORB",     "ELK4",     quo((ct_subcluster %in% c("preCG-excitatory-4", "preCG-excitatory-7")) & clinical_dx == "AD"),
    "precg-ins-nptx-comb",  "RORB",     "NPTX2",    quo((ct_subcluster == "insula-excitatory-2" & clinical_dx == "bvFTD") | (ct_subcluster == "preCG-excitatory-4" & clinical_dx == "AD"))
)

main <- function() { 
    dir.create(out_path_base, recursive = TRUE)
    liger_meta <- readRDS(in_liger_meta)
    sobj <- readRDS(in_seurat_rds)

    genes_subset <- unique(c(test_statements$gene1, test_statements$gene2))
    sobj_gene_subset <- subset(sobj, features = genes_subset)

    cor_tb <- mutate(test_statements, t_test_tb = pmap(test_statements, function(...) {
        cr <- list(...)
        gene_cor(sobj_gene_subset, liger_meta, cr$gene1, cr$gene2, cr$subset_quo)
    }))

    saveRDS(cor_tb, file.path(out_path_base, "cor_tb.rds"))

    cor_csv <- cor_tb %>%
        unnest(t_test_tb) %>%
        mutate(subset_quo = as.character(subset_quo)) %>%
        print(width = Inf) %>%
        write_csv(file.path(out_path_base, "cor_tb.csv"))
    
    cor_nz_tb <- mutate(test_statements, t_test_tb = pmap(test_statements, function(...) {
        cr <- list(...)
        gene_cor_nonzero(sobj_gene_subset, liger_meta, cr$gene1, cr$gene2, cr$subset_quo)
    }))

    saveRDS(cor_nz_tb, file.path(out_path_base, "cor_nz_tb.rds"))

    cor_nz_csv <- cor_nz_tb %>%
        unnest(t_test_tb) %>%
        mutate(subset_quo = as.character(subset_quo)) %>%
        print(width = Inf) %>%
        write_csv(file.path(out_path_base, "cor_nonzero_tb.csv"))

    cor_nz_tb <- cor_nz_tb %>% mutate(plot = pmap(., function(...) {
        cr <- list(...)
        title <- str_glue("liger_subcluster_expr_cor.R\n{cr$name} {cr$gene1} {cr$gene2} n={cr$t_test_tb$n} cor={cr$t_test_tb$pearson.cor}")
        cell_cor_plot(sobj_gene_subset, liger_meta, title, cr)
    }))

    pdf(file.path(out_path_base, "cor_nz_tb.pdf"))
    walk(cor_nz_tb$plot, print)
    graphics.off()
}

gene_cor <- function(sobj, liger_meta, gene1, gene2, cell_subset_quo) {
    cell_subset <- liger_meta %>% filter(!!cell_subset_quo)
    sobj_subset <- subset(sobj, cells = cell_subset$UMI, features = c(gene1, gene2))
    sobj_mat <- GetAssayData(sobj_subset)

    vec1 <- sobj_mat[gene1, ]
    vec2 <- sobj_mat[gene2, ]
    vec_cor <- cor(vec1, vec2, method = "pearson")
    vec_t_test <- tidy(t.test(vec1, vec2))

    res_tbl <- tibble(n = length(vec1), pearson.cor = vec_cor, vec_t_test = vec_t_test) %>% unnest(vec_t_test)
    return(res_tbl) 
}

gene_cor_nonzero <- function(sobj, liger_meta, gene1, gene2, cell_subset_quo) {
    cell_subset <- liger_meta %>% filter(!!cell_subset_quo)
    sobj_subset <- subset(sobj, cells = cell_subset$UMI, features = c(gene1, gene2))
    sobj_mat <- GetAssayData(sobj_subset)

    vec1 <- sobj_mat[gene1, ]
    vec2 <- sobj_mat[gene2, ]
    nonzeroes <- (vec1 > 0 & vec2 > 0)
    vec1 <- vec1[nonzeroes]
    vec2 <- vec2[nonzeroes]
    vec_cor <- cor(vec1, vec2, method = "pearson")
    vec_t_test <- tidy(t.test(vec1, vec2))

    res_tbl <- tibble(n = length(vec1), pearson.cor = vec_cor, vec_t_test = vec_t_test) %>% unnest(vec_t_test)
    return(res_tbl) 
}

cell_cor_plot <- function(sobj, liger_meta, title, cr) {
    cell_subset <- liger_meta %>% filter(!!cr$subset_quo)
    sobj_subset <- subset(sobj, cells = cell_subset$UMI, features = c(cr$gene1, cr$gene2))
    sobj_mat <- GetAssayData(sobj_subset)
    
    vec1 <- sobj_mat[cr$gene1, ]
    vec2 <- sobj_mat[cr$gene2, ]
    nonzeroes <- (vec1 > 0 & vec2 > 0)
    vec1 <- vec1[nonzeroes]
    vec2 <- vec2[nonzeroes]

    plot_tb <- tibble(x = vec1, y = vec2)

    gg <- ggplot(plot_tb, aes(x = x, y = y)) +
        geom_point() +
        geom_abline(slope = cr$t_test_tb$pearson.cor, intercept = cr$t_test_tb$estimate) +
        ggtitle(title) +
        labs(x = "gene1", y = "gene2")
    return(gg)
}

if (!interactive()) main()
