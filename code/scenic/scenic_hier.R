# Hierarchical clustering + heatmap of scenic regulons.

liblist <- c("tidyverse", "ggplot2", "ComplexHeatmap", "circlize", "scales", "patchwork")
l <- lapply(liblist, library, character.only = TRUE, quietly = TRUE)

SCENIC_TABLES <- "../../resources/scenic/rss"
OUT_DIR <- "../../analysis/seurat_lchen/scenic_hier/"
plot_ts <- function() { str_glue("{date()} {system('md5sum scenic_hier.R', intern = T)}")}
k_display <- 25

main <- function() {
    dir.create(OUT_DIR)
    tables <- tibble(path = list.files(SCENIC_TABLES, full.names = T))
    tables <- tables %>% mutate(
        data = map(path, read_csv),
        basename = str_match(basename(path), "(regulon_specificity_score)(.*)(.csv)")[,3],
        celltype = str_match(basename, "_*(.+)_(.+)")[,2],
        region = str_match(basename, "_*(.+)_(.+)")[,3]
    )
    
    tables <- tables %>%
        mutate(regulon_data = pmap(., function(...) { 
            cr <- list(...)
            tb <- rename(cr$data, regulon = `...1`)
            tb <- mutate(tb, regulon = str_match(regulon, "^[a-zA-Z0-9]+")[,1]) # only 1st alphanumeric component, drop _extended and (\\d+g)
            
            # subset to rss rank columns and convert to matrix
            ad_col <- str_glue("AD_{cr$celltype}_rank")
            ftd_col <- str_glue("bvFTD_{cr$celltype}_rank") 
            psp_col <- str_glue("PSP-S_{cr$celltype}_rank") 
            control_col <- str_glue("Control_{cr$celltype}_rank")

            # sort by the minimum rank of any dx, then subset to first k_display rows
            tb <- select(tb, regulon, all_of(c(ad_col, ftd_col, psp_col, control_col)))
            tb$min <- pmin(tb[[ad_col]], tb[[ftd_col]], tb[[psp_col]], tb[[control_col]])
            tb <- tb %>% arrange(min) %>% slice_head(n = k_display) %>% select(-min)

            tb_mat <- tb %>%
                as.data.frame %>% 
                column_to_rownames("regulon") %>%
                as.matrix

            return(tb_mat)
        }))

    pdf(file.path(OUT_DIR, "heatmap.pdf"), height = 10, width = 10) 
    subtitle <- plot_ts()
    pwalk(tables, function(...) {
        cr <- list(...)
        hier <- hclust(dist(cr$regulon_data))

        colormap <- circlize::colorRamp2(
            breaks = quantile(cr$regulon_data, seq(0, 1, by = 0.1)),
            viridisLite::viridis(11, direction = -1)
        )
        tmp_pdf <- function(w, h) { pdf(tempfile(), w, h) }
        heat <- Heatmap(
            cr$regulon_data,
            cluster_rows = hier,
            cluster_columns = F,
            name = "rss rank",
            col = colormap
        )
        heatmap_gtree <- grid.grabExpr(
            draw(heat),
            wrap = TRUE,
            device = tmp_pdf,
            width = unit(20, "cm"),
            height = unit(20, "cm")
        )
        title <- str_glue("{cr$basename} regulon rss rank clustering, top {k_display} regulon")
        print(wrap_plots(heatmap_gtree) + plot_annotation(title = title, subtitle = subtitle))

    })
    graphics.off()
}

main()
