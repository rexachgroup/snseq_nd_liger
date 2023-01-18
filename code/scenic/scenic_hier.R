# Hierarchical clusteirng + heatmap of scenic regulons.

liblist <- c("tidyverse", "ggplot2", "ComplexHeatmap", "circlize", "scales", "patchwork")
l <- lapply(liblist, library, character.only = TRUE, quietly = TRUE)

SCENIC_TABLES <- "../../resources/scenic/rss"
OUT_DIR <- "../../analysis/seurat_lchen/scenic_hier/"
plot_ts <- function() { str_glue("{date()} {system('md5sum scenic_hier.R', intern = T)}")}

main <- function() {
    dir.create(OUT_DIR)
    tables <- tibble(path = list.files(SCENIC_TABLES, full.names = T))
    tables <- tables %>% mutate(
        data = map(path, read_csv),
        basename = str_match(basename(path), "(regulon_specificity_score)(.*)(.csv)")[,3],
        celltype = str_match(basename, "(.+)_(.+)")[,2],
        region = str_match(basename, "(.+)_(.+)")[,3],
        out_path = str_glue("{OUT_DIR}/{basename}_mds_grouped.pdf")
    )
    
    tables <- tables %>%
        mutate(regulon_data = pmap(., function(...) { 
            cr <- list(...)
            #tb <- cr$data[2:5]# AD_[region], bvFTD_[region], Control_[region], PSP-S_[region]
            tb <- rename(cr$data, regulon = `...1`)
            tb <- mutate(tb, regulon = str_match(regulon, "^[a-zA-Z0-9]+")[,1]) # only 1st alphanumeric component, drop _exended and (\\d+g)
            ad_col <- str_glue("AD_{cr$celltype}_rank")
            ftd_col <- str_glue("bvFTD_{cr$celltype}_rank") 
            psp_col <- str_glue("PSP-S_{cr$celltype}_rank") 
            control_col <- str_glue("Control_{cr$celltype}_rank")

            tb <- select(tb, regulon, ad_col, ftd_col, psp_col, control_col)

            # filter to top 20 per dx.
            ad_pts <- slice_col_head(tb, ad_col, 20)
            ftd_pts <- slice_col_head(tb, ftd_col, 20)
            psp_pts <- slice_col_head(tb, psp_col,20)
            ctl_pts <- slice_col_head(tb, control_col, 20)
            all_pts <- bind_rows(ctl_pts, ad_pts, ftd_pts, psp_pts) %>% 
                filter(!duplicated(regulon))

            all_pts %>%
                as.data.frame %>% 
                column_to_rownames("regulon") %>%
                as.matrix
        }))

    tables <- tables %>%
        mutate(regulon_hier = pmap(., function(...) {
            cr <- list(...)
            regulon_dist <- dist(cr$regulon_data)            
            hier <- hclust(regulon_dist)
            return(hier)
        }))

    pdf(file.path(OUT_DIR, "heatmap.pdf"), height = 20, width = 10) 
    subtitle <- plot_ts()
    pwalk(tables, function(...) {
        cr <- list(...)
        colormap <- circlize::colorRamp2(
            breaks = quantile(cr$regulon_data, seq(0, 1, by = 0.1)),
            viridisLite::viridis(11, direction = -1)
        )
        tmp_pdf <- function(w, h) { pdf(tempfile(), w, h) }
        heat <- Heatmap(
            cr$regulon_data,
            cluster_rows = cr$regulon_hier,
            cluster_columns = F,
            name = "regulon rank",
            col = colormap
        )
        heatmap_gtree <- grid.grabExpr(
            draw(heat),
            wrap = TRUE,
            device = tmp_pdf,
            width = unit(20, "cm"),
            height = unit(20, "cm")
        )
        title <- str_glue("{cr$basename} regulon rank clustering, top 20 each dx")
        print(wrap_plots(heatmap_gtree) + plot_annotation(title = title, subtitle = subtitle))

    })
    graphics.off()

}

slice_col_head <- function(tb, col, thresh) {
    tb %>% filter(.data[[col]] < thresh)
}

main()
