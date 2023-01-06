# MDS plot of scenic rss scores + draw ellipses around top 10 terms per dx

liblist <- c("tidyverse", "ggrepel", "ggalt")
l <- lapply(liblist, library, character.only = TRUE, quietly = TRUE)

SCENIC_TABLES <- "../../resources/scenic/rss"
OUT_DIR <- "../../analysis/seurat_lchen/scenic_mds/"

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
        mutate(regulon_coords = map(data, function(tb) {
            mds <- as.data.frame(cmdscale(dist(tb[2:5]), k = 2)) # AD_[region], bvFTD_[region], Control_[region], PSP-S_[region]
            colnames(mds) <- c("mds1", "mds2")
            return(mds)
        }))
    # bind original table with mds results
    tables <- tables %>%
        mutate(plot_tb = pmap(., function(...) {
            cr <- list(...)
            tb <- rename(cr$data, regulon = `...1`)
            tb <- mutate(tb, regulon = str_match(regulon, "^[a-zA-Z0-9]+")[,1]) # only 1st alphanumeric group, drop _exended and (\\d+ g)
            tb <- bind_cols(tb, cr$regulon_coords)
            return(tb)
        }))

    #     tables <- tables %>%
    #         mutate(plot_obj = pmap(., function(...) {
    #             cr <- list(...)
    #             non_grouped <- cr$plot_tb %>% filter(color == "black")
    #             grouped <- cr$plot_tb %>% filter(color != "black")
    #             ggplot(cr$plot_tb, aes(x = mds1, y = mds2, colour = color, fill = color, label = regulon)) +
    #                 scale_color_manual(values = c("AD_top10" = "red", "FTD_top10" = "green", "PSP_top10" = "blue", "CTL_top10" = "purple", "black" = "black")) +
    #                 scale_fill_manual(values = c("AD_top10" = "red", "FTD_top10" = "green", "PSP_top10" = "blue", "CTL_top10" = "purple", "black" = "black")) +
    #                 geom_point() + 
    #                 geom_label_repel(fill = "white") +
    #                 stat_ellipse(data = grouped)
    #         }))

    tables <- tables %>%
        mutate(plot_obj = pmap(., function(...) {
            cr <- list(...)
            ad_pts <- slice_col_head(cr$plot_tb, str_glue("AD_{cr$celltype}_rank"), 10)
            ftd_pts <- slice_col_head(cr$plot_tb, str_glue("bvFTD_{cr$celltype}_rank"), 10)
            psp_pts <- slice_col_head(cr$plot_tb, str_glue("PSP-S_{cr$celltype}_rank"), 10)
            ctl_pts <- slice_col_head(cr$plot_tb, str_glue("Control_{cr$celltype}_rank"), 10)
            
           ggplot(cr$plot_tb, aes(x = mds1, y = mds2, label = regulon)) +
               geom_point() +
               geom_label_repel(fill = "white") +
               geom_encircle(data = ad_pts, color = "red", s_shape = 1, expand = 0) +
               geom_encircle(data = ftd_pts, color = "green", s_shape = 1, expand = 0) +
               geom_encircle(data = psp_pts, color = "blue", s_shape = 1, expand = 0) +
               geom_encircle(data = ctl_pts, color = "purple", s_shape = 1, expand = 0) +
               theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        }))
    
    pwalk(tables, function(...) {
        cr <- list(...)
        pdf(cr$out_path, width = 20, height = 20)
        print(cr$plot_obj)
        dev.off()
    })
}

slice_col_head <- function(tb, col, thresh) {
    tb %>% filter(.data[[col]] < thresh)
}

top10_rank_color_row <- function(...) {
    tb_row <- tibble(...)
    threshold <- 15
    ad_col <- str_subset(colnames(tb_row), "AD_.*_rank")
    ftd_col <- str_subset(colnames(tb_row), "bvFTD_.*_rank")
    psp_col <- str_subset(colnames(tb_row), "PSP-S_.*_rank")
    ctl_col <- str_subset(colnames(tb_row), "Control_.*_rank")
    if (tb_row[[ad_col]] < threshold) {
        return("AD_top10")
    } else if (tb_row[[ftd_col]] < threshold) {
        return("FTD_top10")
    } else if (tb_row[[psp_col]] < threshold) {
        return("PSP_top10")
    } else if (tb_row[[ctl_col]] < threshold) {
        return("CTL_top10")
    } else {
        return("black")
    }
}

main()
