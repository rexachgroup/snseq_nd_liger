# MDS plot of scenic rss scores and group top 10 terms per dx.

liblist <- c("tidyverse", "ggrepel", "ggalt")
l <- lapply(liblist, library, character.only = TRUE, quietly = TRUE)

SCENIC_TABLES <- "../../resources/scenic/rss"
OUT_DIR <- "../../analysis/seurat_lchen/scenic_mds/"
plot_ts <- function() { str_glue("{date()} {system('md5sum scenic_mds.R', intern = T)}")}

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

    # Compute mds on scenic rss scores.
    tables <- tables %>%
        mutate(regulon_coords = pmap(., function(...) {
            cr <- list(...)
            mds <- as.data.frame(cmdscale(dist(cr$data[2:5]), k = 2)) # AD_[region], bvFTD_[region], Control_[region], PSP-S_[region]
            colnames(mds) <- c("mds1", "mds2")
            return(mds)
        }))
    # Bind original table with mds coordinates.
    tables <- tables %>%
        mutate(plot_tb = pmap(., function(...) {
            cr <- list(...)
            tb <- rename(cr$data, regulon = `...1`)
            tb <- mutate(tb, regulon = str_match(regulon, "^[a-zA-Z0-9]+")[,1]) # only 1st alphanumeric component, drop _exended and (\\d+g)
            tb <- bind_cols(tb, cr$regulon_coords)
            return(tb)
        })

    # Format plot. 
    # Draw one geom_encircle for each of the top 10 regulons per dx.
    # Label top regulons with geom_label_repel.
    tables <- tables %>%
        mutate(plot_obj = pmap(., function(...) {
            cr <- list(...)
            ad_pts <- slice_col_head(cr$plot_tb, str_glue("AD_{cr$celltype}_rank"), 10)
            ftd_pts <- slice_col_head(cr$plot_tb, str_glue("bvFTD_{cr$celltype}_rank"), 10)
            psp_pts <- slice_col_head(cr$plot_tb, str_glue("PSP-S_{cr$celltype}_rank"), 10)
            ctl_pts <- slice_col_head(cr$plot_tb, str_glue("Control_{cr$celltype}_rank"), 10)
            all_pts <- bind_rows(ctl_pts, ad_pts, ftd_pts, psp_pts) %>% 
                filter(!duplicated(regulon))
            
           ggplot(cr$plot_tb, aes(x = mds1, y = mds2, label = regulon)) +
               geom_encircle(data = ctl_pts, aes(color = "Control"), size = 3, s_shape = 1, expand = 0.005) +
               geom_encircle(data = ad_pts, aes(color = "AD"), size = 3, s_shape = 1, expand = 0.005) +
               geom_encircle(data = ftd_pts, aes(color = "FTD"), size = 3, s_shape = 1, expand = 0.005) +
               geom_encircle(data = psp_pts, aes(color = "PSP"), size = 3, s_shape = 1, expand = 0.005) +
               geom_point() +
               geom_label_repel(data = all_pts, fill = "white", max.overlaps = 999) +
               theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
               scale_color_manual(name = "dx", breaks = c("Control", "AD", "FTD", "PSP"), values = c("gray25", "red2", "blue2", "gold2"))+ 
               scale_fill_manual(name = "dx", breaks = c("Control", "AD", "FTD", "PSP"), values = c("gray25", "red2", "blue2", "gold2"))+ 
               labs(title = str_glue("{unique(cr$celltype)} {unique(cr$region)} scenic mds"))
        }))
    
    pdf(file.path(OUT_DIR, "scenic_mds.pdf"), width = 15, height = 15)
    ts <- plot_ts()
    pwalk(tables, function(...) {
        cr <- list(...)
        print(cr$plot_obj + labs(subtitle = ts))
    })
    graphics.off()
}

slice_col_head <- function(tb, col, thresh) {
    tb %>% filter(.data[[col]] < thresh)
}

main()
