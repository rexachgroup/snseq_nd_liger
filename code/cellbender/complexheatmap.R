require(ComplexHeatmap)
require(tidyr)
require(tidyselect)
require(dplyr)
require(purrr)

ch_to_grob <- function(ch, ...) {
    tmp_pdf <- function(w, h) { pdf(tempfile(), w, h) }
    pdf(tempfile())
    heatmap_gtree <- grid.grabExpr(
        draw(ch, ...),
        wrap = TRUE,
        device = tmp_pdf,
        width = unit(10, "cm"),
        height = unit(10, "cm")
    )
    dev.off()
    return(heatmap_gtree)
}

symnum_signif <- function(pval_fdr) {
    symnum(pval_fdr, cutpoints = c(0, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", " "))
}

pivot_matrix <- function(tb, cols_from, values_from, rows_from, fill_na = NA) {
    tb_pivot <- tb %>%
        select(all_of(c(cols_from, values_from, rows_from))) %>%
        pivot_wider(names_from = all_of(cols_from), values_from = all_of(values_from), values_fill = fill_na) %>%
        ungroup

    tb_matrix <- select(tb_pivot, -all_of(rows_from)) %>% as.matrix()

    rownames(tb_matrix) <- tb_pivot %>% 
        rowwise() %>%
        summarize(join_rowname = paste(c_across(all_of(rows_from)), collapse = "|"), .groups = "drop") %>%
        pluck("join_rowname")

    return(tb_matrix)
}

ch_layout_size <- function(ch, unit = "inch") {
    pdf(tempfile()) 
    ch_draw <- draw(ch)
    w <-  ComplexHeatmap:::width(ch_draw)
    w <- convertX(w, unit, valueOnly = TRUE)
    h <- ComplexHeatmap:::height(ch_draw)
    h <- convertY(h, unit, valueOnly = TRUE)
    return(c(w, h))
}

