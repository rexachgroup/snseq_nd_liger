# Damon Polioudakis
# 2018-09-10
# Function library

# Cowplot plot_grid and add title
plot_grid_wrapper <- function(plotlist, ncol = 2, title = "", rel_height = 0.1
  , ...) {
  # cowplot plot_grid ...: align = 'v', axis = 'l'
  # Plot grid
  pg <- plot_grid(plotlist = plotlist, ncol = ncol, ...)
  # now add the title
  title <- ggdraw() + draw_label(title)
  # rel_heights values control title margins
  plot_grid(title, pg, ncol = 1, rel_heights = c(rel_height, 1))
}

clean_strings <- function(string_vector){
  print("clean_strings")
  cleaned <- string_vector %>%
    gsub("* ", "_", .) %>%
    gsub("\\.", "_", .) %>%
    gsub("\\(", "_", .) %>%
    gsub("\\)", "_", .) %>%
    gsub("\\+", "and", .) %>%
    gsub("#", "_number", .) %>%
    gsub("_$", "", .) %>%
    gsub("__", "_", .) %>%
    tolower
  return(cleaned)
}

clean_variable_names <- function(data){
  cleaned <- data %>%
    rename_all(
      funs(clean_strings)
    )
  return(cleaned)
}

make_plot_title <- function(title){
  paste0(script_name
    , "\n\n", title
    , "\n", graph_subtitle)
}
