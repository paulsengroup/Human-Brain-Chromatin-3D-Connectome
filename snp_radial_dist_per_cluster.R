library(ggpubr)

input_dir <- ""

files <- list.files(input_dir, pattern = ".tsv$", full.names = TRUE)



for (file in files) {
  base_name <- tools::file_path_sans_ext(basename(file))
  
  coordinates <- read.csv(
    paste0(input_dir, base_name, "_beads.tsv"),
    sep = "\t", check.names = FALSE
  )
  coordinates$dist = sqrt(coordinates$x_1^2 + coordinates$y_1^2 + coordinates$z_1^2)
  coordinates = coordinates[, c(2, 8)]
  clusters <- read.csv(
    file,
    sep = "\t", check.names = FALSE
  )
  colnames(clusters)[1] = "region"
  combined = merge(coordinates, clusters, by = "region")
  
  pdf(paste0(input_dir, base_name, "_distance_from_center.pdf"))
  p = ggplot(combined, aes(x = factor(cluster), y = dist)) +
    geom_boxplot(fill = "white") +
    labs(x = "Cluster",
         y = "Distance from nuclear center") +
    theme_minimal()
  
  print(p)
  dev.off()
  
}
