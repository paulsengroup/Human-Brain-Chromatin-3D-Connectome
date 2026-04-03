library(pheatmap)
library(dplyr)
library(cluster)
library(RColorBrewer)


input_dir <- "dir"
files <- list.files(input_dir, pattern = "top_domains.tsv$", full.names = TRUE)



for (file in files) {
  base_name <- tools::file_path_sans_ext(basename(file))
  cluster_file <- file.path(input_dir, paste0(base_name, ".tsv"))
  heatmap_file <- file.path(input_dir, paste0(base_name, ".pdf"))
  
  distance_matrix = read.csv(file, 
                             sep = "\t", check.names = FALSE)
  
  rownames(distance_matrix) = colnames(distance_matrix)
  
  
  hc <- hclust(as.dist(distance_matrix), method = "complete")
  ph <- pheatmap(distance_matrix, cluster_rows = hc, cluster_cols = hc, silent = TRUE)
  
  sil_width <- c()
  
  for (k in 2:20) {
    if (k >= length(rownames(distance_matrix))) break
    cl <- cutree(ph$tree_row, k = k)
    if (any(table(cl) < 2)) next
    sil <- silhouette(cl, distance_matrix)
    sil_width[k] <- mean(sil[, 3])
  }
  
  
  if (all(is.na(sil_width))) {
    best_k <- 1
  } else {
    best_k <- which.max(sil_width)
  }
  
  cl <- cutree(ph$tree_row, k = best_k)
  cluster_df <- tibble(
    bead_id = names(cl),
    cluster = as.integer(cl)
  )
  
  
  annotation_row <- as.data.frame(cluster_df)
  rownames(annotation_row) <- annotation_row$bead_id
  annotation_row$bead_id <- NULL
  annotation_row <- annotation_row[rownames(distance_matrix), , drop = FALSE]
  annotation_row$cluster <- factor(annotation_row$cluster)
  n_clusters <- nlevels(annotation_row$cluster)
  
  
  distance_matrix_capped <- distance_matrix
  distance_matrix_capped[distance_matrix_capped > 7] <- 7
  
  breaks <- seq(0, 7, length.out = 101)
  colors <- colorRampPalette(c("red", "white", "blue"))(100)
  
  
  
  if (n_clusters == 1) {
    cluster_palette <- "steelblue"
  } else {
    cluster_palette <- grDevices::hcl.colors(n_clusters, "Dynamic")
  }
  
  cluster_colors <- list(
    cluster = setNames(cluster_palette, levels(annotation_row$cluster))
  )
  
  
  pdf(heatmap_file, width = 12, height = 12)
  pheatmap(
    distance_matrix_capped,
    cluster_rows = hc,
    cluster_cols = hc,
    color = colors,
    breaks = breaks,
    annotation_row = annotation_row,
    annotation_col = annotation_row,
    show_rownames = FALSE,
    show_colnames = FALSE,
    annotation_colors = cluster_colors,
    na_col = "grey80"
  )
  dev.off()
  
  write.table(cluster_df,
              cluster_file,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
}
