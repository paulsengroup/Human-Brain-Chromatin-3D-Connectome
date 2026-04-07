library(stringr)

input_dir <- ""
files <- list.files(input_dir, pattern = ".tsv$", full.names = TRUE)

for (file in files) {
  base_name <- tools::file_path_sans_ext(basename(file))
  cluster_file <- file.path(input_dir, paste0(base_name, ".tsv"))
  clusters = read.csv(cluster_file, 
                      sep = "\t", check.names = FALSE)
  
  bead_file <- file.path(input_dir, paste0(base_name, "_beads.tsv"))
  beads = read.csv(bead_file, 
                   sep = "\t", check.names = FALSE)
  
  colnames(clusters)[1] = "region"
  cluster_bead = merge(clusters, beads, by="region")[, c(1, 2, 3)]
  
  
  category_size <- tibble::tibble(
    trait_2 = c("AD", "PD", "MS", "IBD"),
    total_size = c(84, 152, 200, 192)
  )
  
  category_norm <- cluster_bead %>%
    separate_rows(trait_2, sep = ",") %>%
    filter(trait_2 %in% c("AD", "PD", "MS", "IBD")) %>%    
    dplyr::count(cluster, trait_2) %>%
    left_join(category_size, by = "trait_2") %>%
    mutate(norm = n / total_size) %>%        # normalize by gwas size
    group_by(cluster) %>%
    mutate(norm_share = norm / sum(norm))

  trait_order = c("IBD", "MS", "AD", "PD")
  category_norm$trait_2 <- factor(category_norm$trait_2, levels = trait_order)
  
  pdf(file.path(input_dir, paste0("trait_compos/", base_name, "_gwas_share.pdf")), width = 8, height = 8)
  
  print(ggplot(category_norm, aes(x = factor(cluster), y = norm_share, fill = trait_2)) +
          geom_col() +
          labs(x = "Cluster", y = "Share", fill = "GWAS") +
          theme_minimal()
  )
  dev.off()

}
