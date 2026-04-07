library(scatterpie)

input_dir <- ""
files <- list.files(input_dir, pattern = "*.tsv$", full.names = TRUE)

all_results <- list()

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
  
  trait_map <- tibble::tibble(
    trait_2 = c("AD","PD","MS","IBD"),
    category = c(
      rep("neurodegenerative", 2),
      rep("immune", 2)
    )
  )
  
  category_size <- tibble::tibble(
    category = c("neurodegenerative", "immune"),
    total_size = c(236, 392)
  )
  
  category_norm <- cluster_bead %>%
    separate_rows(trait_2, sep = ",") %>%
    left_join(trait_map, by = "trait_2") %>%
    filter(!is.na(category)) %>%
    dplyr::count(cluster, category) %>%
    left_join(category_size, by = "category") %>%
    mutate(norm = n / total_size) %>%  
    group_by(cluster) %>%
    mutate(norm_share = norm / sum(norm))
  
  result <- category_norm %>%
    group_by(cluster) %>%
    summarise(
      assigned_category = if (any(norm_share > 0.6666667)) {
        category[which.max(norm_share)]
      } else {
        "immune/neurodegenerative"
      },
      .groups = "drop"
    )
  
  GO = read.xlsx(paste0(input_dir, "GO/GO_",base_name, ".xlsx"))
  
  result <- result %>%
    mutate(cluster = as.character(cluster))

  GO <- GO %>%
    dplyr::left_join(result, by = c("Cluster" = "cluster")) %>%
    dplyr::mutate(Cluster = assigned_category) %>%
    dplyr::select(-assigned_category)  %>%
    group_by(Cluster) %>%
    slice_min(order_by = p.adjust, n = 2, with_ties = FALSE) %>%
    ungroup()
  
  GO$sample = base_name
  
  all_results[[base_name]] <- GO

}

final_df <- bind_rows(all_results)


plot_df <- final_df %>%
  mutate(value = 1) %>% 
  group_by(Cluster, Description, sample) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = value, values_fill = 0)

plot_df <- plot_df %>%
  mutate(
    x = as.numeric(factor(Cluster,
                          levels = c("immune", "neurodegenerative", "immune/neurodegenerative")
    )) * 4,  
    y = as.numeric(factor(Description))
  )


pdf("GO_summary.pdf", width = 10, height = 10)
ggplot() +
  geom_scatterpie(
    data = plot_df,
    aes(x = x, y = y),
    cols = unique(final_df$sample),
    pie_scale = 2.5
  ) +
  scale_x_continuous(
    breaks = c(4, 8, 12),
    labels = c("immune", "neurodegenerative", "immune/neurodegenerative")
  ) +
  scale_y_continuous(
    breaks = plot_df$y,
    labels = stringr::str_wrap(plot_df$Description, width = 40)
  ) +
  theme_minimal() +
  labs(x = "Cluster", y = "Description") + coord_fixed() +
  scale_fill_manual(
    values = scales::hue_pal()(4),
    labels = c(
      "GABA_adolescence" = "GABA",
      "GLU_adolescence" = "GLU",
      "MgAs_adolescence" = "MgAs",
      "Olig_adolescence" = "Olig"
    )
  )

dev.off()
