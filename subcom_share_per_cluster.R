input_dir1 = ""
input_dir2 = ""

sub_com <- read.delim(paste0(input_dir2, "model_chrom3d_subcompartments.tsv"), stringsAsFactors = FALSE, header = FALSE)
cluster <- read.delim(paste0(input_dir1, "model_top30.tsv"), stringsAsFactors = FALSE)

colnames(sub_com) = c("bead_id", "sub_com")
combined = merge(cluster, sub_com, by = "bead_id")
combined$cluster <- factor(combined$cluster, levels = sort(unique(combined$cluster)))


colors_8 <- rev(c('#3b4cc0', '#5f7fe8', '#88abfd', '#c7d7f0',
              '#f2cbb7', '#d1493f', '#ca3b37', '#b40426'))

sub_com_order <- rev(c('B3','B2','B1','B0','A0','A1','A2','A3'))

combined$sub_com <- factor(combined$sub_com, levels = sub_com_order)

df_plot <- combined %>%
  dplyr::count(cluster, sub_com) %>%
  dplyr::group_by(cluster) %>%
  dplyr::mutate(prop = n / sum(n))

pdf(paste0(input_dir1, "model_sub_com_share.pdf"))
ggplot(df_plot, aes(x = cluster, y = prop, fill = sub_com)) +
  geom_col() +
  ylab("Proportion") +
  scale_fill_manual(values = colors_8, breaks = sub_com_order) +
  labs(fill = "sub_com") +
  theme_minimal()
dev.off()

