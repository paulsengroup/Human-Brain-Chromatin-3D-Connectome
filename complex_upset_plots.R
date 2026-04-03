library(tidyverse)
library(dplyr)
library(ComplexUpset)

input_dir = "dir"
df1 <- read.delim(paste0(input_dir, "GABA_adolescence_25_top30_beads.tsv"), stringsAsFactors = FALSE)
df2 <- read.delim(paste0(input_dir, "GLU_adolescence_6_top30_beads.tsv"), stringsAsFactors = FALSE)
df3 <- read.delim(paste0(input_dir, "Olig_adolescence_100_top30_beads.tsv"), stringsAsFactors = FALSE)
df4 <- read.delim(paste0(input_dir, "MgAs_adolescence_40_top30_beads.tsv"), stringsAsFactors = FALSE)

df1$cell_type <- "GABA"
df2$cell_type <- "GLU"
df3$cell_type <- "Olig"
df4$cell_type <- "MgAs"

df <- bind_rows(df1, df2, df3, df4)

df_long <- df %>%
  separate_rows(trait_2, sep = ",") %>%
  mutate(trait_2 = str_trim(trait_2)) %>%
  filter(trait_2 %in% traits) %>%
  distinct(region, cell_type, trait_2)



df_bin <- df_long %>%
  mutate(value = TRUE) %>%
  pivot_wider(
    names_from = trait_2,
    values_from = value,
    values_fill = FALSE
  )

pdf(paste0(input_dir, "snp_stats_per_cell_type.pdf"))
upset(
  df_bin,
  intersect = traits,
  name = "Traits",
  
  base_annotations = list(
    "Intersection size" = intersection_size()
  ),
  
  annotations = list(
    "cell_type composition" = (
      ggplot(mapping = aes(fill = cell_type)) +
        geom_bar(position = "fill") +
        ylab("Proportion")
    )
  ),
  set_sizes = FALSE 
)
dev.off()


