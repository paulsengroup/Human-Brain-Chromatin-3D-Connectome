library(tidyr)
library(tidyverse)

groups <- c(
  "GABA_infancy", 
  "GABA_childhood", 
  "GABA_adolescence", 
  "GABA_adulthood", 
  "GABA_lateadulthood", 
  "GLU_infancy", 
  "GLU_childhood", 
  "GLU_adolescence",
  "GLU_adulthood", 
  "GLU_lateadulthood", 
  "MgAs_infancy", 
  "MgAs_childhood", 
  "MgAs_adolescence",
  "MgAs_adulthood", 
  "MgAs_late_adulthood", 
  "Olig_infancy", 
  "Olig_childhood", 
  "Olig_adolescence", 
  "Olig_adulthood", 
  "Olig_late_adulthood"
)

base_path <- "sim_results/"

# Target subcompartment order
target <- c("A3", "A2", "A1", "A0", "B0", "B1", "B2", "B3")

# Number of random samples
number_of_samples <- 100

process_group <- function(group_name) {
  file_path <- file.path(base_path, paste0(group_name, "_distance_to_center_subcom.tsv"))
  
  
  df <- read.csv(file_path, sep = "\t", header = FALSE)
  df <- separate(df, col = V1, into = c("chromosome", "start", "end"))[, c(1, 2, 3, 4, 5)]
  colnames(df)[4:5] <- c("sub_compartment", "Distance from the nucleus center")
  df[is.na(df)] <- "NA"
  
  # Compute median distance for each subcompartment across samples
  group_per_sample <- data.frame()
  step <- floor(nrow(df) / number_of_samples)
  
  for (i in seq(1, nrow(df), by = step)) {
    sample <- df[i:min(i + step - 1, nrow(df)), ]
    sample <- sample[order(sample$`Distance from the nucleus center`, decreasing = FALSE), ]
    sample_sub <- sample %>%
      group_by(sub_compartment) %>%
      summarise(median_pos = median(`Distance from the nucleus center`))
    sample_sub <- sample_sub[match(target, sample_sub$sub_compartment), ]
    group_per_sample <- rbind(group_per_sample, sample_sub)
  }
  
  # Median and SD per subcompartment
  sub_median <- group_per_sample %>%
    group_by(sub_compartment) %>%
    summarise(median_pos = median(median_pos, na.rm = TRUE))
  
  sub_sd <- group_per_sample %>%
    group_by(sub_compartment) %>%
    summarise(sd_pos = sd(median_pos, na.rm = TRUE))
  
  sub_median <- sub_median[match(target, sub_median$sub_compartment), ]
  sub_sd <- sub_sd[match(target, sub_sd$sub_compartment), ]
  
  list(median = sub_median, sd = sub_sd)
}

# Process all samples
results <- lapply(groups, process_group)
names(results) <- groups

# Combine all medians into one dataframe for comparison
median_df <- map_dfr(results, ~.x$median, .id = "group")





df <- median_df %>%
  separate(group, into = c("Cell_type", "age"), sep = "_", extra = "merge")

df$age <- gsub("lateadulthood", "late_adulthood", df$age)

df$sub_compartment <- factor(df$sub_compartment,
                             levels = c("A3","A2","A1","A0","B0","B1","B2","B3"))

pdf(file.path(base_path, "sub_com_radial_pos_all_groups_divided_by_cell_type.pdf"), width = 8, height = 8)

ggplot(df, aes(x = sub_compartment, y = median_pos, 
               color = age, group = age)) +
  geom_point(position = position_dodge(0.3)) +
  geom_line(position = position_dodge(0.3)) +
  facet_wrap(~ Cell_type) +
  theme_bw() +
  labs(title = "Median radial position by sub-compartment",
       y = "Median radial distance")

dev.off()


df$age <- factor(df$age,
                             levels = c("infancy","childhood","adolescence","adulthood","late_adulthood"))

pdf(file.path(base_path, "sub_com_radial_pos_all_groups_divided_by_age.pdf"), width = 12, height = 8)

ggplot(df, aes(x = sub_compartment, y = median_pos, 
               color = Cell_type, group = Cell_type)) +
  geom_point(position = position_dodge(0.3)) +
  geom_line(position = position_dodge(0.3)) +
  facet_wrap(~ age) +
  theme_bw() +
  labs(title = "Median radial position by sub-compartment",
       y = "Median radial distance")

dev.off()

