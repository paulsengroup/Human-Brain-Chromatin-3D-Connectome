#Rscript plot_loss_combined.R adolescence adulthood childhood infancy lateadulthood


args <- commandArgs(trailingOnly = TRUE)
base_names <- args
base_path <- getwd()

n_samples <- length(base_names)

colors <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00")

all_time <- list()
all_mean_loss <- list()

global_min <- Inf
global_max <- -Inf

for (s in seq_along(base_names)) {

  base_name <- base_names[s]

  runs <- lapply(1:5, function(i) {
    file_path <- file.path(base_path, paste0(base_name, "_", i, ".log"))
    read.csv(file_path, skip = 12, header = FALSE, sep = " ")
  })

  time <- runs[[1]]$V1
  loss_matrix <- sapply(runs, function(x) log2(x$V2))

  mean_loss <- rowMeans(loss_matrix, na.rm = TRUE)

  all_time[[s]] <- time
  all_mean_loss[[s]] <- mean_loss

  global_min <- min(global_min, mean_loss, na.rm = TRUE)
  global_max <- max(global_max, mean_loss, na.rm = TRUE)
}

padding <- 0.05 * (global_max - global_min)
ylim_range <- c(global_min - padding, global_max + padding)

cell_type <- sub("_.*", "", base_names[1])
pdf_name <- paste0(cell_type, "_all_samples.pdf")

pdf(pdf_name, width = 7, height = 5)

plot(all_time[[1]],
     all_mean_loss[[1]],
     type = "l",
     col = colors[1],
     lwd = 2,
     ylim = ylim_range,
     xlab = "i",
     ylab = "Loss score (log2)",
     main = "")

for (s in 2:n_samples) {
  lines(all_time[[s]],
        all_mean_loss[[s]],
        col = colors[s],
        lwd = 2)
}

legend("topright",
       legend = base_names,
       col = colors[1:n_samples],
       lwd = 2,
       bty = "n")

dev.off()

