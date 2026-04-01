library(igraph)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tools)

dir_path <- "sim_input/"

files <- list.files(path = dir_path, pattern = "hicexplorer\\.bed$", full.names = TRUE)

all_centrality <- list()
for (file in files) {
  sample_name <- tools::file_path_sans_ext(basename(file))
  sample_name <- gsub("_fdr0\\.01_cis2_trans2_hicexplorer", "", sample_name)
  sample_name <- gsub("merge", "", sample_name)
  
  data <- read.csv(file, sep = "\t", header = FALSE)
  
  data$ID1 <- paste0("chr", data$V1, "_", data$V2, "_", data$V3)
  data$ID2 <- paste0("chr", data$V4, "_", data$V5, "_", data$V6)
  
  data <- data[, c("ID1", "ID2")]
  
  g <- graph_from_edgelist(as.matrix(data), directed = FALSE)
  
  
  V(g)$degree <- degree(g)
  V(g)$eig <- eigen_centrality(g)$vector
  V(g)$closeness <- closeness(g)
  V(g)$betweenness <- betweenness(g)
  
  clique_size_file <- read.csv(sub("\\.bed$", "_max_clique_size.tsv", file), sep = "\t")
  clique_size_file$node = paste0("chr", clique_size_file$node)
  
  V(g)$max_clique <- clique_size_file$max_clique_size[
    match(V(g)$name, clique_size_file$node)
  ]
  
  comm <- cluster_louvain(g)
  
  global_clust <- transitivity(g, type = "global")
  avg_local_clust <- transitivity(g, type = "average")
  assort_deg <- assortativity_degree(g)
  
  diam <- diameter(g, directed = FALSE, weights = NA)
  avg_path <- average.path.length(g, directed = FALSE, unconnected = TRUE)
  
  
  centrality <- data.frame(
    node       = V(g)$name,
    degree     = log2(V(g)$degree),
    closeness  = log2(V(g)$closeness),
    betweenness = log2(V(g)$betweenness),
    max_clique = V(g)$max_clique,
    eigenvector = log2(V(g)$eig),
    community   = length(comm),
    
    diameter = diam,
    avg_path_length = avg_path,
    
    global_clust   = global_clust,
    avg_local_clust   = avg_local_clust,
    assortativity_degree   = assort_deg,
    
    
    
    
    sample     = sample_name
  )
  
  all_centrality[[sample_name]] <- centrality
}

centrality_all <- bind_rows(all_centrality)

centrality_long <- centrality_all %>%
  tidyr::pivot_longer(cols = c(degree, closeness, betweenness, eigenvector, 
                               community, global_clust, avg_local_clust, 
                               diameter, avg_path_length, max_clique),
                      names_to = "metric", values_to = "value")
centrality_long <- centrality_all %>%
  tidyr::pivot_longer(cols = c(community, global_clust, avg_local_clust),
                      names_to = "metric", values_to = "value")

centrality_long$sample <- factor(centrality_long$sample, levels = 
                                   c("GABA_infancy", "GABA_childhood", "GABA_adolescence", "GABA_adulthood", "GABA_lateadulthood",
                                     "GLU_infancy", "GLU_childhood", "GLU_adolescence", "GLU_adulthood", "GLU_lateadulthood",
                                     "MgAs_infancy", "MgAs_childhood", "MgAs_adolescence", "MgAs_adulthood", "MgAs_lateadulthood",
                                     "Olig_infancy", "Olig_childhood", "Olig_adolescence", "Olig_adulthood", "Olig_lateadulthood"))


centrality_long = centrality_long[order(centrality_long$sample), ]
centrality_long$cell_type <- sub("_.*", "", centrality_long$sample)
centrality_long$sample <- gsub("_", " ", centrality_long$sample)

box_stats <- centrality_long %>%
  group_by(sample, metric) %>%
  summarise(
    middle = median(value, na.rm = TRUE),
    std = sd(value[is.finite(value)], na.rm = TRUE),
    .groups = "drop"
  )

box_stats$sample = sub("lateadulthood", "late adulthood", box_stats$sample)
box_stats$sample <- factor(box_stats$sample, levels = 
                             c("GABA infancy", "GABA childhood", "GABA adolescence", "GABA adulthood", "GABA late adulthood",
                               "GLU infancy", "GLU childhood", "GLU adolescence", "GLU adulthood", "GLU late adulthood",
                               "MgAs infancy", "MgAs childhood", "MgAs adolescence", "MgAs adulthood", "MgAs late adulthood",
                               "Olig infancy", "Olig childhood", "Olig adolescence", "Olig adulthood", "Olig late adulthood"))

box_stats$sample <- factor(box_stats$sample, levels = 
                             c("GABA infancy", "GLU infancy", "MgAs infancy", "Olig infancy",
                               "GABA childhood", "GLU childhood", "MgAs childhood", "Olig childhood",
                               "GABA adolescence", "GLU adolescence", "MgAs adolescence", "Olig adolescence",
                               "GABA adulthood", "GLU adulthood", "MgAs adulthood", "Olig adulthood",
                               "GABA late adulthood", "GLU late adulthood", "MgAs late adulthood", "Olig late adulthood"))

box_stats = box_stats[order(box_stats$sample), ]
box_stats <- box_stats %>%
  separate(sample, into = c("cell_type", "age_group"), sep = " ", extra = "merge")

box_stats$age_group <- factor(box_stats$age_group,
                              levels = c("infancy", "childhood", "adolescence", "adulthood", "late adulthood"))


cis = c(5478,2168,2703,1094,985,1910,2891,5876,2259,2306,
        339,546,702,415,355,528,2866,1231,731,7065)
trans = c(10693,11467,21323,15582,7451,10885,12322,12303,8953,
          12314,4228,7404,8745,7865,7227,4894,4540,7534,4365,1271)


nchg_cis = data.frame(cell_type = c(rep("GABA", 5), rep("GLU", 5), rep("MgAs", 5), rep("Olig", 5)),
                      age_group= rep(c("infancy","childhood", "adolescence","adulthood", "late adulthood"), 4),
                      metric = rep("cis", 20),
                      middle = cis/(cis + trans),
                      std = rep(0, 20))
nchg_trans = data.frame(cell_type = c(rep("GABA", 5), rep("GLU", 5), rep("MgAs", 5), rep("Olig", 5)),
                        age_group= rep(c("infancy","childhood", "adolescence","adulthood", "late adulthood"), 4),
                        metric = rep("trans", 20),
                        middle = trans/(cis + trans),
                        std = rep(0, 20))
box_stats = rbind(box_stats, nchg_cis, nchg_trans)

metric_order <- c(
  "cis", "trans", "degree", "closeness", "betweenness", "eigenvector",
  "avg_local_clust", "global_clust", "community", "avg_path_length", "diameter", "max_clique"
)
metric_order <- c(
  "cis", "trans",
  "avg_local_clust", "global_clust", "community"
)

box_stats$metric <- factor(box_stats$metric, levels = metric_order, labels = c("Cis share", "Trans share", "Avg local clustering", "Global clustering", "Community"))


write.table(centrality_all, "graph_properties_complete_table.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
