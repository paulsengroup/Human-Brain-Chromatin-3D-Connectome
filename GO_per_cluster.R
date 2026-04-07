library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(openxlsx)

input_dir  <- ""
output_dir <- ""

files <- list.files(
  input_dir,
  pattern = ".tsv$",
  full.names = TRUE
)

ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror = "useast"
)

get_genes_region <- function(chr, start, end, cluster_id) {
  
  res <- getBM(
    attributes = c(
      "ensembl_gene_id",
      "external_gene_name",
      "entrezgene_id",
      "gene_biotype",
      "chromosome_name",
      "start_position",
      "end_position"
    ),
    filters = c("chromosome_name", "start", "end"),
    values = list(chr, start, end),
    mart = ensembl
  )
  
  res$cluster <- cluster_id
  return(res)
}

for (f in files) {
  
  base_name <- basename(f)
  base_name <- sub("\\.tsv$", "", base_name)
  print(base_name)
  cluster <- read.csv(f, sep = "\t")
  
  cluster <- cluster %>%
    separate(bead_id, into = c("chromosome", "positions"), sep = ":", remove = FALSE) %>%
    separate(positions, into = c("start", "end"), sep = "-")
  
  results <- bind_rows(
    mapply(
      get_genes_region,
      cluster$chromosome,
      cluster$start,
      cluster$end,
      cluster$cluster,
      SIMPLIFY = FALSE
    )
  )
  
  gene_list <- results %>%
    filter(!is.na(external_gene_name)) %>%
    group_by(cluster) %>%
    summarise(genes = list(unique(external_gene_name))) %>%
    deframe()
  
  ego <- compareCluster(
    geneCluster = gene_list,
    fun = "enrichGO",
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  pdf(file.path(output_dir, paste0("GO_", base_name, ".pdf")),
      width = 8, height = 12)
  print(dotplot(ego, showCategory = 3))
  dev.off()
  
  write.xlsx(
    ego@compareClusterResult,
    file.path(output_dir, paste0("GO_", base_name, ".xlsx"))
  )
}
