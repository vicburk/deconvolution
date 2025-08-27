# Download gene markers from PanglaoDB

file <- "data/PanglaoDB_markers_27_Mar_2020.tsv.gz"
f <- read.table(
  file,
  header = TRUE,
  sep = "\t",
  quote = ''
)

#####################
# Skin cells
#####################
# Schwann cell
# melanocyte
# keratinocyte
# endothelial cell
# fibroblast
# pericyte
#####################
# Immune cells
#####################
# plasma cell
# mast cell
# macrophage
# dendritic cell
# natural killer cell
# helper T cell
# regulatory T cell
# cytotoxic T cell
# innate lymphoid cell

# Unique cell types in dataset
f["cell.type"] %>% unique

# Immune cells of interest
immune <- c(
  "Plasma cells",
  "Mast cells",
  "Macrophages",
  "Dendritic cells",
  "Natural killer T cells",
  "T helper cells",
  "T regulatory cells",
  "T cytotoxic cells"
)
# immune <- c(
#   "Plasma cells",
#   "Mast cells",
#   "Macrophages",
#   "Dendritic cells",
#   "Natural killer T cells",
#   "T helper cells",
#   "T regulatory cells",
#   "T cytotoxic cells",
#   "Monocytes",
#   "Schwann cells",
#   "Chondrocytes",
#   "T memory cells",
#   "Basal cells"
# )

# Genes associated with immune cells
immune_genes <- f["official.gene.symbol"][unlist(f["cell.type"]) %in% immune,]
immune_genes <- immune_genes %>% unique

# Skin cells of interest
skin <- c(
  "Schwann cells",
  "Keratinoctyes",
  "Fibroblasts",
  "Endothelial cells",
  "Pericytes"
)

# skin <- c(
#   "Keratinoctyes",
#   "Fibroblasts",
#   "Endothelial cells",
#   "Pericytes"
# )

# Genes associated with skin cells
skin_genes <- f["official.gene.symbol"][unlist(f["cell.type"]) %in% skin,]
skin_genes <- skin_genes %>% unique

# Remove intersecting genes
immune_genes1 <- setdiff(immune_genes, skin_genes)
skin_genes1 <- setdiff(skin_genes, immune_genes)

########################
###### test
########################
###
###temp <- counts[, counts$cell_type %in% groups[["C5"]]]
###temp1 <- counts[, counts$cell_type %in% groups[["C4"]]]
###
###a <- t.test(counts(temp[1,]), counts(temp1[1,]))
###
###t_genes <- lapply(
###  seq_len(nrow(temp)),
###  function(x) {
###    t.test(counts(temp[x, ]), counts(temp1[x, ]))$p.value
###  }
###)
###
###library(parallel)
###t_genes <- mclapply(
###  seq_len(nrow(temp)),
###  function(x) {
###    t.test <- t.test(counts(temp[x, ]), counts(temp1[x, ]))
###    p <- t.test$p.value
###    t.stat <- t.test$statistic
###    if (i %% 1000 == 0) {
###      print(i)
###    }
###    return(c(p, t.stat))
###  },
###  mc.cores = 10
###)
###
###
###t_genes <- NULL
###for (i in seq_len(nrow(temp))) {
###  t <- t.test(counts(temp[i, ]), counts(temp1[i, ]))$p.value
###  t_genes <- c(t_genes, t)
###  if (i %% 100 == 0) {
###    print(i)
###  }
###}
###
###
###library(doParallel)
###ncores <- detectCores() - 1
###cl <- makeCluster(ncores, type = "SOCK")
###registerDoParallel(cl)
###
###t_genes <- foreach(i = seq_len(nrow(temp)), .combine = rbind) %dopar% {
###  t_test <- t.test(counts(temp[i, ]), counts(temp1[i, ]))
###  t <- t_test$statistic
###  p <- t_test$p.value
###  t_genes <- c(t.statistic = t, p.value = p)
###  if (i %% 100 == 0) {
###    print(i)
###  }
###  return(t_genes)
###}
###
###stopCluster(cl)
