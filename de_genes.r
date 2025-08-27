library(doParallel)

immune <- groups[[length(groups)]]
skin <- groups[[length(groups) - 1]]

counts_immune <- counts[, counts$cell_type %in% immune]
counts_skin <- counts[, counts$cell_type %in% skin]

nrow(counts_immune) == nrow(counts_skin)

cl <- makeCluster(10)
registerDoParallel(cl)

t_genes <- foreach(i = seq_len(nrow(counts_immune)), .combine = rbind) %dopar% {
  d0 <- SingleCellExperiment::counts(counts_immune[i, ])
  d1 <- SingleCellExperiment::counts(counts_skin[i, ])
  t_test <- t.test(d0, d1)
  t <- t_test$statistic
  p <- t_test$p.value
  t_genes <- c(t.statistic = t, p.value = p)
  return(t_genes)
}

stopCluster(cl)

rownames(t_genes) <- rownames(counts_immune)

p_adj <- p.adjust(t_genes[, "p.value"], method = "fdr")

t_genes <- cbind(t_genes, p.adj = p_adj)
