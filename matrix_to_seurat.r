library(Seurat)
library(dplyr)
# BiocManager::install("EnsDb.Hsapiens.v79")
# 1. Convert from ensembl.gene to gene.symbol
library(EnsDb.Hsapiens.v79)
source("gene_expression_pipeline.R")

ge <- 2^gene_exp
ge <- rownames(ge)

rm(list = ls()[which(ls() != "ge")])
gc()

read_rds <- function(name) {
  if (!file.exists(paste0(name, ".rds"))) {
    path <- paste0("matrix_files/", name)
    raw_data <- Read10X(path)
    gc()
    file <- gzfile(paste0(path, "/metadata.csv.gz"))
    metadata <- read.csv(file)
    gc()
    data_so <- CreateSeuratObject(counts = raw_data, meta.data = metadata)
    saveRDS(data_so, paste0(name, ".rds"))
    return(data_so)
  } else {
    return(readRDS(paste0(name, ".rds")))
  }
}

clean_data <- function(rds) {
  cell_types <- rds$cell_type
  counts <- rds@assays$RNA$counts
  ensembl_genes <- rownames(counts)

  gene_ids <- ensembldb::select(EnsDb.Hsapiens.v79,
                                keys = ensembl_genes,
                                keytype = "GENEID",
                                columns = c("SYMBOL", "GENEID"))
  counts <- counts[rownames(counts) %in% gene_ids$GENEID, ]
  rownames(counts) <- gene_ids$SYMBOL
  colnames(counts) <- cell_types
  counts <- counts[rowSums(counts) > 0, ]
  counts <- counts[, colSums(counts) > 0]
  counts <- counts[rownames(counts) %in% ge, ]
  return(counts)
}

options(future.globals.maxSize = 7 * 1024^3)

name0 <- "5c9ab5a5-04a9-4282-9320-3b4d7b95131c"
name1 <- "4d3469a7-339f-40b3-92a3-22f7043545f8"
name2 <- "26f6625b-e76c-490a-beb1-aea16933cd6d"

count0 <- read_rds(name0)
gc()
counts0 <- clean_data(count0)

rm("count0")
gc()

count1 <- read_rds(name1)
counts1 <- clean_data(count1)
rm("count1")
gc()

count2 <- read_rds(name2)
counts2 <- clean_data(count2)
rm("count2")
gc()

counts1 <- counts1[rownames(counts1) %in% rownames(counts0), ]
counts0 <- counts0[rownames(counts0) %in% rownames(counts1), ]

all(rownames(counts1) %in% rownames(counts0))
all(rownames(counts0) %in% rownames(counts1))

dim(counts1)
dim(counts0)

any(duplicated(rownames(counts1)))
any(duplicated(rownames(counts0)))

duplicated <- rownames(counts0)[duplicated(rownames(counts0))]

# Remove duplicated gene
counts0 <- counts0[!rownames(counts0) %in% duplicated, ]
counts1 <- counts1[!rownames(counts1) %in% duplicated, ]

ge <- ge[rownames(ge) %in% rownames(counts0), ]

ge %>% rownames %>% head

rownames(counts0)[order(match(rownames(counts0), rownames(ge)))] %>% head

counts0 <- counts0[order(match(rownames(counts0), rownames(ge))), ]
counts1 <- counts1[order(match(rownames(counts1), rownames(ge))), ]

all(rownames(counts0) == rownames(counts1))

counts0[1:10, 1:10]
counts1[1:10, 1:10]

cells_to_remove_0 <- c("leukocyte")

counts0 <- counts0[, !colnames(counts0) %in% cells_to_remove]

cells_to_remove_1 <- c("cell of skeletal muscle",
  "chondrocyte",
  "Langerhans cell",
  "mast cell",
  "natural killer cell",
  "plasma cell",
  "plasmacytoid dendritic cell",
  "smooth muscle cell"
)

old_new_counts_0 <- rbind(
  c("skin fibroblast", "fibroblast")
)

old_new_counts_1 <- rbind(
  c("endothelial cell of lymphatic vessel", "endothelial cell"),
  c("endothelial cell of vascular tree", "endothelial cell"),
  c("suprabasal keratinocyte", "keratinocyte"),
  c("melanocyte of skin", "melanocyte"),
  c("CD8-positive, alpha-beta memory T cell, CD45RO-positive", "memory T cell")
)

rename_cell_types <- function(data, old_new) {
  for (i in seq_len(nrow(old_new))) {
    colnames(data) <- data %>%
      colnames() %>%
      gsub(old_new[i, 1], old_new[i, 2], x = .)
  }
  return(data)
}

counts0 <- rename_cell_types(counts0, old_new_counts_0)
counts1 <- rename_cell_types(counts1, old_new_counts_1)

sum(rownames(counts1) %in% rownames(counts0))
nrow(counts1)

sum(rownames(counts0) %in% rownames(counts1))
nrow(counts0)

counts0 <- counts0[rownames(counts0) %in% rownames(counts1),]
counts1 <- counts1[rownames(counts1) %in% rownames(counts0),]

counts0 <- counts0[order(match(rownames(counts0),rownames(counts1))),]

all(rownames(counts0) == rownames(counts1))

gc()

# Assuming `sparse_matrix` is a sparse matrix (e.g., of class dgCMatrix)
con <- file("scRNA_combined.txt", "w")

col_names <- c("GeneSymbol", colnames(counts1), colnames(counts0))
col_names <- col_names |> as.matrix() |> t()
write.table(col_names,
            con,
            row.names = FALSE,
            col.names = FALSE,
            append = TRUE,
            quote = FALSE,
            sep = "\t")

n <-  nrow(counts1)
chunk_size <- 500

gc()

d <- ncol(counts1) + ncol(counts0)

for (i in seq(1, n, by = chunk_size)) {
  if (any(rownames(counts1) != rownames(counts0))) {
    error(
    )
  }
  end_row <-  min(i + chunk_size - 1, n)
  dense_chunk1 <- as.matrix(counts1[i:end_row, , drop = FALSE])
  dense_chunk2 <- as.matrix(counts0[i:end_row, , drop = FALSE])
  # Convert one row at a time
  dense_chunk <- cbind(dense_chunk1, dense_chunk2)
  dense_chunk <- cbind.data.frame(GeneSymbol = rownames(counts1)[i:end_row],
                                  dense_chunk)
  write.table(dense_chunk,
              con,
              row.names = FALSE,
              col.names = FALSE,
              append = TRUE,
              quote = FALSE,
              sep = "\t")
  gc()
  if (i %% d / chunk_size / 10 == 0) {
    cat("*", sep = "")
  }
}

close(con)
