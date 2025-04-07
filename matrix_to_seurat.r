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

# Function to find the intersection
get_intersect <- function(sets, type = c("row", "col", "vector")) {
  # This will throw an error if the type is not valid
  type <- match.arg(type)
  if (type == "row") {
    sets <- lapply(sets, rownames)
  } else if (type == "col") {
    sets <- lapply(sets, colnames)
  }
  Reduce(intersect, sets)
}

# Function to extract and order elements based on the intersection
extract_and_order <- function(original, intersection, type = c("row", "col", "vector")) {
  type <- match.arg(type)
  if (type == "row") {
    rows <- rownames(original)
    return(original[match(intersection, rows), ])
  } else if (type == "col") {
    cols <- colnames(original)
    return(original[, match(intersection, cols)])
  } else {
    return(original[match(intersection, original)])
  }
}

sets <- list(
  ge = ge,
  counts0 = counts0,
  counts1 = counts1,
  counts2 = counts2
)

inter_sect <- get_intersect(setsm, type = "row")

reordered_intersect <- lapply(
  sets,
  function(x) {
    extract_and_order(
      x,
      inter_sect,
      type = "row"
    )
  }
)

counts0 <- reordered_intersect$counts0
counts1 <- reordered_intersect$counts1
counts2 <- reordered_intersect$counts2

cells_to_remove0 <- c(
  "cell of skeletal muscle",
  "chondrocyte",
  "mast cell",
  "plasma cell",
  "smooth muscle cell",
  "Schwann cell"
)

cells_to_remove1 <- c("leukocyte")

cells_to_recove2 <- c(
  "Schwann cell",
  "plasma cell",
  "mast cell"
)

old_new_counts0 <- rbind(
  c("endothelial cell of lymphatic vessel", "endothelial cell"),
  c("endothelial cell of vascular tree", "endothelial cell"),
  c("suprabasal keratinocyte", "keratinocyte"),
  c("melanocyte of skin", "melanocyte"),
  c("CD8-positive, alpha-beta memory T cell, CD45RO-positive", "memory T cell"),
  c("basal cell of epidermis", "basal cell"),
  c("plasmacytoid dendritic cell", "dendritic cell"),
  c("conventional dendritic cell", "dendritic cell")
)

old_new_counts1 <- rbind(
  c("skin fibroblast", "fibroblast")
)

old_new_counts2 <- rbind(
  c("conventional dendritic cell", "dendritic cell"),
  c("monocyte-derived dendritic cell", "dendritic cell"),
  c("inflammatory macrophage", "macrophage"),
  c("endothelial cell of vascular tree", "endothelial cell"),
  c("endothelial cell of lymphatic vessel", "endothelial cell"),
  c("skin fibroblast", "fibroblast")
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
counts2 <- rename_cell_types(counts2, old_new_counts_2)

gc()

# Assuming `sparse_matrix` is a sparse matrix (e.g., of class dgCMatrix)
con <- file("scRNA_combined.txt", "w")

col_names <- c(
  "GeneSymbol",
  colnames(counts0),
  colnames(counts1),
  colnames(counts2)
)

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

check <- any(rownames(counts1) != rownames(counts0) | any(rownames(counts1 != rownames(counts2))))

for (i in seq(1, n, by = chunk_size)) {
  if (check) {
    error(
    )
  }
  end_row <-  min(i + chunk_size - 1, n)
  dense_chunk0 <- as.matrix(counts0[i:end_row, , drop = FALSE])
  dense_chunk1 <- as.matrix(counts1[i:end_row, , drop = FALSE])
  dense_chunk2 <- as.matrix(counts2[i:end_row, , drop = FALSE])
  # Convert one row at a time
  dense_chunk <- cbind(dense_chunk0, dense_chunk1, dense_chunk2)
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
