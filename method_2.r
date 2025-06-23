library(Seurat)
library(dplyr)
# BiocManager::install("EnsDb.Hsapiens.v79")
# 1. Convert from ensembl.gene to gene.symbol
library(EnsDb.Hsapiens.v79)
source("gene_expression_pipeline.R")

dictionary <- read.xlsx("processed_data/dictionary.xlsx")

ge <- 2^gene_exp

rm(list = ls()[which(ls() != "ge")])
gc()

remove_rename <- rjson::fromJSON(file = "remove_rename.json")

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
  counts <- counts[rownames(counts) %in% rownames(ge), ]
  return(counts)
}

options(future.globals.maxSize = 7 * 1024^3)

counts <- list()
for (i in seq_along(remove_rename)) {
  name <- remove_rename[[i]]$name
  count <- read_rds(name)
  counts[[name]] <- clean_data(count)
}

rm("count")
gc()

#########################################################

rename_cell_types <- function(data, rename) {
  renamed_data <- list()
  for (i in seq_along(rename)) {
    name <- rename[[i]]$name
    dataframe <- data[[name]]
    for (j in seq_along(rename[[i]]$rename)) {
      original <- names(rename[[i]]$rename)[j]
      new <- rename[[i]]$rename[[j]]
      colnames(dataframe) <- dataframe %>%
        colnames() %>%
        gsub(original, new, x = .)
    }
    renamed_data[[name]] <- dataframe
  }
  return(renamed_data)
}

renamed_data <- rename_cell_types(counts, remove_rename)

#############################
### Downsample the data ######
### #########################

tolerance <- 4e8

downsample <- function(data, tolerance = 4e8) {
  downsample <- list()
  for (i in seq_along(data)) {
    if (prod(dim(data[[i]])) > tolerance) {
      size <- prod(dim(data[[i]]))
      c <- tolerance / size
      m <- round(ncol(data[[i]]) * c)
      downsample[[i]] <- data[[i]][, sample(ncol(data[[i]]), m)]
    } else {
      downsample[[i]] <- data[[i]]
    }
  }
  names(downsample) <- names(data)
  return(downsample)
}

downsampled_data <- downsample(renamed_data)

#############################
### Save the data in chunks ###
#############################

save_matrix <- function(
  data,
  chunk_size
) {
  for (i in seq_along(data)) {
    file <- paste0("scRNA_", names(data)[i], ".txt")
    con <- file(file, "w")
    col_names <- colnames(data[[i]])
    col_names <- c("GeneSymbol", col_names)
    col_names <- col_names %>% as.matrix() %>% t()
    write.table(
      col_names,
      con,
      row.names = FALSE,
      col.names = FALSE,
      append = TRUE,
      quote = FALSE,
      sep = "\t"
    )
    for (j in seq(1, nrow(data[[i]]), by = chunk_size)) {
      end_row <- min(j + chunk_size - 1, nrow(data[[i]]))
      dense_chunk <- data[[i]][j:end_row, ]
      dense_chunk <- cbind.data.frame(
        GeneSymbol = rownames(dense_chunk),
        dense_chunk
      )
      write.table(
        dense_chunk,
        con,
        row.names = FALSE,
        col.names = FALSE,
        append = TRUE,
        quote = FALSE,
        sep = "\t"
      )
    }
    close(con)
  }
}

chunk_size <- 500
save_matrix(downsampled_data, chunk_size)
