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

counts[["ge"]] <- ge

inter_sect <- get_intersect(counts, type = "row")

reordered_intersect <- lapply(
  seq_along(counts),
  function(x) {
    extract_and_order(
      counts[[x]],
      inter_sect,
      type = "row"
    )
  }
)

names(reordered_intersect) <- names(counts)

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

renamed_data <- rename_cell_types(reordered_intersect, remove_rename)

# Check if all rownames are equal
rn <- lapply(renamed_data, rownames)
all(sapply(rn, identical, rn[[1]]))

gc()

# Downsampling

downsample <- function(data, n) {
  cell_types <- lapply(data, colnames) %>% unlist
  cell_types_unique <- cell_types %>%
    unique
  downsample_index <- NULL
  for (i in seq_along(cell_types_unique)) {
    index <- which(cell_types == cell_types_unique[i])
    if (length(index) >= n) {
      index <- sample(index, n)
      downsample_index <- c(downsample_index, index)
    }
  }
  lengths <- lapply(data, ncol)
  cutoffs <- NULL
  init <- 0
  for (i in seq_along(lengths)) {
    init <- init + lengths[[i]]
    cutoffs <- c(cutoffs, init)
  }
  index_split <- NULL
  copy <- downsample_index
  for (i in seq_along(cutoffs)) {
    temp <- copy - cutoffs[i]
    id <- which(temp <= 0)
    if (i == 1) {
      index_split[[i]] <- copy[id]
      copy <- copy[-id]
    } else {
      new_id <- copy[id] - cutoffs[i - 1]
      index_split[[i]] <- new_id
      copy <- copy[-id]
    }
  }
  return(index_split)
}

save_matrix <- function(
  data,
  file,
  chunk_size,
  downsample_n,
  pseudobulk = FALSE
) {
  con <- file(file, "w")
  index <- downsample(data, downsample_n)
  col_names <- lapply(
    seq_along(data),
    function(x) {
      colnames(data[[x]])[index[[x]]]
    }
  )
  col_names <- unlist(col_names)
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
  d <- lapply(data, rownames)
  check <- all(sapply(d, identical, d[[1]]))
  if (check == FALSE) {
    stop("Row names are not equal")
  }
  n <- nrow(data[[1]])
  for (i in seq(1, n, by = chunk_size)) {
    end_row <-  min(i + chunk_size - 1, n)
    dense_chunk <- NULL
    for (j in seq_along(data)) {
      dense_chunk[[j]] <- as.matrix(
        data[[j]][i:end_row,
        index[[j]],
        drop = FALSE]
      )
    }
    dense_chunk <- Reduce(cbind, dense_chunk)
    # Convert one row at a time
    dense_chunk <- cbind.data.frame(GeneSymbol = rownames(data[[1]])[i:end_row],
                                    dense_chunk)
    write.table(dense_chunk,
                con,
                row.names = FALSE,
                col.names = FALSE,
                append = TRUE,
                quote = FALSE,
                sep = "\t")
  }
  close(con)
  if (pseudobulk == TRUE) {
    cells <- lapply(
      seq_along(data),
      function(x) {
        colnames(data[[x]])[-index[[x]]]
      }
    )
    cells <- unlist(cells)
    cells <- unique(cells)
    newdata <- NULL
    for (i in seq_along(cells)) {
      tempdata <- NULL
      for (j in seq_along(data)) {
        index <- which(colnames(data[[j]]) == cells[i])
        temp <- data[[j]][, index]
        cbind(tempdata, temp)
      }
      tempdata <- rowSums(tempdata)
      newdata <- cbind(newdata, tempdata)
    }
    colnames(newdata) <- cells
    pseudo <- newdata/rowSums(newdata)*1e6
    pseudo <- log2(pseudo + 1)
    return(pseudo)
  }
}

sample_sizes <- scan("sample_size.txt")

for (i in seq_along(sample_sizes)) {
  save_matrix(
    data = renamed_data,
    file = paste0("scRNA_combined_", sample_sizes[i], ".txt"),
    chunk_size = 500,
    downsample_n = sample_sizes[i]
  )
}

save_matrix(
  data = renamed_data,
  file = "scRNA_combined_2000.txt",
  chunk_size = 500,
  downsample_n = 2000,
  pseudobulk = TRUE
)
