# Install SCDC package

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("GEOquery"))
  install.packages("GEOquery")

if (!require("remotes")) {
  install.packages("remotes")
}
remotes::install_github("renozao/xbioc")

if (!require("SingleCellExperiment"))
  BiocManager::install("SingleCellExperiment")
library("SingleCellExperiment")

if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("meichendong/SCDC")

################################################
##### Load SCDC package
################################################

library(SCDC)

################################################
##### Load data
################################################

library(EnsDb.Hsapiens.v79)
source("gene_expression_pipeline.R")

ge <- 2^gene_exp

rm(list = ls()[which(ls() != "ge")])
gc()

remove_rename <- rjson::fromJSON(file = "remove_rename.json")

# remove_rename1 <- list()
# j <- 1
# for (i in seq_along(remove_rename)) {
#   if (remove_rename[[i]]$name == "26f6625b-e76c-490a-beb1-aea16933cd6d"){
#     next
#   } else {
#     remove_rename1[[j]] <- remove_rename[[i]]
#     j <- j + 1
#   }
# }
# remove_rename <- remove_rename1

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
  counts <- rds@assays$RNA$counts
  ensembl_genes <- rownames(counts)

  gene_ids <- ensembldb::select(EnsDb.Hsapiens.v79,
                                keys = ensembl_genes,
                                keytype = "GENEID",
                                columns = c("SYMBOL", "GENEID"))
  counts <- counts[rownames(counts) %in% gene_ids$GENEID, ]
  rownames(counts) <- gene_ids$SYMBOL
  counts <- counts[rowSums(counts) > 0, ]
  index <- colSums(counts) > 0
  counts <- counts[, index]
  counts <- counts[rownames(counts) %in% rownames(ge), ]
  # Remove duplicate rows
  counts <- counts[!duplicated(rownames(counts)), ]
  # counts <- SingleCellExperiment(assays = list(counts = counts))
  counts <- as.matrix(counts)
  pData <- cbind.data.frame(
    cell_type = rds$cell_type[index],
    donor_id = rds$donor_id[index]
  )
  counts <- ExpressionSet(
    assayData = counts,
    phenoData = AnnotatedDataFrame(pData)
  )
  return(counts)
}

options(future.globals.maxSize = 7 * 1024^3)

names_counts <- c(
  "palm_sole_hip",
  "inflammatory_skin_disease",
  "basal_cell_carcinoma"
)

lapply(
  seq_along(remove_rename),
  function(i) {
    name <- remove_rename[[i]]$name
    print(paste("Processing", name))
    count <- read_rds(name)
    print("Cleaning data")
    count <- clean_data(count)
    print("Running Quality Control for SCDC")
    count <- SCDC_qc(
      count,
      ct.varname = "cell_type",
      sample = "donor_id",
      scsetname = names_counts[i],
      ct.sub = unique(count$cell_type),
      qcthreshold = 0.7
    )
    print("Saving data")
    saveRDS(count, paste0(names_counts[i], ".rds"))
    gc()
  }
)
# names(counts) <- sapply(remove_rename, function(x) x$name)

gc()

########################################
### Import bulk data
########################################

bulk <- read.table(
  "gene_expression.txt",
  header = TRUE,
  sep = "\t"
)

rownames(bulk) <- bulk[, 1]
bulk <- bulk[, -1]
bulk <- as.matrix(bulk)

bulk <- Biobase::ExpressionSet(assayData = bulk)

#########################################
### Try SCDC
#########################################

counts
names_counts <- c(
  "palm_sole_hip",
  "inflammatory_skin_disease",
  "basal_cell_carcinoma"
)
lapply(
  seq_along(counts),
  function(i) {
    counts1_qc <- SCDC_qc(
      counts[[i]],
      ct.varname = "cell_type",
      sample = "donor_id",
      scsetname = names_counts[i],
      ct.sub = unique(counts[[i]]$cell_type),
      qcthreshold = 0.7
    )
    saveRDS(counts1_qc, paste0(names_counts[i], ".rds"))
  }
)


bulk_sc_counts1 <- SCDC_prop(
  bulk.eset = bulk,
  sc.eset = counts1_qc$sc.eset.qc,
  ct.varname = "cell_type",
  sample = "donor_id",
  ct.sub = unique(counts1_qc$sc.eset.qc$cell_type)
)

########################################
### Combine scRNA data
########################################

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
      dataframe$cell_type <- dataframe$cell_type %>%
        gsub(original, new, x = .)
    }
    renamed_data[[name]] <- dataframe
  }
  return(renamed_data)
}

renamed_data <- rename_cell_types(reordered_intersect, remove_rename)

column_bind <- function(list)  {
  for (i in seq_along(list)) {
    if (i == 1) {
      df <- SingleCellExperiment::counts(list[[i]])
      ct <- list[[i]]$cell_type
      di <- list[[i]]$donor_id
    } else {
      temp <- list[[i]]
      cell_type <- temp$cell_type
      donor_id <- temp$donor_id
      temp_df <- SingleCellExperiment::counts(temp)
      df <- MatrixExtra::cbind2(df, temp_df)
      ct <- c(ct, cell_type)
      di <- c(di, donor_id)
    }
  }
  df <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = df))
  df$cell_type <- ct
  df$donor_id <- di
  return(df)
}

counts <- column_bind(renamed_data)


