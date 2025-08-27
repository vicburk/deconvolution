library(MASS)
library(Seurat)
library(dplyr)
library(nnls)
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


resample_to_proportions <- function(data, target_props, column = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Handle different input types
  if (is.data.frame(data)) {
    if (is.null(column)) {
      stop("Column must be specified when input is a data frame")
    }
    if (!column %in% names(data)) {
      stop("Specified column not found in data frame")
    }
    # Store original data frame and extract column for sampling
    original_df <- data
    sampling_vector <- data[[column]]
  } else {
    sampling_vector <- data
  }
  
  # Input validation
  if (sum(target_props) > 1) {
    stop("Sum of target proportions cannot exceed 1")
  }
  
  if (is.null(names(target_props))) {
    stop("target_props must be a named vector")
  }
  
  # Convert data to table
  freq_table <- table(sampling_vector)
  
  # Check if all names in target_props exist in the data
  if (!all(names(target_props) %in% names(freq_table))) {
    stop("Some categories in target_props don't exist in the data")
  }
  
  # Handle remaining proportion
  remaining_prop <- 1 - sum(target_props)
  if (remaining_prop > 0) {
    unnamed_cats <- setdiff(names(freq_table), names(target_props))
    if (length(unnamed_cats) == 0) {
      stop("No categories available for remaining proportion")
    }
    remaining_props <- rep(remaining_prop / length(unnamed_cats), length(unnamed_cats))
    names(remaining_props) <- unnamed_cats
    target_props <- c(target_props, remaining_props)
  }
  
  # Reorder freq_table to match target_props order
  freq_table <- freq_table[names(target_props)]
  
  # Calculate current proportions
  current_props <- prop.table(freq_table)
  
  # Calculate scaling factor based on limiting category
  scaling_factor <- min(freq_table/sum(freq_table) / target_props)
  
  # Calculate total possible sample size
  total_size <- sum(freq_table) * scaling_factor
  
  # Calculate target sizes for each category
  target_sizes <- round(total_size * target_props)
  names(target_sizes) <- names(target_props)
  
  # Check if we have enough data
  if (any(target_sizes > freq_table)) {
    warning("Not enough data in some categories to achieve desired proportions")
    target_sizes <- pmin(target_sizes, freq_table)
  }
  
  # Sample from each category
  indices <- lapply(names(target_sizes), function(category) {
    category_indices <- which(sampling_vector == category)
    sample(category_indices, size = target_sizes[category], replace = FALSE)
  })
  
  # Combine all indices
  selected_indices <- sort(unlist(indices))
  
  # Return results based on input type
  if (is.data.frame(data)) {
    result <- original_df[selected_indices, ]
    # Add attributes to the data frame
    attr(result, "original_props") <- current_props
    attr(result, "target_props") <- target_props
    attr(result, "achieved_props") <- prop.table(table(result[[column]]))
  } else {
    result <- sampling_vector[selected_indices]
    attr(result, "original_props") <- current_props
    attr(result, "target_props") <- target_props
    attr(result, "achieved_props") <- prop.table(table(result))
  }
  
  return(result)
}


save_matrix <- function(
  data,
  downsample_n,
  pseudobulk = FALSE,
  pseudo_n = 3
) {
  index <- downsample(data, downsample_n)
  col_names <- lapply(
    seq_along(data),
    function(x) {
      colnames(data[[x]])[index[[x]]]
    }
  )
  col_names <- unlist(col_names)
  reference_matrix <- NULL
  for (i in seq_along(unique(col_names))) {
    cell_type <- unique(col_names)[i]
    data_by_cell_type <- 0
    for (j in seq_along(data)) {
      df <- data[[j]][, colnames(data[[j]]) == cell_type]
      df <- rowSums(df)
      data_by_cell_type <- data_by_cell_type + df
    }
    data_by_cell_type <- data_by_cell_type/sum(data_by_cell_type)*1e6
    data_by_cell_type <- log2(data_by_cell_type + 1)
    reference_matrix <- cbind(reference_matrix, data_by_cell_type)
  }
  colnames(reference_matrix) <- unique(col_names)
  if (pseudobulk == TRUE) {
    inverse <- lapply(
      seq_along(data),
      function(x) {
        data[[x]][, -index[[x]]]
      }
    )
    lengths <- lapply(inverse, ncol)
    ncells <- lengths %>% unlist %>% sum
    cells <- lapply(
      inverse,
      colnames
    ) %>% unlist
    cells <- cbind.data.frame(
      cells = cells,
      index = seq_len(ncells)
    )
    scramble <- sample(
      seq_len(ncells),
    )
    group <- scramble %% pseudo_n
    group <- group + 1
    cells_by_group <- split(cells, group)
    target_prop <- cbind(
      "keratinocyte" = c(0.85, 0.6, 0.8),
      "fibroblast" = c(0.05, 0.2, 0.1),
      "endothelial cell" = c(0.05, 0.1, 0.01),
      "smooth muscle cell" = c(0, 0, 0),
      "cell of skeletal muscle" = c(0, 0, 0)
    )
    cell_group_index <- lapply(
      seq_along(cells_by_group),
      function(x) {
        resample_to_proportions(
          cells_by_group[[x]],
          target_props = target_prop[x, ],
          column = "cells"
        )
      }
    )
 
    cutoffs <- NULL
    init <- 0
    for (i in seq_along(lengths)) {
      init <- init + lengths[[i]]
      cutoffs <- c(cutoffs, init)
    }

    index_split <- list()
    for (i in seq_along(cell_group_index)) {
      index_split[[i]] <- list()
      index <- cell_group_index[[i]]$index
      copy <- index
      for (j in seq_along(cutoffs)) {
        temp <- copy - cutoffs[j]
        id <- which(temp <= 0)
        if (j == 1) {
          index_split[[i]][[j]] <- copy[id]
          copy <- copy[-id]
        } else {
          new_id <- copy[id] - cutoffs[j - 1]
          index_split[[i]][[j]] <- new_id
          copy <- copy[-id]
        }
      }
    }

    pseudo <- NULL
    for (i in seq_along(index_split)) {
      bulk <- NULL
      for (j in seq_along(index_split[[i]])) {
        samp <- inverse[[j]][, index_split[[i]][[j]]]
        samp <- samp %>% rowSums
        bulk <- cbind(bulk, samp)
      }
      samp <- bulk %>% rowSums
      samp <- samp/sum(samp)*1e6
      samp <- log2(samp + 1)
      pseudo[[i]] <- samp
    }
    pseudo <- Reduce(
      cbind,
      pseudo
    )
    true_prop <- lapply(
      seq_along(cell_group_index),
      function(x) {
        attr(cell_group_index[[x]], "achieved_props")
      }
    )
    true_prop <- Reduce(
      cbind,
      true_prop
    )
    return(
      list(
        reference_matrix = reference_matrix,
        pseudo = pseudo,
        true_prop = true_prop
      )
    )
  }
  return(reference_matrix)
}

#

processed_data <- save_matrix(
  data = renamed_data,
  downsample_n = 2000,
  pseudobulk = TRUE
)

data <- processed_data$reference_matrix
pseudo <- processed_data$pseudo
true_prop <- processed_data$true_prop
true_prop <- true_prop[colnames(data), ]

########################################
# Load bulk data
########################################

bulk <- read.table(
  "gene_expression.txt",
  header = TRUE,
  sep = "\t"
)
rownames(bulk) <- bulk[, 1]
bulk <- bulk[, -1]
bulk <- as.matrix(bulk)

bulk <- bulk[match(rownames(data), rownames(bulk)), ]

#######################################
# Load signature matrix from CIBERSORT
#######################################

file <- "pseudo_output_2000/CIBERSORTx_cell_type_sourceGEP.txt"

cibersort <- read.table(
  file,
  header = TRUE,
  sep = "\t"
)
rownames(cibersort) <- cibersort[, 1]
cibersort <- cibersort[, -1]
cibersort <- as.matrix(cibersort)

reorder <- match(rownames(data), rownames(cibersort))

cibersort <- cibersort[reorder, ]

cibersort <- na.omit(cibersort)

pseudo_0 <- pseudo[rownames(cibersort), ]
bulk_0 <- bulk[rownames(cibersort), ]

##################################
# NNLS Deconvolution
##################################

# Manually calculated signature matrix

# Pseudo Bulk Deconvolution

custom_deconvolution <- function(signature, bulk, method, cell_types = NULL, seed = 1) {
  if (!is.null(cell_types)) {
    colnames(signature) <- gsub("\\.", " ", colnames(signature))
    signature <- signature[, cell_types]
  }
  if (!any(c("NNLS","lasso") %in% method)) {
    stop("method must be NNLS or lasso")
  }

  if (method == "NNLS") {
    props <- lapply(
      seq_len(ncol(bulk)),
      function(i) {
        model <- nnls(signature, bulk[, i])
        props <- coef(model)
        props <- props/sum(props)
        return(props)
      }
    )
    props <- Reduce(
      cbind,
      props
    )
    rownames(props) <- colnames(signature)
  }

  if (method == "lasso") {
    set.seed(seed)
    props <- lapply(
      seq_len(ncol(bulk)),
      function(i) {
        model <- cv.glmnet(signature, bulk[, i], alpha = 1, lower.limit = 0, intercept = FALSE)
        props <- coef(model)
        props <- props/sum(props)
        return(props)
      }
    )
    props <- Reduce(
      cbind,
      props
    )
  }
  return(props)
}

# props_pseudo_nnls
ppn <- custom_deconvolution(signature = data, bulk = pseudo, method = "NNLS")

# Bulk Deconvolution

# props_bulk_nnls

pbn <- custom_deconvolution(signature = data, bulk = bulk, method = "NNLS")

# CIBERSORT signature matrix

# Pseudo Bulk Deconvolution

# props_pseudo_nnls_cibersort

ppnc <- custom_deconvolution(signature = cibersort, bulk = pseudo_0, method = "NNLS")

# Bulk Deconvolution

# props_bulk_nnls_cibersort

pbnc <- custom_deconvolution(signature = cibersort, bulk = bulk_0, method = "NNLS")

#################################################
### Lasso
#################################################

library(glmnet)

# Manually calculated signature matrix

# Pseudo Bulk Deconvolution

# props_pseudo_lasso

ppl <- custom_deconvolution(signature = data, bulk = pseudo, method = "lasso")

# Bulk Deconvolution

# props_bulk_lasso

pbl <- custom_deconvolution(signature = data, bulk = bulk, method = "lasso")

# CIBERSORT signature matrix

# Pseudo Bulk Deconvolution

# props_pseudo_lasso_cibersort

pplc <- custom_deconvolution(signature = cibersort, bulk = pseudo_0, method = "lasso")

# Bulk Deconvolution

# props_bulk_lasso_cibersort

pblc <- custom_deconvolution(signature = cibersort, bulk = bulk_0, method = "lasso")

########################################
### No Fibroblast
########################################

cells <- colnames(data)

no_fib <- cells[cells != "fibroblast"]

##################################
# NNLS Deconvolution
##################################

# Manually calculated signature matrix

# Pseudo Bulk Deconvolution

# props_pseudo_nnls_no_fib

ppn_no_fib <- custom_deconvolution(signature = data, bulk = pseudo, method = "NNLS", cell_types = no_fib)

# Bulk Deconvolution

# props_bulk_nnls_no_fib

pbn_no_fib <- custom_deconvolution(signature = data, bulk = bulk, method = "NNLS", cell_types = no_fib)

# CIBERSORT signature matrix

# Pseudo Bulk Deconvolution

# props_pseudo_nnls_cibersort_no_fib

ppnc_no_fib <- custom_deconvolution(signature = cibersort, bulk = pseudo_0, method = "NNLS", cell_types = no_fib)

# Bulk Deconvolution

# props_bulk_nnls_cibersort_no_fib

pbnc_no_fib <- custom_deconvolution(signature = cibersort, bulk = bulk_0, method = "NNLS", cell_types = no_fib)

#################################################
### Lasso
#################################################

# Manually calculated signature matrix

# Pseudo Bulk Deconvolution

# props_pseudo_lasso_no_fib

ppl_no_fib <- custom_deconvolution(signature = data, bulk = pseudo, method = "lasso", cell_types = no_fib)

# Bulk Deconvolution

# props_bulk_lasso_no_fib

pbl_no_fib <- custom_deconvolution(signature = data, bulk = bulk, method = "lasso", cell_types = no_fib)

# CIBERSORT signature matrix

# Pseudo Bulk Deconvolution

# props_pseudo_lasso_cibersort_no_fib

pplc_no_fib <- custom_deconvolution(signature = cibersort, bulk = pseudo_0, method = "lasso", cell_types = no_fib)

# Bulk Deconvolution

# props_bulk_lasso_cibersort_no_fib

pblc_no_fib <- custom_deconvolution(signature = cibersort, bulk = bulk_0, method = "lasso", cell_types = no_fib)

########################################
### Clustering
########################################

distance <- 1 - cor(data, method = "spearman")
distance <- as.dist(distance)

hc <- hclust(distance, method = "complete")

plot(hc)

################

dist_cibersort <- 1 - cor(cibersort, method = "spearman")
dist_cibersort <- as.dist(dist_cibersort)

hc_cibersort <- hclust(dist_cibersort, method = "complete")

plot(hc_cibersort)

##################################
### No Fibroblast, Pericyte, melanocyte, cytotoxic T cell, or Endothelial
##################################

exclude <- c("fibroblast", "pericyte", "endothelial cell")
no_fpe <- cells[!(cells %in% exclude)]

ppn_no_fpe <- custom_deconvolution(signature = data, bulk = pseudo, method = "NNLS", cell_types = no_fpe)

pbn_no_fpe <- custom_deconvolution(signature = data, bulk = bulk, method = "NNLS", cell_types = no_fpe)
