#################################################
##### Install MuSiC and dependencies
#################################################

packages <- c(
  "devtools",
  "TOAST",
  "reshape2",
  "cowplot",
  "RColorBrewer",
  "openxlsx",
  "readxl",
  "rjson",
  "Seurat",
  "MatrixExtra"
)

need <- ! unlist(lapply(packages, require, character.only = TRUE))

if (any(need)) {
  install.packages(packages[need])
}

if (!require("MuSiC")) {
  devtools::install_github("xuranw/MuSiC")
}

lapply(
  c(
    packages,
    "MuSiC",
    "SingleCellExperiment" # Without this package loaded, `music_prop()` will not work
                           # because it uses the `music_basis()` function, which uses
                           # the the `counts()` function from SingleCellExperiment.
                           # `counts()` is also a function in the BiocGenerics package.
  ),
  require,
  character.only = TRUE
)

source("MuSiC_functions.r")

# BiocManager::install("EnsDb.Hsapiens.v79")
# 1. Convert from ensembl.gene to gene.symbol
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
  cell_types <- rds$cell_type
  counts <- rds@assays$RNA$counts
  ensembl_genes <- rownames(counts)

  gene_ids <- ensembldb::select(EnsDb.Hsapiens.v79,
                                keys = ensembl_genes,
                                keytype = "GENEID",
                                columns = c("SYMBOL", "GENEID"))
  counts <- counts[rownames(counts) %in% gene_ids$GENEID, ]
  rownames(counts) <- gene_ids$SYMBOL
  counts <- counts[rowSums(counts) > 0, ]
  counts <- counts[, colSums(counts) > 0]
  counts <- counts[rownames(counts) %in% rownames(ge), ]
  counts <- SingleCellExperiment(assays = list(counts = counts))
  counts$cell_type <- cell_types
  counts$donor_id <- rds$donor_id
  return(counts)
}

options(future.globals.maxSize = 7 * 1024^3)

counts <- list()
counts <- lapply(
  seq_along(remove_rename),
  function(i) {
    name <- remove_rename[[i]]$name
    count <- read_rds(name)
    clean_data(count)
  }
)
names(counts) <- sapply(remove_rename, function(x) x$name)

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

# Clustering
count_basis <- music_basis(
  counts,
  clusters = "cell_type",
  samples = "donor_id",
)

# Plot the dendrogram of design matrix and cross-subject mean of realtive abundance
par(mfrow = c(1, 2))
d <- dist(t(log(count_basis$Disgn.mtx + 1e-6)), method = "euclidean")
# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete" )
# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1, main = 'Cluster log(Design Matrix)')
d <- dist(t(log(count_basis$M.theta + 1e-8)), method = "euclidean")
# Hierarchical clustering using Complete Linkage
# hc2 <- hclust(d, method = "complete" )
hc2 <- hclust(d, method = "complete")

png("MuSiC_cluster.png", width = 3000, height = 3000, res = 300)
# Plot the obtained dendrogram
plot(hc2, cex = 0.6, hang = -1, main = 'Cluster log(Mean of RA)')
dev.off()

######################
### Create Groups
######################

# groups <- list(
#   C1 = "cell of skeletal muscle",
#   C2 = "smooth muscle cell",
#   C3 = "leukocyte",
#   C4 = c(
#     "melanocyte",
#     "keratinocyte",
#     "endothelial cell",
#     "fibroblast",
#     "pericyte"
#   ),
#   C5 = c(
#     "monocyte",
#     "Schwann cell",
#     "chondrocyte",
#     "memory T cell",
#     "basal cell",
#     "plasma cell",
#     "mast cell",
#     "macrophage",
#     "dendritic cell",
#     "natural killer cell",
#     "helper T cell",
#     "regulatory T cell",
#     "cytotoxic T cell",
#     "innate lymphoid cell"
#   )
# )

groups <- list(
  C1 = "cell of skeletal muscle",
  C2 = "smooth muscle cell",
  C3 = "monocyte",
  C4 = "chondrocyte",
  C5 = c(
    "memory T cell",
    "basal cell"
  ),
  C6 = "leukocyte",
  C7 = "Langerhans cell",
  C8 = c(
    "Schwann cell",
    "melanocyte",
    "keratinocyte",
    "endothelial cell",
    "fibroblast",
    "pericyte"
  ),
  C9 = c(
    "plasma cell",
    "mast cell",
    "macrophage",
    "dendritic cell",
    "natural killer cell",
    "helper T cell",
    "regulatory T cell",
    "cytotoxic T cell",
    "innate lymphoid cell"
  )
)

cl_type <- counts$cell_type

for (i in seq_along(groups)) {
  cl_type[cl_type %in% groups[[i]]] <- names(groups)[i]
}

counts$clusterType <- as.factor(cl_type)

# source("markers.r")

null_markers <- paste0("C", seq(1,length(groups) - 1))
markers <- NULL
markers <- lapply(
  null_markers,
  function(x) {
    markers[[x]] <- NULL
    return(markers)
  }
)
names(markers) <- null_markers

# markers[[paste0("C", length(groups) - 1)]] <- skin_genes1
# markers[[paste0("C", length(groups))]] <- immune_genes1

markers[[paste0("C", length(groups))]] <- rownames(counts)

bulk <- read.table(
  "gene_expression.txt",
  header = TRUE,
  sep = "\t"
)

rownames(bulk) <- bulk[, 1]
bulk <- bulk[, -1]
bulk <- as.matrix(bulk)
ge <- log2(ge)

Est.skin.bulk = music_prop.cluster(
  bulk.mtx = bulk,
  sc.sce = counts,
  group.markers = markers,
  clusters = 'cell_type',
  groups = 'clusterType',
  samples = 'donor_id',
  clusters.type = groups
)

deconvolution <- cbind.data.frame(
  id = rownames(Est.skin.bulk$Est.prop.weighted.cluster),
  Est.skin.bulk$Est.prop.weighted.cluster
)

pseudo_bulk <- read.table(
  "pseudobulk_scRNA_combined_2000.txt",
  header = TRUE,
  sep = "\t"
)

rownames(pseudo_bulk) <- pseudo_bulk[, 1]
pseudo_bulk <- pseudo_bulk[, -1]
pseudo_bulk <- as.matrix(pseudo_bulk)

Est.pseudo.bulk = music_prop.cluster(
  bulk.mtx = pseudo_bulk,
  sc.sce = counts,
  group.markers = markers,
  clusters = 'cell_type',
  groups = 'clusterType',
  samples = 'donor_id',
  clusters.type = groups
)
pseudo_deconvolution <- cbind.data.frame(
  id = rownames(Est.pseudo.bulk$Est.prop.weighted.cluster),
  Est.pseudo.bulk$Est.prop.weighted.cluster
)

write.xlsx(deconvolution, "MuSiC_deconvolution.xlsx")
write.xlsx(pseudo_deconvolution, "MuSiC_pseudo_deconvolution.xlsx")
