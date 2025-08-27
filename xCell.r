# BiocManager::install(c("GSVA","GSEABase"))

# devtools::install_github('dviraran/xCell') 

# Load packages

library(xCell)
library(dplyr)
library(openxlsx)

#################################

# Load data

data <- "gene_expression.txt"
expr <- read.table(data, header = TRUE, row.names = 1, as.is = TRUE, sep = "\t")

pseudo <- read.table("pseudobulk_scRNA_combined_2000.txt", header = TRUE, sep = "\t")
rownames(pseudo) <- pseudo[, 1]
pseudo <- pseudo[, -1]

pseudo2 <- read.table("pseudobulk_scRNA_combined2_2000.txt", header = TRUE, sep = "\t")
rownames(pseudo2) <- pseudo2[, 1]
pseudo2 <- pseudo2[, -1]

###############################
# Select cell types
###############################

all_cell_types <- xCell.data$spill$K %>% colnames

all_cell_types

cell_types <- c(
  "Fibroblasts",
  "Keratinocytes",
  "Pericytes",
  "Melanocytes",
  "Tregs",
  "DC",
  "pDC",
  "Monocytes",
  "Smooth muscle",
  "Skeletal muscle",
  "Chondrocytes",
  "Endothelial cells",
  "NK cells"
)

cds <- grep("CD.*", all_cell_types, value = TRUE)
macrophages <- grep("Macrophages.*", all_cell_types, value = TRUE)
# bcells <- grep(".*B-cells", all_cell_types, value = TRUE)
th <- grep("Th.*", all_cell_types, value = TRUE)

cell_types <- c(cell_types, cds, macrophages, th)

cell_types

#####################################
##### xCell analysis
#####################################

xCell_result_scaled <- xCellAnalysis(expr, rnaseq = FALSE, cell.types.use = cell_types,
                                     scale = TRUE)


xCell_result <- xCellAnalysis(expr, rnaseq = FALSE, cell.types.use = cell_types,
                              scale = FALSE)

######################################
### Redefine cell types
######################################

xCell_pipeline <- function(expr, cell_type, scale = TRUE) {
  if (scale == TRUE) {
    xCell_result <- xCellAnalysis(expr,
                                  rnaseq = FALSE,
                                  cell.types.use = cell_type,
                                  scale = TRUE)
  } else {
    xCell_result <- xCellAnalysis(expr,
                                  rnaseq = FALSE,
                                  cell.types.use = cell_type,
                                  scale = FALSE)
  }
  helper <- c(
    "CD4+ T-cells",
    "CD4+ naive T-cells",
    "Th1 cells",
    "Th2 cells"
  )

  cytotoxic <- c(
    "CD8+ T-cells",
    "CD8+ naive T-cells"
  )

  memory <- c(
    "CD4+ memory T-cells",
    "CD4+ Tcm", # (Central memory T cells)
    "CD4+ Tem", # (Effector memory T cells)
    "CD8+ Tcm", # (CD8+ Central memory T cells)
    "CD8+ Tem" # (CD8+ Effector memory T cells)
  )

  dc <- c(
    "DC",
    "pDC"
  )

  # Sum rows for helper T cells
  helper_mean <- apply(xCell_result[helper, ], 2, mean)

  # Sum rows for cytotoxic T cells
  cytotoxic_mean <- apply(xCell_result[cytotoxic, ], 2, mean)

  # Sum rows for memory T cells
  memory_mean <- apply(xCell_result[memory, ], 2, mean)

  # Sum rows for Macrophages
  macrophages_mean <- apply(xCell_result[macrophages, ], 2, mean)

  # Sum rows for Dendritic cells
  dendritic_mean <- apply(xCell_result[dc, ], 2, mean)

  # Create a new data frame with the summed rows
  index <- ! rownames(xCell_result) %in% c(helper,
                                           cytotoxic,
                                           memory,
                                           macrophages,
                                           dc)
  xCell_result_new <- xCell_result[index, ]
  rownames(xCell_result_new) <- tolower(rownames(xCell_result_new))
  rownames(xCell_result_new) <- gsub("s$", "", rownames(xCell_result_new))

  xCell_result_new <- rbind(xCell_result_new,
                            helper.T.cell = helper_mean,
                            cytotoxic.T.cell = cytotoxic_mean,
                            memory.T.cell = memory_mean,
                            macrophage = macrophages_mean,
                            dendritic.cell = dendritic_mean)

  xCell_result_new <- t(xCell_result_new)

  xCell_result_new <- xCell_result_new %>% as.data.frame

  return(xCell_result_new)
}

###############

xCell_result <- xCell_pipeline(expr = expr, cell_type = cell_types, scale = FALSE)
t(xCell_result/rowSums(xCell_result))


xCell_pseudo_result <- xCell_pipeline(expr = pseudo, cell_type = cell_types, scale = FALSE)
xCell_pseudo_result/rowSums(xCell_pseudo_result)

xCell_pseudo2_result <- xCell_pipeline(expr = pseudo2, cell_type = cell_types, scale = FALSE)
t(xCell_pseudo2_result/rowSums(xCell_pseudo2_result))

xCell_pseudo3_result <- xCellAnalysis(
  expr = pseudo2,
  cell.types.use = c(
    "Keratinocytes",
    "Endothelial cells",
    "Fibroblasts"
  ),
  scale = FALSE
)
xCell_pseudo3_result <- t(xCell_pseudo3_result)
t(xCell_pseudo3_result/rowSums(xCell_pseudo3_result))


xCell_pseudo3_result_all <- xCellAnalysis(
  expr = pseudo2,
  scale = FALSE
)
xCell_pseudo3_result_all <- xCell_pseudo3_result_all[-seq(nrow(xCell_pseudo3_result_all)-2,nrow(xCell_pseudo3_result_all)),]
xCell_pseudo3_result_all <- t(xCell_pseudo3_result_all)
t(xCell_pseudo3_result_all/rowSums(xCell_pseudo3_result_all))



