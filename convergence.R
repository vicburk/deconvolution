#####################################
# Load cibersortx fraction estimates
#####################################

files <- list.files(pattern = "output*")

files <- paste0(files,"/CIBERSORTx_Adjusted.txt")

data  <- lapply(
  files,
  read.table,
  sep = "\t",
  skip = 1
)

column_names  <- lapply(
  files,
  scan,
  sep = "\t",
  what = "character",
  nlines = 1
)

row_names <- lapply(seq_along(data), function(x) data[[x]][, 1])

######################################
# Clean data
######################################

for (i in seq_along(data)) {
  names(data[[i]]) <- column_names[[i]]
  n <- ncol(data[[i]])
  data[[i]] <- data[[i]][, -seq(n-2, n)]
  data[[i]] <- data[[i]][, -1]
  nonzero <- colSums(data[[i]]) != 0
  data[[i]] <- data[[i]][, nonzero]
}

# Get cell types that are in all datasets

common <- Reduce(intersect, lapply(data,colnames))

for (i in seq_along(data)) {
  data[[i]] <- subset(data[[i]], select = common)
}

#######################################
# Convergence analysis
#######################################

corr <- list()
for (i in seq_along(data)[- length(data)]) {
  flatten0 <- unlist(data[[i]])
  flatten1 <- unlist(data[[i + 1]])
  corr[[i]] <- cor(flatten0, flatten1, method = "spearman")
}

patient_wise <- list()
for (i in seq_along(data)[-length(data)]) {
  corrs <- NULL
  for (j in seq_len(nrow(data[[i]]))) {
    v0 <- unlist(data[[i]][j, ])
    v1 <- unlist(data[[i + 1]][j, ])
    coR <- cor(v0, v1, method = "spearman")
    corrs <- append(corrs, coR)
  }
  patient_wise[[i]] <- cbind(
    mean = mean(corrs),
    std = sd(corrs)
  )
}

cell_wise <- list()
for (i in seq_along(data)[-length(data)]) {
  corrs <- NULL
  for (j in seq_len(ncol(data[[i]]))) {
    v0 <- unlist(data[[i]][, j])
    v1 <- unlist(data[[i + 1]][, j])
    coR <- cor(v0, v1, method = "spearman")
    corrs <- append(corrs, coR)
  }
  cell_wise[[i]] <- cbind(
    mean = mean(corrs),
    std = sd(corrs)
  )
}
