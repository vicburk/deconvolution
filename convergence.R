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

# Extract numeric values from files

values <- gsub("[^0-9]", "", files)
values <- as.numeric(values)

names(data) <- values

names(column_names) <- values

row_names <- lapply(seq_along(data), function(x) data[[x]][, 1])

######################################
# Clean data
######################################

neworder <- data %>%
  names %>%
  as.numeric() %>%
  sort() %>%
  as.character()

for (i in seq_along(neworder)) {
  index <- neworder[i]
  names(data[index][[1]]) <- column_names[index][[1]]
  n <- ncol(data[index][[1]])
  data[index][[1]] <- data[index][[1]][, -seq(n-2, n)]
  data[index][[1]] <- data[index][[1]][, -1]
  nonzero <- colSums(data[index][[1]]) != 0
  data[index][[1]] <- data[index][[1]][, nonzero]
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
  flatten0 <- unlist(data[neworder[i]][[1]])
  flatten1 <- unlist(data[neworder[i + 1]][[1]])
  title <- paste(neworder[i], neworder[i + 1],sep = "-")
  corr[title] <- cor(flatten0, flatten1, method = "spearman")
}

patient_wise <- list()
for (i in seq_along(data)[-length(data)]) {
  corrs <- NULL
  for (j in seq_len(nrow(data[[i]]))) {
    v0 <- unlist(data[neworder[i]][[1]][j, ])
    v1 <- unlist(data[neworder[i + 1]][[1]][j, ])
    coR <- cor(v0, v1, method = "spearman")
    corrs <- append(corrs, coR)
  }
  title <- paste(neworder[i], neworder[i + 1],sep = "-")
  patient_wise[[title]] <- cbind(
    mean = mean(corrs),
    std = sd(corrs)
  )
}

cell_wise <- list()
for (i in seq_along(data)[-length(data)]) {
  corrs <- NULL
  for (j in seq_len(ncol(data[[i]]))) {
    v0 <- unlist(data[neworder[i]][[1]][, j])
    v1 <- unlist(data[neworder[i + 1]][[1]][, j])
    coR <- cor(v0, v1, method = "spearman")
    corrs <- append(corrs, coR)
  }
  title <- paste(neworder[i], neworder[i + 1],sep = "-") 
  cell_wise[[title]] <- cbind(
    mean = mean(corrs),
    std = sd(corrs)
  )
}
