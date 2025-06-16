pseudo_bulk <- read.table(
  "pseudo_scRNA_combined_2000.txt",
  sep = "\t",
  header = TRUE,
  skip = 1,
  row.names = 1)

pseudo_bulk <- pseudo_bulk %>% as.matrix
pseudo_columns <- scan(
  "pseudo_scRNA_combined_2000.txt",
  what = "character",
  nlines = 1,
  sep = "\t"
)

pseudo_columns <- pseudo_columns[-1]
colnames(pseudo_bulk) <- pseudo_columns

true_proportions <- pseudo_bulk %>%
  colnames %>%
  table %>%
  prop.table

cells <- pseudo_columns %>%
  unique

pseudo <- list()
for (i in seq_along(cells)) {
  index <- which(pseudo_columns == cells[i])
  temp <- pseudo_bulk[, index]
  temp <- rowSums(temp)
  pseudo[[i]] <- temp
}
pseudo <- Reduce(cbind, pseudo)
colnames(pseudo) <- cells

pseudo <- pseudo/rowSums(pseudo)*1e6
pseudo <- log2(pseudo + 1)

main_cells <- c("fibroblast",
                "helper T cell",
                "keratinocyte",
                "pericyte",
                "regulatory T cell")

proportions <- diag(5)
proportions <- cbind(proportions,
                     c(0.5,0,0.5,0,0),
                     c(0.5,0.5,0,0,0),
                     c(0,0,0,0.5,0.5),
                     c(0.01,0,1,0,0),
                     c(0.05,0,0,1,0),
                     c(0.1,1,0,0,0),
                     c(0.15,0,0,0,1))

total_prop <- proportions/sum(proportions)

n_samp <- round(apply(total_prop/rowSums(total_prop),2,function(x) x * head(sort(table(cells),decreasing = TRUE),5)))

n_samp |> rowSums()
cells |> table() |> sort(decreasing = TRUE) |> head(5)

n_samp |> sum()

sum(!cells %in% main_cells)/ncol(proportions)
rbind(n_samp,rep(round(sum(!cells %in% main_cells)/ncol(proportions)),ncol(proportions)))
prop.table(rbind(n_samp,rep(round(sum(!cells %in% main_cells)/ncol(proportions)),ncol(proportions))),2)

rbind(n_samp,rep(round(sum(!cells %in% main_cells)/ncol(proportions)),ncol(proportions))) |> rowSums()


cells1 <- cells
cell_type_list <-  list()
set.seed(42)
for (i in seq_along(main_cells)) {
  master_index <- which(cells1 %in% main_cells[i])
  sample_list <- list()
  index <- seq_along(master_index)
  for (j in 1:ncol(n_samp)) {
    samp <- sample(index, n_samp[i,j])
    sample_list[[j]] <- master_index[samp]
    index <- index[!index %in% samp]
  }
  cell_type_list[[i]] <- sample_list
}

names(cell_type_list) <- main_cells

other_cells <- unique(cells)[! unique(cells) %in% main_cells]

##

master_index <- which(!cells1 %in% main_cells)
index <- seq_along(master_index)
samp_n <- floor(length(master_index)/ncol(proportions))
samp_l <- list()
bulk_cell_types <- list()
set.seed(42)
for (i in 1:ncol(n_samp)) {
  samp <- sample(index, samp_n)
  samp_l[[i]] <- master_index[samp]
  bulk_cell_types[[i]] <- cells[master_index[samp]]
  index <- index[!index %in% samp]
}

total_counts_other_cells <- lapply(bulk_cell_types,table)

counts_other_cells <- NULL
for (j in seq_along(samp_l)) {
  if (length(samp_l[[j]] > 0)) {
    a <- split(samp_l[[j]],samp_l[[j]] <= ncol(counts))
    a1 <- unlist(a$`FALSE`)
    a1 <- a1 - ncol(counts)
    a2 <- unlist(a$`TRUE`)
    data_a1 <- counts_other[,a1]
    data_a2 <- counts[,a2]
    data <- cbind(data_a1,data_a2)
    sample_counts <- rowSums(data)
  } else {
    sample_counts <- rep(0,nrow(counts))
  }
  counts_other_cells <- cbind(counts_other_cells,sample_counts) 
}
counts_other_cells

##

counts_by_cell_type <- list()
i_j_counts <- NULL
for (i in seq_along(cell_type_list)) {
  counts_matrix <- NULL
  for (j in seq_along(cell_type_list[[i]])) {
    if (length(cell_type_list[[i]][[j]] > 0)) {
      a <- split(cell_type_list[[i]][[j]],cell_type_list[[i]][[j]] <= ncol(counts))
      a1 <- unlist(a$`FALSE`)
      a1 <- a1 - ncol(counts)
      a2 <- unlist(a$`TRUE`)
      data_a1 <- counts_other[,a1]
      data_a2 <- counts[,a2]
      data <- cbind(data_a1,data_a2)
      sample_counts <- rowSums(data)
    } else {
      sample_counts <- rep(0,nrow(counts))
    }
    counts_matrix <- cbind(counts_matrix,sample_counts)  
    title <- paste(main_cells[i],j)
    temp <- c(length(cell_type_list[[i]][[j]]))
    names(temp) <- title
    i_j_counts <- append(i_j_counts, temp)
  }
  counts_by_cell_type[[i]] <- counts_matrix
}

counts_by_cell_type

counts_by_cell_type[[length(counts_by_cell_type) + 1]] <- counts_other_cells

counts_by_cell_type


pseudo_bulk <- Reduce(function(x,y) x + y, counts_by_cell_type)

pseudo_bulk <- pseudo_bulk/colSums(pseudo_bulk)*1e6

pseudo_bulk <- log2(pseudo_bulk + 1)

count_table <-  matrix(i_j_counts,ncol = ncol(proportions),byrow = TRUE)
temp <- Reduce(cbind,total_counts_other_cells)

count_table <- rbind(count_table,temp)

colnames(count_table) <- NULL
rownames(count_table) <- c(main_cells,rownames(temp))

count_table
actual_proportions <- prop.table(count_table,2)

colSums(prop.table(count_table,2)[-(1:5),])

prop.table(rbind(n_samp,rep(round(sum(!cells %in% main_cells)/ncol(proportions)),ncol(proportions))),2)

colnames(pseudo_bulk) <- 1:ncol(proportions)

pseudo_bulk <- cbind.data.frame(GeneSymbol = rownames(pseudo_bulk),
                                pseudo_bulk)


