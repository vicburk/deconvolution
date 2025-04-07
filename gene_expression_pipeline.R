source("packages.R")

name <- "processed_data/subset_gene_expression.xlsx"
subset_gene_expression <- read.xlsx(name) |> as.matrix()
subset_gene_meta <- read.xlsx("processed_data/subset_gene_meta.xlsx")
dictionary <- read.xlsx("processed_data/dictionary.xlsx")

# Note, subject ID 54 is not in the protein dataset.
rownames(subset_gene_expression) <- subset_gene_meta$SUBJECT_ID
index <- which(rownames(subset_gene_expression) != "54")
subset_gene_expression <- subset_gene_expression[index, ]
subset_gene_meta <- subset_gene_meta[subset_gene_meta$SUBJECT_ID != "54", ]

# Filter probes to probes with known gene names.
index <- colnames(subset_gene_expression) %in% dictionary$PROBE
subset_gene_expression <- subset_gene_expression[, index]

# Translate probe ID to gene symbol.
translate <- function(data) {
  probes <- rownames(data)
  dict <- dictionary[dictionary$PROBE %in% probes, ]
  symbol <- dict[match(probes, dict$PROBE), 2]
  entrez <- dict[match(probes, dict$PROBE), 3]
  dup <- unique(symbol[duplicated(symbol)])
  clean <- matrix(nrow = length(dup),
                  ncol = ncol(data))
  colnames(clean) <- colnames(data)
  for (i in seq_along(dup)) {
    clean[i, ] <- apply(data[which(symbol == dup[i]), ], 2, mean)
  }
  rownames(clean) <- dup
  temp <- data[!symbol %in% dup,]
  rownames(temp) <- symbol[!symbol %in% dup]
  clean <- rbind(clean, temp)
  return(list(data = clean,
              symbol = symbol,
              entrez = entrez))
}

gene_exp <- translate(gene_expression_clean)
gene_exp <- gene_exp$data

rownames(gene_exp) <- gsub("[.]", "-", rownames(gene_exp))

dim(gene_exp)
dim(na.omit(gene_exp))

gene_exp <- na.omit(gene_exp)

rm(list = c("subset_gene_expression",
            "subset_gene_meta",
            "index",
            "name"))

gc()