if (! file.exists("processed_data/dictionary.xlsx")) {
  source("dictionary.R")
  name <- "processed_data/gene_expression_clean.xlsx"
  gene_expression <- read.xlsx(name, rowNames = TRUE) |> as.matrix()
  gene_meta <- read.xlsx("processed_data/gene_meta.xlsx")
  dictionary <- read.xlsx("processed_data/dictionary.xlsx")
} else {
  source("packages.R")
  name <- "processed_data/gene_expression_clean.xlsx"
  gene_expression <- read.xlsx(name, rowNames = TRUE) |> as.matrix()
  gene_meta <- read.xlsx("processed_data/gene_meta.xlsx")
  dictionary <- read.xlsx("processed_data/dictionary.xlsx")
}

probes <- rownames(gene_expression)

# Filter probes to probes with known gene names.
index <- probes %in% dictionary$PROBE
gene_expression <- gene_expression[index, ]


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

gene_exp <- translate(gene_expression)
gene_exp <- gene_exp$data

rownames(gene_exp) <- gsub("[.]", "-", rownames(gene_exp))

dim(gene_exp)
dim(na.omit(gene_exp))

gene_exp <- na.omit(gene_exp)

rm(list = c("gene_expression",
            "gene_meta",
            "index",
            "name"))

gc()

write.table(cbind(Gene = rownames(gene_exp), gene_exp),
            file = "gene_expression.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)
