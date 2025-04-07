source("packages.R")

##################################
# Import gene expression metadata
##################################

gene_expression_meta <- readxl::read_xlsx("data/BabyMeta_7_27_2021.xlsx",
                                          sheet = 1)
gene_expression_meta <- as.data.frame(gene_expression_meta)

##################################
# Import gene expression data
##################################

file <- unz("data/GSS2656_signal_noOutlier.zip", "GSS2656_signal_noOutlier.tsv")
gene_expression <- read.delim(file)
rownames(gene_expression) <- gene_expression$probeset_id

##################################

index0 <- which(gene_expression_meta$SAMP_ID %in% colnames(gene_expression))
gene_expression_meta <- gene_expression_meta[index0, ]

gene_expression <- gene_expression[, gene_expression_meta$SAMP_ID]

meta_list <- split(gene_expression_meta,
                   gene_expression_meta$`UseforALLBabies_vs_Adults20 years`)

adult_meta <- meta_list$Adult
baby_meta <- meta_list$Baby

baby_meta <- baby_meta[which(as.numeric(baby_meta$TA) < 26), ]
baby_meta <- baby_meta[baby_meta$SUBJECT_ID != "54", ]

meta_gene <- rbind(adult_meta, baby_meta)

adult_gene <- gene_expression[, adult_meta$SAMP_ID]
baby_gene <- gene_expression[, baby_meta$SAMP_ID]

colnames(adult_gene) <- adult_meta$SUBJECT_ID
colnames(baby_gene) <- baby_meta$SUBJECT_ID

gene_expression_clean <- cbind(adult_gene, baby_gene)

gene_expression_clean <- log2(gene_expression_clean)

write.xlsx(gene_expression_clean,
           "processed_data/gene_expression_clean.xlsx")