headers <- scan(
  "pseudo_output2_2000/CIBERSORTx_Results.txt",
  sep = "\t",
  what = "character",
  nlines = 1
)

pseudo_bulk <- read.table(
  "pseudo_output2_2000/CIBERSORTx_Results.txt",
  sep = "\t",
  skip = 1
)

names(pseudo_bulk) <- headers

remove <- c(
  1,
  seq(
    ncol(pseudo_bulk) - 2,
    ncol(pseudo_bulk)
  )
)

pseudo_bulk <- pseudo_bulk[, -remove]
pseudo_bulk <- unlist(pseudo_bulk)


estimated_proportions <- cbind.data.frame(
  cell_types = names(pseudo_bulk),
  estimated_proportions = pseudo_bulk
)

rownames(estimated_proportions) <- NULL

true_proportions <- read.table(
  "true_proportions_2000.txt",
  sep = "\t"
)

names(true_proportions) <- c(
  "cell_type",
  "true_proportion"
)

estimated_proportions
true_proportions

estimated_proportions <- estimated_proportions[match(
  true_proportions$cell_type,
  estimated_proportions$cell_types
), ]

proportions  <- cbind.data.frame(
  true_proportions,
  estimated_proportions = estimated_proportions$estimated_proportions
)

proportions <- proportions[order(
  proportions$true_proportion,
  decreasing = TRUE
), ]
proportions
