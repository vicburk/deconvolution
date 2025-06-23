header <- scan(
  "output_5c9ab5a5-04a9-4282-9320-3b4d7b95131c/CIBERSORTx_Adjusted.txt",
  sep = "\t",
  what = "character",
  nlines = 1
)

estimated_proportions <- read.table(
  "output_5c9ab5a5-04a9-4282-9320-3b4d7b95131c/CIBERSORTx_Adjusted.txt",
  sep = "\t",
  skip = 1
)

colnames(estimated_proportions) <- header

estimated_proportions
