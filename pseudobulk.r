headers <- scan(
  "output_5c9ab5a5-04a9-4282-9320-3b4d7b95131c/CIBERSORTx_Adjusted.txt",
  sep = "\t",
  what = "character",
  nlines = 1
)
pseudo_bulk <- read.table(
  "output_5c9ab5a5-04a9-4282-9320-3b4d7b95131c/CIBERSORTx_Adjusted.txt",
  sep = "\t",
  skip = 1
)

colnames(pseudo_bulk) <- headers

headers1 <- scan(
  "output_26f6625b-e76c-490a-beb1-aea16933cd6d/CIBERSORTx_Adjusted.txt",
  sep = "\t",
  what = "character",
  nlines = 1
)
pseudo_bulk1 <- read.table(
  "output_26f6625b-e76c-490a-beb1-aea16933cd6d/CIBERSORTx_Adjusted.txt",
  sep = "\t",
  skip = 1
)

colnames(pseudo_bulk1) <- headers1

headers2 <- scan(
  "output_4d3469a7-339f-40b3-92a3-22f7043545f8/CIBERSORTx_Adjusted.txt",
  sep = "\t",
  what = "character",
  nlines = 1
)

pseudo_bulk2 <- read.table(
  "output_4d3469a7-339f-40b3-92a3-22f7043545f8/CIBERSORTx_Adjusted.txt",
  sep = "\t",
  skip = 1
)

colnames(pseudo_bulk2) <- headers2

View(pseudo_bulk)
View(pseudo_bulk1)
View(pseudo_bulk2)
