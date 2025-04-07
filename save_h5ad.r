##############################################################

options(timeout = 300)

load_data <- function(url_path, file) {
  if (!file.exists(file)) {
    url <- paste0(url_path, file)
    download.file(url = url, destfile = file, mode = "wb")
    print("File downloaded")
  } else {
    print("File already exists")
  }
}

file0 <- "4d3469a7-339f-40b3-92a3-22f7043545f8.h5ad"
file1 <- "26f6625b-e76c-490a-beb1-aea16933cd6d.h5ad"
file2 <- "5c9ab5a5-04a9-4282-9320-3b4d7b95131c.h5ad"
url_path <- "https://datasets.cellxgene.cziscience.com/"

load_data(url_path, file0)
load_data(url_path, file1)
load_data(url_path, file2)