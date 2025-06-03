# If using ubuntu, install r2u to download R package binaries.

if (!require("devtools")) {
  install.packages("devtools")
}

using <- function(packages, versions) {
  libs <- packages[, 1]
  req <- unlist(lapply(libs, require, character.only = TRUE))
  need <- packages[req == FALSE, ]
  need <- t(as.matrix(need)) 
  if (length(need)) {
    versions <- need[, 2]
    need <- need[, 1]
    lapply(
      seq_along(need),
      function(x) {
        install_version(need[x], versions[x])
      }
    )
    lapply(need, require, character.only = TRUE)
  }
}

packages <- rbind(
  c("BiocManager", "1.30.25"),
  c("readxl", "1.4.3"),
  c("openxlsx", "4.2.8"),
  c("dplyr", "1.1.4"),
  c("Seurat", "5.2.1")
)

using(packages)

using_bioconductor <- function(packages, versions) {
  libs <- packages[, 1]
  req <- unlist(lapply(libs, require, character.only = TRUE))
  need <- packages[req == FALSE, ]
  need <- matrix(need, ncol = 2, byrow = TRUE)
  if (length(need)) {
    versions <- need[, 2]
    need <- need[, 1]
    lapply(
      seq_along(need),
      function(x) {
        BiocManager::install(
          pkgs = need[x]
        )
      }
    )
    lapply(need, require, character.only = TRUE, force = TRUE)
  }
}

bioconductor <- rbind(
  c("HDO.db", "1.0.0"),
  c("biomaRt", "2.62.0"),
  c("hgu219.db", "3.2.3"),
  c("Homo.sapiens", "1.3.1"),
  c("SeuratObject", "5.0.2"),
  c("EnsDb.Hsapiens.v79", "2.99.0")
)

using_bioconductor(bioconductor)


