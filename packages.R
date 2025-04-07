using <- function(packages, versions) {
  libs <- packages[, 1]
  req <- unlist(lapply(libs, require, character.only = TRUE))
  need <- packages[req == FALSE, ]
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
  c("dplyr", "1.1.4")
)

using(packages)

using_bioconductor <- function(packages, versions) {
  libs <- packages[, 1]
  req <- unlist(lapply(libs, require, character.only = TRUE))
  need <- packages[req == FALSE, ]
  if (length(need)) {
    versions <- need[, 2]
    need <- need[, 1]
    lapply(
      seq_along(need),
      function(x) {
        BiocManager::install(
          pkgs = need[x],
          version = versions[x]
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
  c("Seurat", "5.2.1"),
  c("EnsDb.Hsapiens.v79", "2.99.0")
)

using_bioconductor(bioconductor)
