using = function(...) {
  libs = unlist(list(...))
  req = unlist(lapply(libs,require,character.only=TRUE))
  need = libs[req==FALSE]
  if(length(need)){ 
    install.packages(need)
    lapply(need,require,character.only=TRUE)
  }
}

packages <- c(
  "BiocManager",
  "readxl",
  "openxlsx",
  "dplyr"
)

using(packages)

using_bioconductor = function(...) {
  libs = unlist(list(...))
  req = unlist(lapply(libs,require,character.only=TRUE))
  need = libs[req==FALSE]
  if(length(need)){ 
    BiocManager::install(need)
    lapply(need,require,character.only=TRUE,force=TRUE)
  }
}

bioconductor <- c("HDO.db",
                 "biomaRt",
                 "hgu219.db",
                 "Homo.sapiens",
                 "SeuratObject",
                 "Seurat",
                 "EnsDb.Hsapiens.v79")

using_bioconductor(bioconductor)
