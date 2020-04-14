packages <- c("BSgenome",
    "BSgenome.Hsapiens.UCSC.hg19",
    "devtools",
    "DT",
    "GenomicRanges",
    "ggplot2",
    "gbm",
    "gridExtra",
    "IRanges",
    "nnls",
    "pheatmap",
    "reshape2",
    "Rmisc",
    "VariantAnnotation")

if(getRversion() > "3.6.0"){
  if(grepl('3.6', getRversion())){
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(packages)
  }
}
if(getRversion() < "3.6.0" & getRversion() >= "3.5.0"){
  source("https://bioconductor.org/biocLite.R")
  biocLite(packages)
}

