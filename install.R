if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("BSgenome",
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
"shinyccsloaders",
"Rmisc",
"VariantAnnotation"))
