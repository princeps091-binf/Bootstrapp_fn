library(renv)

renv::init(bioconductor = TRUE)
renv::install("tidyverse")
renv::install("furrr")
renv::install("valr")
renv::install("vroom")
renv::install("parallel")
renv::install("data.tree")
renv::install("igraph")
renv::install("crayon")
renv::install("bedr")
renv::install("UpSetR")
renv::install("svglite")
renv::install("bedr")
renv::install("pkgload")

renv::install("PhanstielLab/bedtoolsr")

renv::install("bioc::GenomicRanges")
renv::install("bioc::TxDb.Hsapiens.UCSC.hg19.knownGene")
renv::install("bioc::ChIPseeker")
renv::install("bioc::org.Hs.eg.db")
renv::install("bioc::rtracklayer")


renv::snapshot()
