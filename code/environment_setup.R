library(renv)

renv::init(bioconductor = TRUE)
renv::install("tidyverse")
renv::install("furrr")
renv::install("valr")
renv::install("vroom")
renv::install("parallel")
renv::install("data.tree")
renv::install("igraph")

renv::install("bioc::GenomicRanges")
renv::install("bioc::TxDb.Hsapiens.UCSC.hg19.knownGene")
renv::install("bioc::ChIPseeker")
renv::install("bioc::org.Hs.eg.db")

renv::snapshot()