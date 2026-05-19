#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(sommer)

# args <- c("Genotype_Matrix.tsv",
#           "Phenotype_data.tsv",
#           "output.txt"
#          )

args <- commandArgs(trailingOnly = TRUE)

# load genotype matrix
genotype_matrix <- readr::read_tsv(args[1])

# load phenotpe data
phenotype_data <- data.table::fread(args[2]) %>%
  na.omit() %>%
  as.data.frame()
names(phenotype_data) <- c("strain", "value")

A <- sommer::A.mat(t(genotype_matrix %>% dplyr::select(-CHROM, -REF, -ALT, -POS)))
h2_res <- sommer::mmes(value ~ 1, random = ~sommer::vsm(sommer::ism(strain), Gu = A), data = phenotype_data)
h2 <- as.numeric(sommer::vpredict(h2_res, h2 ~ (V1) / (V1+V2))[[1]][1])

writeLines(paste(h2), con=args[3])
