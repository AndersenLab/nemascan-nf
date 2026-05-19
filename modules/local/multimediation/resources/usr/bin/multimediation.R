#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(tibble)
library(MultiMed)

# load arguments
# 1 - genotype_matrix
# 2 - transcript_expression
# 3 - subset eQTL
# 4 - normalized phenotypes
# 5 - chrom
# 6 - peakPOS
# 7 - output fname

args <- commandArgs(trailingOnly = TRUE)

# GWAS qtl infor
gwas_chr = args[5]
gwas_peak = args[6] %>% as.numeric()

# load genotype matrix
Genotype_Matrix <- readr::read_tsv(args[1])

# load pheno data
trait_phenotype <- readr::read_tsv(args[4])
names(trait_phenotype) = c('strain', 'trait')

# load eqtl data
eqtl_infor <- read.delim(args[3], stringsAsFactors=FALSE)

transcript_list <- eqtl_infor %>% 
  dplyr::select(trait) %>% 
  dplyr::distinct()

# transcript level
texpression_pheno_raw <- data.table::fread(args[2])

texpression_pheno <- texpression_pheno_raw %>% 
  tidyr::gather(trait2,value,-strain) %>% 
  dplyr::mutate(trait=sub("(^X)(.*)","\\2",trait2)) %>% 
  dplyr::select(strain,trait,value) %>% 
  dplyr::filter(strain %in% trait_phenotype$strain) %>% 
  dplyr::filter(trait %in% transcript_list$trait) %>% 
  na.omit() 

# processed pheno data
trait_pheno_all <- trait_phenotype%>% 
  dplyr::filter(strain %in% texpression_pheno$strain)

# get the genotype at the peak marker
gwas_g_all <- Genotype_Matrix %>% 
  dplyr::filter(CHROM==gwas_chr & POS == gwas_peak) %>% 
  dplyr::select(-(1:4)) %>% 
  tidyr::gather(strain,geno) %>% 
  dplyr::filter(strain %in% texpression_pheno$strain) %>% 
  dplyr::arrange(strain) %>%
  na.omit()

# step through genes and test for significant mediation
multimed_trait_list=list()
for(trss in unique(texpression_pheno$trait)) {  
  # get gene expression for strains with phenotype data
  t_lgmtpm_gwas_all <- texpression_pheno %>% 
    dplyr::filter(trait==trss,
                  strain %in% gwas_g_all$strain)
  
  # reformat expression
  t_lgmtpm_gwas <- t_lgmtpm_gwas_all %>% 
    tidyr::spread(trait,value) %>% 
    dplyr::arrange(strain)  %>% 
    dplyr::select(-strain) 
  
  # limit genotypes and phenotypes to strains in expression data
  gwas_g <- gwas_g_all %>% 
    dplyr::filter(strain %in% t_lgmtpm_gwas_all$strain)
   
  trait_pheno <- trait_pheno_all%>% 
    dplyr::filter(strain %in% t_lgmtpm_gwas_all$strain)
  
  # if there are enough values, continue with test  
  if( length(unique(gwas_g$geno)) == 2 & length(unique(t_lgmtpm_gwas_all$value)) > 1 ){
    
    # pick transcripts with variation in expression
    t_lgmtpm_gwas_vari <- t_lgmtpm_gwas[vapply(t_lgmtpm_gwas, function(x) length(unique(x)) > 1, logical(1L))]
    exp_matr_transcript <- as.matrix(t_lgmtpm_gwas_vari)
    
    # perform mediation test
    mt_multi_transcript <- medTest(gwas_g$geno, exp_matr_transcript, trait_pheno$trait, nperm = 1000)
    df_multi_transcript <- as.data.frame(mt_multi_transcript) 
    row.names(df_multi_transcript) <- (colnames(exp_matr_transcript))
    df_multi_transcript2 <- df_multi_transcript %>% 
      tibble::rownames_to_column(var="gene") 
    
    # add results to list
    multimed_trait_list[[trss]] <- df_multi_transcript2
  }  
}

df_multi_med <- dplyr::bind_rows(multimed_trait_list)  %>% 
  dplyr::arrange(p)

# Join original columns with mediation test values
df_multi_med <- df_multi_med %>%
  dplyr::left_join(eqtl_infor,by=c("gene"="trait"))

readr::write_tsv(df_multi_med,
                 path = glue::glue("{args[7]}_medmulti.tsv"),
                 col_names = T)

# Filter results by P-value and 99th quantile
q99_S = quantile(df_multi_med$S, probs = 0.99)[[1]]

gene_qtl_genelist <- df_multi_med %>% 
  dplyr::filter(p<0.05 | S > q99_S) %>% 
  dplyr::select(gene) %>% 
  dplyr::distinct()

 # save gene list 
if(nrow(gene_qtl_genelist)>0){
  readr::write_tsv(gene_qtl_genelist, 
                   path = glue::glue("{args[7]}_genes.tsv"),
                   col_names = F)
}



