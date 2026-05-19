#!/usr/bin/env python3

import re
import sys
import numpy as np
import statsmodels.api as sm
from statsmodels.stats.mediation import Mediation

# python modules/local/mediation/resources/usr/bin/mediation.py roi_gm.tsv ../NemaScan/input_data/c_elegans/phenotypes/expression/tx5291exp_st207.tsv pheno.tsv V 2717620 ../NemaScan/input_data/c_elegans/phenotypes/expression/eQTL6545forMed.tsv C08E3.11 testvun

# Load arguments
args = sys.argv[1:]
genotype_matrix_fname = args[0]
texpression_pheno_fname = args[1]
phenotype_fname = args[2]
gwas_chr = args[3]
gwas_peak = int(args[4])
eqtl_fname = args[5]
gene_name = args[6]
out_fname = args[7]

# Load genotype matrix
fs = open(genotype_matrix_fname)
header = fs.readline().rstrip().split('\t')
gt_strains = np.array(header[4:], 'U8')
strain_set = set(header[4:])
for line in fs:
    line = line.rstrip().split('\t')
    if int(line[1]) == gwas_peak:
        genotype = np.array([(gt_strains[i], x) if x != "NA" else (gt_strains[i], np.nan) for i, x in enumerate(line[4:])],
                            dtype=np.dtype([('strain', 'U8'), ('gt', np.float32)]))
fs.close()
genotype = genotype[np.argsort(genotype['strain'])]

# Load pheno data
phenotype = np.loadtxt(phenotype_fname, skiprows=1, dtype=np.dtype([('strain', 'U8'), ('value', np.float32)]))
phenotype['value'] -= np.nanmean(phenotype['value'])
phenotype['value'] /= np.nanstd(phenotype['value'], ddof=1)
phenotype = phenotype[np.array([i for i in range(phenotype.shape[0])
                                            if str(phenotype['strain'][i]) in strain_set])]
phenotype = phenotype[np.argsort(phenotype['strain'])]

# Transcript level
fs = open(texpression_pheno_fname)
header = fs.readline().rstrip().split()
index = header.index(gene_name)
expression = []
for line in fs:
    line = line.rstrip().split("\t")
    if line[0] not in strain_set:
        continue
    if line[index] == 'NA':
        line[index] = np.nan
    expression.append((line[0], line[index]))
expression = np.array(expression, dtype=np.dtype([('strain', 'U8'), ('exp', np.float32)]))
expression['exp'] -= np.nanmean(expression['exp'])
expression['exp'] /= np.nanstd(expression['exp'])
expression = expression[np.argsort(expression['strain'])]
strain_set = set([str(strain) for strain in expression['strain']])

# Filter genotypes by valid strains
gt_valid_strains = np.array([i for i in range(gt_strains.shape[0]) if str(gt_strains[i]) in strain_set])
genotype = genotype[gt_valid_strains]
ph_valid_strains = np.array([i for i in range(phenotype.shape[0]) if str(phenotype['strain'][i]) in strain_set])
phenotype = phenotype[ph_valid_strains]

data = np.concatenate((genotype['gt'].reshape(-1, 1),
                       expression['exp'].reshape(-1, 1)), axis=1)

mediator_model = sm.OLS(expression['exp'], genotype['gt'], missing='drop')
outcome_model = sm.OLS(phenotype['value'], data, missing='drop')
med = Mediation(outcome_model, mediator_model, [0, 0], 1).fit()
results = med.summary()
results.to_csv(out_fname, sep="\t", index=True)
