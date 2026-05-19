#!/usr/bin/env python3

import re
import sys
import numpy as np
from medTest import medTest

# Load arguments
args = sys.argv[1:]
genotype_matrix_fname = args[0]
texpression_pheno_fname = args[1]
trait_phenotype_fname = args[2]
gwas_chr = args[3]
gwas_peak = int(args[4])
eqtl_fname = args[5]
out_prefix = args[6]

# Load genotype matrix
fs = open(genotype_matrix_fname)
header = fs.readline().rstrip().split('\t')
gm_strains = np.array(header[4:], 'U8')
strain_set = set(header[4:])
data = []
for line in fs:
    line = line.rstrip().split('\t')
    data.append(tuple(line[:4] + [tuple([x if x != "NA" else 0 for x in line[4:]])]))
fs.close()
genotype_matrix = np.array(data, dtype=np.dtype([('chr', f'U{max([len([0]) for x in data])}'), ('pos', np.int32), ('ref', 'U1'),
                                                 ('alt', 'U1'), ('genotype', np.int8, (gm_strains.shape[0],))]))
genotype_matrix = genotype_matrix[np.lexsort((genotype_matrix['pos'], genotype_matrix['chr']))]
gstart = genotype_matrix['pos'][0]
gstop = genotype_matrix['pos'][-1]

# Load pheno data
trait_phenotype = np.loadtxt(trait_phenotype_fname, skiprows=1, dtype=np.dtype([('strain', 'U8'), ('value', np.float32)]))

trait_phenotype['value'] -= np.nanmean(trait_phenotype['value'])
trait_phenotype['value'] /= np.nanstd(trait_phenotype['value'], ddof=1)
strain_set = strain_set.intersection(set([str(trait_phenotype['strain'][i]) for i in range(trait_phenotype.shape[0])]))

# Load eqtl data
fs = open(eqtl_fname)
header = fs.readline().rstrip().split('\t')
data = []
for i in range(len(header)):
    data.append([])
for line in fs:
    line = line.rstrip().split('\t')
    for i, item in enumerate(line):
        data[i].append(item)
fs.close()
dtypes = []
for i in range(len(data)):
    dtype_int = min([2 if re.fullmatch("-?[0-9]*\.[0-9]*", x) else int(x.isdigit()) for x in data[i]])
    if dtype_int == 0:
        dtype_val = f"U{max([len(x) for x in data[i]])}"
    elif dtype_int == 1:
        dtype_val = np.int32
    else:
        dtype_val = np.float32
    dtypes.append((header[i], dtype_val))
    data[i] = np.array(data[i], dtype=dtype_val)
eqtl_infor = np.zeros(data[0].shape[0], dtype=np.dtype(dtypes))
for i in range(len(data)):
    eqtl_infor[header[i]] = data[i]
eqtl_infor = eqtl_infor[np.where(np.logical_and(eqtl_infor['e_chr'] == gwas_chr, np.logical_and(eqtl_infor['e_end'] >= gstart - 1e6,
                                                                                                eqtl_infor['e_start'] <= gstop + 1e6)))[0]]
transcript_list = set([str(t) for t in eqtl_infor['trait']])

# Transcript level
fs = open(texpression_pheno_fname)
header = fs.readline().rstrip().split()
fs.close()
cols = [0] + [i for i in range(1, len(header)) if header[i] in transcript_list]
texpression_pheno = np.genfromtxt(texpression_pheno_fname,
                                  skip_header=1,
                                  usecols=cols,
                                  missing_values="NA",
                                  filling_values=np.nan,
                                  dtype=np.dtype([(header[0], 'U8')] + [(header[i].replace('.', '_'), np.float32) for i in cols[1:]]))
strain_set = strain_set.intersection(set([str(texpression_pheno['strain'][i]) for i in range(texpression_pheno.shape[0])]))

# Get the genotype at the peak marker
gm_valid_strains = np.array([i for i in range(gm_strains.shape[0]) if str(gm_strains[i]) in strain_set])
genotypes = genotype_matrix['genotype'][np.where(genotype_matrix['pos'] == gwas_peak)[0][0], gm_valid_strains]

ph_valid_strains = np.array([i for i in range(trait_phenotype.shape[0]) if str(trait_phenotype['strain'][i]) in strain_set])
phenotypes = trait_phenotype['value'][ph_valid_strains]

exp_valid_strains = np.array([i for i in range(texpression_pheno.shape[0]) if str(texpression_pheno['strain'][i]) in strain_set])
texpression_pheno = texpression_pheno[exp_valid_strains]

# List to store results for each transcript
multimed_trait_list = []
Y=[]
# Loop through unique transcripts
for trss in np.unique(eqtl_infor['trait']):
    trss2 = trss.replace('.', '_')
    valid = np.where(np.logical_and(genotypes != 0, np.logical_and(np.logical_not(np.isnan(texpression_pheno[trss2])),
                                                                   np.logical_not(np.isnan(phenotypes)))))[0]
    
    # Run medTest
    mt_multi_transcript = medTest(
        E=genotypes[valid],
        M=texpression_pheno[trss2][valid],
        Y=phenotypes[valid],
        nperm=1000
    )
    # print(trss, np.sum(genotypes[valid]), np.sum(texpression_pheno[trss2][valid]),np.sum(phenotypes[valid]))
    Y.append(f"{mt_multi_transcript[0, 1]:0.3f}")

    multimed_trait_list.append((trss, mt_multi_transcript[0, 0], mt_multi_transcript[0, 1]))

tmp=",".join(Y)
# print(f"Y=c({tmp})")
# Combine all results
multimed_traits = np.array(multimed_trait_list, dtype=np.dtype([('gene', f"U{max(len(row[0]) for row in multimed_trait_list)}"),
                                                                ('S', np.float32), ('p', np.float32)]))

# Calculate 99th percentile of S
q99_S = np.quantile(multimed_traits['S'], 0.99)

# Filter significant genes
sig_genes = np.where(np.logical_or(multimed_traits['p'] < 0.05,
                                   multimed_traits['S'] > q99_S))[0]

if sig_genes.shape[0] > 0:
    output = open(f"{out_prefix}_medmulti.tsv", 'w')
    tmp = "\t".join(['gene', 'S', 'p'] + list(eqtl_infor.dtype.names[1:]))
    output.write(f"{tmp}\n")
    for i in range(multimed_traits.shape[0]):
        where = np.where(eqtl_infor['trait'] == multimed_traits['gene'][i])[0][0]
        tmp = "\t".join([str(x) for x in (list(multimed_traits[i]) + list(eqtl_infor[where])[1:])])
        output.write(f"{tmp}\n")
    output.close()

    output = open(f"{out_prefix}_siggenes.tsv", 'w')
    tmp = "\n".join([str(multimed_traits['gene'][x]) for x in sig_genes])
    output.write(f"{tmp}\n")
    output.close()
else:
    print("No significant genes found.")