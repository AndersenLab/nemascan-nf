#!/usr/bin/env python

import sys
import argparse

import numpy

def main():
    # usage: aggregate_mappings.py [options] genotype_matrix traits -o out_prefix
    parser = argparse.ArgumentParser()
    parser.add_argument('gwas_mapping')
    parser.add_argument('-i', '--independent_tests', type=float)
    parser.add_argument('-s', '--snp_grouping', default=1000, type=int)
    parser.add_argument('-c', '--CI_size', default=150, type=int)
    parser.add_argument('-t', '--significance_threshold', default=None)
    parser.add_argument('-m', '--method', choices=['inbred', 'loco'])
    parser.add_argument('--trait')
    parser.add_argument('--chromosome_names')
    parser.add_argument('-o', '--output', required=True)

    args = parser.parse_args()

    # Load chromosome name mapping
    if args.chromosome_names is not None:
        int2chr = {int(line.rstrip().split("\t")[1]): line.split("\t")[0] for line in open(args.chromosome_names) if not line.startswith('chromosome')}
    else:
        int2chr = {x:x for x in range(numpy.unique(map_data['CHR']).shape[0])}
    
    # Load marker mapping data, using "method" to determine format
    if args.method == "inbred":
        map_data = numpy.loadtxt(args.gwas_mapping, skiprows=1, dtype=numpy.dtype([
            ('CHR', 'U5'), ('marker', 'U12'), ('POS', numpy.int32), ('A1', 'U1'), ('A2', 'U1'), ('N', numpy.int32),
            ('AF1', numpy.float32), ('BETA', numpy.float32), ('SE', numpy.float32),
            ('log10p', numpy.float32)]))
    else:
        map_data = numpy.loadtxt(args.gwas_mapping, skiprows=1, dtype=numpy.dtype([
            ('CHR', 'U5'), ('marker', 'U12'), ('POS', numpy.int32), ('A1', 'U1'), ('A2', 'U1'),
            ('AF1', numpy.float32), ('BETA', numpy.float32), ('SE', numpy.float32),
            ('log10p', numpy.float32)]))
    map_data['log10p'] = -numpy.log10(map_data['log10p'])
    map_data = map_data[numpy.where(map_data['log10p'] != 0)]
    map_data = map_data[numpy.lexsort((map_data['POS'], map_data['CHR']))]
    for i in range(map_data.shape[0]):
        tmp = map_data['marker'][i].split(':')
        map_data['marker'][i] = ":".join([int2chr[int(tmp[0])], tmp[1]])
        map_data['CHR'][i] = int2chr[int(map_data['CHR'][i])]

    # Determine the QTL cutoff based on the type of significance test
    significance_threshold = args.significance_threshold
    if significance_threshold == "EIGEN":
        QTL_cutoff = args.independent_tests
        BF = -numpy.log10(0.05 / QTL_cutoff)
    elif significance_threshold == "BF":
        QTL_cutoff = None
        BF = -numpy.log10(0.05 / numpy.where(numpy.logical_not(numpy.isnan(map_data['log10p'])))[0].shape[0])
    else:
        QTL_cutoff = float(significance_threshold)
        BF = QTL_cutoff
    aboveBF = numpy.zeros(map_data.shape[0], bool)
    where = numpy.where(numpy.logical_not(numpy.isnan(map_data['log10p'])))[0]
    aboveBF[where] = map_data['log10p'][where] >= BF
    
    snp_grouping = int(args.snp_grouping)
    CI_size = int(args.CI_size)

    # Make sure there is at least one significant SNP
    intervals = []
    if numpy.sum(aboveBF) > 0:
    
        # Make sure that there aren't too many significant SNPs
        sig_indices = numpy.where(aboveBF)[0]
        sigSnps = map_data[sig_indices]
        if sigSnps.shape[0] > map_data.shape[0] * 0.15:
            raise RuntimeError("Too many significant SNPs detected!")

        # Find neighborhoods around significant SNPs
        start = 0
        stop = 1
        chr_indices = numpy.r_[0, numpy.where(map_data['CHR'][1:] != map_data['CHR'][:-1])[0] + 1, map_data.shape[0]]
        chroms = list(map_data['CHR'][chr_indices[:-1]])
        group = 1
        while start < sigSnps.shape[0]:
            while stop < sigSnps.shape[0] and sigSnps['CHR'][start] == sigSnps['CHR'][stop] and sig_indices[stop] - sig_indices[start] < snp_grouping:
                stop += 1
            chrint = chroms.index(sigSnps['CHR'][start])
            s, e = chr_indices[chrint:chrint + 2]
            upstream_index = max(s, sig_indices[start] - CI_size)
            downstream_index = min(e - 1, sig_indices[stop - 1] + CI_size)
            best_pval = numpy.argmax(sigSnps['log10p'][start:stop]) + start
            intervals.append((sigSnps['CHR'][start], sigSnps['marker'][best_pval], sigSnps['log10p'][best_pval],
                              map_data['POS'][upstream_index], sigSnps['POS'][best_pval],
                              map_data['POS'][downstream_index], group))
            group += 1
            start = stop
            stop = start + 1

    # Write significant SNPs
    output = open(args.output, "w")
    tmp = "\t".join(["CHROM", "marker", "log10p", "startPOS", "peakPOS", "endPOS", "peak_id", "method", "trait"])
    output.write(f"{tmp}\n")
    if len(intervals) > 0:
        for i in range(len(intervals)):
            tmp = "\t".join([str(x) for x in intervals[i]])
            output.write(f"{tmp}\t{args.method}\t{args.trait}\n")
    output.close()
            

main()
