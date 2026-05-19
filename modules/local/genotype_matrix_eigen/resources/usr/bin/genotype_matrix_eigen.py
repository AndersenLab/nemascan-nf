#!/usr/bin/env python

import sys

import numpy

# Based on https://www.nature.com/articles/6800717

def main():
    header = open(sys.argv[1]).readline().rstrip().split()
    data = numpy.loadtxt(sys.argv[1], skiprows=1, usecols=[x for x in range(4, len(header))], dtype=numpy.int8)
    # Remove rows with too low or too high minor allele frequency
    MAF = numpy.sum(data > 0, axis=1) / numpy.sum(numpy.abs(data), axis=1)
    data = data[numpy.where(numpy.logical_and(MAF >= 0.05, MAF <= 0.95))[0]]
    significant_values = 0
    corr_matrix = numpy.corrcoef(data)
    del data
    eigenvalues = numpy.linalg.eigvals(corr_matrix).astype(numpy.float64)
    significant_values = numpy.sum(eigenvalues >= 1) + numpy.sum(eigenvalues - numpy.floor(eigenvalues))
    print(f"eigen\n{significant_values}")

def load_data(fname, chrom):
    with open(fname) as fs:
        header = fs.readline().strip().split("\t")
        data = []
        chroms = []
        for line in fs:
            line = line.strip().split("\t")
            if chrom is None:
                chroms.append(line[0])
                data.append(line[4:])
            elif chrom == line[0]:
                data.append(line[4:])
    data = numpy.array(data, numpy.int8)
    chroms = numpy.array(chroms)
    return data, chroms

main()
