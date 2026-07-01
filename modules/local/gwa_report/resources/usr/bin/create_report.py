#!/usr/bin/env python 
"""
Generic HTML Report Generator
Creates professional reports with text, charts, and table of contents
"""

import argparse
import os
from datetime import datetime
from typing import List, Dict, Any, Optional
from dataclasses import dataclass, field
import math
import re
import json
import gzip

import numpy as np


def main():
    parser = make_arg_parser()
    args = parser.parse_args()
    data = parse_config(args.config, args.trait)
    template = open(args.template).read()
    template = template.replace("ALLDATA", render_large_json(data))
    # print_data(data)

    output = open(args.output, 'w')
    output.write(template)
    output.close()


def print_data(data, indent=0):
    if indent == 0:
        print("alldata:")
    if type(data) == dict:
        for key in data.keys():
            print(f"{'  ' * indent}- {key}:", end="")
            if type(data[key]) in [list, tuple]:
                if len(data[key]) == 0:
                    continue
                if type(data[key][0]) in [list, tuple]:
                    print('')
                    for line in data[key][:min(5, len(data[key]))]:
                        tmp = ", ".join([str(x) for x in line[:min(5, len(line))]])
                        if len(line) > 5:
                            print(f"{'  ' * (indent + 1)}- [{tmp}, ...]")
                        else:
                            print(f"{'  ' * (indent + 1)}- [{tmp}]")
                    if len(data[key]) > 5:
                        print(f"{'  ' * (indent + 1)}...")
                else:
                    tmp = ", ".join([str(x) for x in data[key][:min(5, len(data[key]))]])
                    if len(data[key]) > 5:
                        print(f" [{tmp}, ...]")
                    else:
                        print(f" [{tmp}]")
            elif type(data[key]) in [str, int, float]:
                print(f" {data[key]}")
            else:
                print("")
                print_data(data[key], indent+1)


def render_large_json(data):
    print(data['highlight_strains'])
    text = toJson(data)
    # text = json.dumps(data)#.replace(' ', '')
    split_text = []
    start = 0
    end = start + 1
    size = 0
    while end < len(text):
        if size + len(text[end]) > 2000:
            split_text.append(''.join(text[start:end]))
            start = end
            end = start + 1
            size = len(text[start])
        else:
            size += len(text[end])
            end += 1
    split_text.append(''.join(text[start:end]))
    return "\n".join(split_text)


def toJson(data):
    text = []
    if type(data) == dict:
        keys = list(data.keys())
        text.append("{")
        for i, key in enumerate(keys):
            if i == 0:
                text.append(f'"{key}":')
            else:
                text.append(f',"{key}":')
            text += toJson(data[key])
        text.append("}")
    elif type(data) in [list, tuple]:
        text.append('[')
        for i, entry in enumerate(data):
            if i > 0:
                text += ','
            text += toJson(entry)
        text.append(']')
    elif type(data) == int:
        text.append(f"{data}")
    elif type(data) == float:
        if np.isnan(data):
            text.append("NaN")
        else:
            text.append(f"{data:0.3e}")
    elif type(data) == str:
        text.append('"' + data + '"')
    return text

    
def make_arg_parser():
    parser = argparse.ArgumentParser(
        description='Create an interactive Manhattan plot from GWAS results'
    )
     
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output report filename'
    )

    parser.add_argument(
        '--template',
        required=True,
    )
    
    parser.add_argument(
        '-c', '--config',
        required=True,
        help='Cofiguration file'
    )
    
    parser.add_argument(
        '-t', '--trait',
        required=True,
        help='Trait name'
    )
    
    return parser


def parse_config(fname, trait):
    data = {'broad': {'maxY': 0}, 'params': {}, 'stats': {}, 'methods': [], 'highlight_strains': ["none"]}
    config = {}
    with open(fname, 'r') as f:
        for line in f:
            line = line.rstrip().split("\t")
            if len(line) < 2:
                continue
            key = line[0]
            if key == 'trait':
                strains, phenotype = load_phenotypes(line[1])
                data['broad']['phenotype'] = phenotype
                data['broad']['strain'] = strains
            elif key == 'broad_gwa':
                method, gwa_fname = line[1:3]
                if method not in data['methods']:
                    data['methods'].append(method)
                data['broad'].setdefault(method, {})
                marker, pval, logpval = load_gwa(gwa_fname, fine=False)
                if 'marker' not in data['broad']: 
                    data['broad']['marker'] = marker
                data['broad'][method]['pval'] = pval
                data['broad'][method]['logpval'] = logpval
                for chrom in logpval.keys():
                    data['broad']['maxY'] = max(data['broad']['maxY'], max([float(x) for x in logpval[chrom]]))
            elif key == 'broad_qtl':
                method, qtl_fname = line[1:3]
                data['broad'].setdefault(method, {})
                data['broad'][method].setdefault('QTL', {})
                with open(qtl_fname) as fs:
                    _ = fs.readline()
                    for line in fs:
                        chrom, marker, logpval, start, peak, end = line.rstrip().split("\t")[:6]
                        data['broad'][method]['QTL'][marker] = [chrom, int(start), int(end), int(peak), float(logpval)]
            elif key == 'finemap_gwa':
                method, chrom, start, end, peak, gwa_fname = line[1:7]
                name = f"{chrom}:{peak}"
                marker, pval, logpval, variants, LD = load_gwa(gwa_fname, fine=True)
                data.setdefault('fine', {})
                data['fine'].setdefault(method, {})
                data['fine'][method].setdefault(name, {})
                if name not in data['fine']:
                    data['fine'][name] = {}
                    data['fine'][name]['marker'] = marker
                    data['fine'][name]['LD'] = LD
                    if len(variants) > 0:
                        data['fine'][name]['variants'] = variants
                else:
                    if marker[0] < data['fine'][name]['marker'][0]:
                        index = np.searchsorted(np.array(marker), data['fine'][name]['marker'][0])
                        data['fine'][name]['marker'] = marker[:index] + data['fine'][name]['marker']
                        data['fine'][name]['LD'] = LD[:index] + data['fine'][name]['LD']
                        if len(variants) > 0:
                            data['fine'][name]['variants'] = variants[:index] + data['fine'][name]['variants']
                    if marker[-1] > data['fine'][name]['marker'][-1]:
                        index = np.searchsorted(np.array(marker), data['fine'][name]['marker'][-1], side='right')
                        data['fine'][name]['marker'] += marker[index:]
                        data['fine'][name]['LD'] += LD[index:]
                        if len(variants) > 0:
                            data['fine'][name]['variants'] += variants[index:]
                data['fine'][method][name]['logpval'] = logpval
            elif key == 'independent_tests':
                data['stats'][key] = float(line[1])
            elif key == 'narrow_h2':
                data['stats'][key] = float(open(line[1]).readline().rstrip().split("\t")[-1])
            elif key == 'mediation':
                method, chrom, start, end, peak, mediation_fname = line[1:7]
                name = f"{chrom}:{start}-{end}"
                data.setdefault('mediation', {})
                if name not in data['mediation']:
                    data['mediation'][name] = load_mediation(mediation_fname) 
            elif key == 'issues':
                data[key] = [line.rstrip() for line in open(line[1]).readlines()]
            elif key == 'genotype_matrix':
                config[key] = line[1]
            elif key == "haplotypes":
                config[key] = line[1]
            elif key == 'highlight_strains':
                data[key] = line[1].split(',')
            elif key == 'params':
                for pair in line[1:]:
                    k, v = pair.split("=")
                    data[key][k] = v
            elif key == 'genes':
                config[key] = line[1]
                
    data['params']['trait'] = trait
    data['methods'].sort()

    # need to wait to load these data until all markers are identified
    all_markers = {}
    for method in data['methods']:
        for marker, marker_data in data['broad'][method]['QTL'].items():
            if marker not in all_markers:
                all_markers[marker] = marker_data[0:3]
            else:
                all_markers[marker][1] = min(all_markers[marker][1], marker_data[1])
                all_markers[marker][2] = max(all_markers[marker][2], marker_data[2])
    strains = data['broad']['strain']
    marker_genotypes = load_genotype_matrix(config['genotype_matrix'], all_markers, strains)
    data['broad']['qtl_genotype'] = {}
    for key, value in marker_genotypes.items():
        data['broad']['qtl_genotype'][key] = value # {ref, alt, [genotypes]}
    if 'fine' in data and 'haplotypes' in config:
        marker_haplotypes = load_haplotypes(config['haplotypes'], all_markers, strains)
        for key, value in marker_haplotypes.items():
            data['fine'][key]['haplotype'] = value # [[(start, end, name, color), ...], [(start, end, name, color), ...], ...]
    if 'fine' in data and 'genes' in config and len(all_markers) > 0:
        data['genes'] = load_genes(config['genes'], all_markers)
    
    # Calculate position offsets for plotting all chroms
    chroms = [(chrom.replace('MtDNA', 'ZMtDNA'), str(chrom)) for chrom in data['broad']['marker'].keys()]
    chroms.sort()
    chroms = [chrom[1] for chrom in chroms]
    starts = [int(data['broad']['marker'][chrom][0]) for chrom in chroms]
    ends = [int(data['broad']['marker'][chrom][-1]) for chrom in chroms]
    spacer = int(np.round(sum([ends[i] for i in range(len(chroms))]) * 0.1 / (len(chroms) - 1)))
    offsets = [0] + [sum(ends[:i+1]) + spacer * (i + 1) for i in range(len(chroms) - 1)]
    data['chrom_regions'] = {'chroms': chroms, 'starts': starts, 'ends': ends, 'offsets': offsets}
    
    # Find the significance cutoff
    n = sum([len(data['broad']['marker'][chrom]) for chrom in data['chrom_regions']['chroms']])
    if data['params']['significance_threshold'] == 'BF':
        data['params']['significance_cutoff'] = [('BF', float(data['params']['alpha']) / n, 'dash', '#D41159')]
    elif data['params']['significance_threshold'] == 'EIGEN':
        data['params']['significance_cutoff'] = [('BF', float(data['params']['alpha']) / n, 'dash', '#D41159'),
                                                 ('EIGEN', float(data['params']['alpha']) / data['stats']['independent_tests'], 'dotted', '#D35FB7')]
    else:
        data['params']['significance_cutoff'] = [('USER', float(data['params']['significance_threshold']), 'dash', '#D41159')]

    # Calculate QTL LD if multiple QTL are present
    for method in data['methods']:
        if len(data['broad'][method]['QTL']) > 1:
            data['broad'][method]['LD_matrix'] = find_QTL_LD(data, method)
    return data


def load_phenotypes(filename):
    """
    Read phenotype data from a file
    
    Parameters:
    -----------
    filename : str
        Path to the phenotype file
    
    Returns:
    --------
    dict mapping strain names to phenotype values
    """
    strains = []
    values = []
    with open(filename, 'r') as f:
        _ = f.readline()
        for line in f:
            line = line.rstrip().split("\t")
            strains.append(line[0])
            if line[1] == 'NA':
                values.append(float(np.nan))
            else:
                values.append(float(line[1]))
    return strains, values


def load_gwa(filename, fine=False):
    """
    Read GWAS summary statistics from a file
    
    Parameters:
    -----------
    filename : str
        Path to the GWAS results file
    pos_col : str
        Column name for position
    
    Returns:
    --------
    list of dicts containing GWAS data
    """
    markers = {}
    pvals = {}
    logpvals = {}
    if fine:
        variants = {}
        LDs = {}
    with open(filename, 'r') as f:
        header = [x.upper() for x in f.readline().rstrip().split("\t")]
        chrom_col = header.index("CHR")
        if "POS" in header:
            pos_col = header.index("POS")
        else:
            pos_col = header.index("BP")
        pval_col = header.index('P')
        if fine:
            ld_col = header.index('LD')
            if "WBGENE" in header:
                annotations = True
                wbgene_col = header.index('WBGENE')
                gene_col = header.index('GENE_NAME')
                consequence_col = header.index('CONSEQUENCE')
                aa_col = header.index('AA')
            else:
                annotations = False

        for line in f:
            row = line.rstrip().split("\t")
            pval = float(row[pval_col])
            if pval <= 0 or pval > 1:
                log_pval = float(np.nan)
            else:
                log_pval = float(-np.log10(pval))

            chrom = row[chrom_col].replace('chr', '').replace('Chr', '')
            pos = int(row[pos_col])
            markers.setdefault(chrom, [])
            pvals.setdefault(chrom, [])
            logpvals.setdefault(chrom, [])
            markers[chrom].append(pos)
            pvals[chrom].append(pval)
            logpvals[chrom].append(log_pval)
            if fine:
                ld = float(row[ld_col])
                LDs.setdefault(chrom, [])
                LDs[chrom].append(ld)
                if annotations:
                    wbgene = row[wbgene_col]            
                    gene = row[gene_col]            
                    consequence = row[consequence_col]
                    aa = row[aa_col]
                    ref = row[pos_col + 1]
                    alt = row[pos_col + 2]
                    variants.setdefault(chrom, [])
                    variants[chrom].append((wbgene, gene, consequence, aa, ref, alt))
    if fine:
        chrom = list(markers.keys())[0]
        markers = markers[chrom]
        pvals = pvals[chrom]
        logpvals = logpvals[chrom]
        LDs = LDs[chrom]
        if annotations:
            variants = variants[chrom]
        return markers, pvals, logpvals, variants, LDs
    return markers, pvals, logpvals


def load_mediation(filename):
    data = []
    for line in open(filename):
        line = line.rstrip().split("\t")
        for i in range(len(line)):
            try:
                line[i] = int(line[i])
                continue
            except:
                pass
            try:
                line[i] = float(line[i])
                continue
            except:
                pass
        data.append(line)
    return data


def load_genotype_matrix(filename, markers, strains):
    """
    Load genotype matrix from a file
    
    Parameters:
    -----------
    filename : str
        Path to the genotype matrix file
    
    Returns:
    --------
    dict mapping markers to genotype arrays
    """
    genotype_matrix = {}
    with open(filename, 'r') as f:
        header = f.readline().rstrip().split("\t")
        cols = [header.index(strain) for strain in strains]
        for line in f:
            line = line.rstrip().split("\t")
            chrom, pos, ref, alt = line[0:4]
            marker = f"{chrom}:{pos}"
            if marker not in markers:
                continue
            genotype_matrix[marker] = {'ref':ref, 'alt': alt, 'genotype': tuple([int(line[x]) for x in cols])}
    return genotype_matrix


def load_haplotypes(filename, markers, strains):
    """
    Load haplotype data from a file
    
    Parameters:
    -----------
    filename : str
        Path to the haplotype file
    
    Returns:
    --------
    list of strain HDRs, each entry being a list genomic intervals
    """
    all_haps = {}
    with open(filename, 'r') as f:
        for line in f:
            line = line.rstrip().split("\t")
            chrom, start, end, strain, name, _, _, color = line
            all_haps.setdefault(chrom, {})
            all_haps[chrom].setdefault(strain, [])
            all_haps[chrom][strain].append((int(start), int(end), name, color))
    for chrom in all_haps.keys():
        for strain in all_haps[chrom].keys():
            all_haps[chrom][strain].sort()
            all_haps[chrom][strain] = np.array(all_haps[chrom][strain], dtype=np.dtype([('start', int), ('end', int), ('name', 'U8'), ('color', 'U7')]))

    haps = {}
    for marker, (chrom, start, stop) in markers.items():
        haps[marker] = []
        for strain in strains:
            start_i = np.searchsorted(all_haps[chrom][strain]['end'], start)
            stop_i = np.searchsorted(all_haps[chrom][strain]['start'], end)
            haps[marker].append([(int(all_haps[chrom][strain]['start'][i]),
                                  int(all_haps[chrom][strain]['end'][i]),
                                  str(all_haps[chrom][strain]['name'][i]),
                                  str(all_haps[chrom][strain]['color'][i])) for i in range(start_i, stop_i)])
    return haps


def find_QTL_LD(data, method):
    LD = []
    markers = [(marker.replace('MtDNA', 'ZMtDNA'), str(marker)) for marker in data['broad'][method]['QTL'].keys()]
    markers.sort()
    markers = [marker[1] for marker in markers]
    for i, marker1 in enumerate(markers):
        LD.append([])
        for j, marker2 in enumerate(markers):
            LD[-1].append(1)
    for i, marker1 in enumerate(markers):
        for j, marker2 in enumerate(markers):
            if j <= i:
                continue
            LD[i][j] = float(np.corrcoef(np.array(data['broad']['qtl_genotype'][marker1]['genotype']),
                                         np.array(data['broad']['qtl_genotype'][marker2]['genotype']))[0, 1]) ** 2
            LD[j][i] = LD[i][j]
    return {'marker': markers, 'LD': LD}


def load_genes(fname, markers):
    data = {}
    with gzip.open(fname, 'rb') as fs:
        for line in fs:
            line = line.decode('utf8').rstrip().split("\t")
            chrom = line[0]
            gene = line[3].split(':')[-1].split('.')
            while not gene[1].isnumeric():
                gene[1] = gene[1][:-1]
            gene = f"{gene[0]}.{gene[1]}"
            esizes = [int(x) for x in line[10].rstrip(',').split(',')]
            estarts = [int(x) + int(line[1]) for x in line[11].rstrip(',').split(',')]
            data.setdefault(chrom, {})
            if gene not in data[chrom]:
                data[chrom][gene] = [int(line[1]), int(line[2]), line[5], int(line[6]), int(line[7]), set()]
            else:
                data[chrom][gene][0] = min(data[chrom][gene][0], int(line[1]))
                data[chrom][gene][1] = max(data[chrom][gene][1], int(line[2]))
                data[chrom][gene][3] = min(data[chrom][gene][3], int(line[6]))
                data[chrom][gene][4] = max(data[chrom][gene][4], int(line[7]))
            for i in range(int(line[9])):
                data[chrom][gene][5].add((estarts[i], estarts[i] + esizes[i]))
    for chrom in data.keys():
        genes = []
        for gene in data[chrom].keys():
            exons = [list(x) for x in data[chrom][gene][-1]]
            exons.sort()
            for i in range(1, len(exons))[::-1]:
                if exons[i - 1][1] >= exons[i][0]:
                    exons[i - 1][1] = exons[i][1]
                    del exons[i]
            genes.append(data[chrom][gene][:5] + [exons, gene])
        genes.sort()
        data[chrom] = genes
        # start, stop, strand, thickstart, thickstop, exons, name
    genes = {}
    layers = []
    for marker, (chrom, start, stop) in markers.items():
        genes[marker] = {'bounds': [0, 0], 'genes': []}
        for (gstart, gstop, strand, thickstart, thickstop, exons, name) in data[chrom]:
            curr_level = 0
            if (gstop < start or gstart > stop):
                continue
            while curr_level < len(layers) and gstart < layers[curr_level]:
                curr_level += 1
            if curr_level == len(layers):
                layers.append(gstop + 1000)
            else:
                layers[curr_level] = gstop + 1000
            Y0 = curr_level
            s = gstart / 1e6
            e = gstop / 1e6
            offset = 0.2
            if thickstart == thickstop:
                s1 = s
                e1 = e
                offset = 0.05
            else:
                s1 = thickstart / 1e6
                e1 = thickstop / 1e6
            X = [s]
            Y = [offset]
            exons = np.array(exons) / 1e6
            for i in range(exons.shape[0]):
                if s1 > s and s1 < exons[i, 0] and (i == 0 or s1 > exons[i - 1, 1]):
                    X += [s1, s1]
                    Y += [offset, 0.05]
                    offset = 0.05
                X += [exons[i, 0], exons[i, 0]]
                Y += [offset, 0.4]
                if s1 < exons[i, 1] and e1 > exons[i, 1]:
                    offset = 0.05
                else:
                    offset = 0.2
                X += [exons[i, 1], exons[i, 1]]
                Y += [0.4, offset]
                if e1 > exons[i, 1] and (i == exons.shape[0] - 1 or e1 < exons[i + 1, 0]):
                    X += [e1, e1]
                    Y += [offset, 0.2]
                    offset = 0.2
            X += [e]
            Y += [offset]
            X = np.array(X)
            Y = np.array(Y)
            X = np.r_[X, X[::-1], X[0]]
            Y = -(np.r_[Y, -Y[::-1], Y[0]] + Y0)
            genes[marker]['bounds'][0] = float(min(genes[marker]['bounds'][0], np.amin(Y)))
            genes[marker]['bounds'][1] = float(max(genes[marker]['bounds'][1], np.amax(Y)))
            genes[marker]['genes'].append([strand, name, [float(x) for x in X], [float(x) for x in Y]])
    return genes


if __name__ == "__main__":
    main()
