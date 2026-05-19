#!/usr/bin/env python

import sys
import gzip


def main():
    gwa_fname, annotation_fname, out_fname = sys.argv[1:4]

    gwa = load_gwa(gwa_fname)
    append_annotations(annotation_fname, gwa)
    header = "\t".join(gwa['header'])
    del gwa['header']
    markers = list(gwa.keys())
    markers.sort()
    output = open(out_fname, 'w')
    output.write(f"{header}\n")
    for marker in markers:
        tmp = "\t".join(gwa[marker])
        output.write(f"{tmp}\n")
    output.close()

def load_gwa(fname):
    data = {}
    fs = open(fname)
    header = fs.readline().upper().rstrip().split("\t")
    if 'bp' in header:
        pos_index = header.index("BP")
    else:
        pos_index = header.index("POS")
    chr_index = header.index("CHR")
    data['header'] = header
    for line in fs:
        line = line.rstrip().split('\t')
        data[(line[chr_index], int(line[pos_index]))] = line + [[]]
    return data

def append_annotations(fname, gwa):
    consequences = {
        'N/A'          : (0, 'Intergenic'),
        'intron'       : (1, 'Intronic'),
        'intronic'     : (1, 'Intronic'),
        'intron_variant': (1, 'Intronic'),
        'missense'     : (5, 'Missense'),
        '*missense'     : (5, 'Missense'),
        'missense_variant': (5, 'Missense'),
        'synonymous'   : (4, 'Synonymous'),
        '*synonymous'   : (4, 'Synonymous'),
        'synonymous_variant': (4, 'Synonymous'),
        '3_prime_utr'  : (2, "3'-UTR"),
        '3_prime_UTR_variant': (2, "3'-UTR"),
        'UTR3'         : (2, "3'-UTR"),
        '5_prime_utr'  : (2, "5'-UTR"),
        '5_prime_UTR_variant': (2, "5'-UTR"),
        'UTR5'         : (2, "5'-UTR"),
        'frameshift'   : (10, 'Frameshift'),
        'frameshift_variant': (10, 'Frameshift'),
        'inframe_deletion': (7, 'Nonframeshift deletion'),
        'non_coding'   : (3, 'Non-coding'),
        'splicing'     : (6, 'Splice region'),
        'splice_region': (6, 'Splice region'),
        'splice_region_variant': (6, 'Splice region'),
        'splice_donor' : (6, 'Splice region'),
        'splice_donor_variant': (6, 'Splice region'),
        'splice_donor_5th_base_variant': (6, 'Splice region'),
        'splice_donor_region_variant': (6, 'Splice region'),
        'splice_polypyrimidine_tract_variant': (6, 'Splice region'),
        'splice_acceptor': (6, 'Splice region'),
        'splice_acceptor_variant': (6, 'Splice region'),
        'start_lost'   : (12, 'Start lost'),
        'start_gained' : (8, 'Start gained'),
        'stop_retained': (4, 'Stop retained'),
        'stop_retained_variant': (4, 'Stop retained'),
        'stop_lost'    : (9, 'Stop lost'),
        'stop_gained'  : (11, 'Stop gained'),
        'upstream\\x3bdownstream': (0, 'Intergenic'),
        'upstream'     : (0, 'Intergenic'),
        'upstream_gene_variant': (0, 'Intergenic'),
        'downstream'   : (0, 'Intergenic'),
        'downstream_gene_variant': (0, 'Intergenic'),
        'ncRNA_splicing': (5, 'ncRNA Splicing'),
        'ncRNA_exonic' : (5, 'ncRNA Exonic'),
        'ncRNA_intronic' : (5, 'ncRNA Intronic'),
        'intergenic'   : (0, 'Intergenic'),
        'intergenic_variant': (0, 'Intergenic'),
        'intragenic_variant': (0, "Intronic"),
        'non_coding_transcript_variant': (3, 'Non-coding'),
        'non_coding_transcript_exon_variant': (3, 'Non-coding'),
        'exonic'       : {
            'synonymous_SNV'        : (4, 'Synonymous'),
            'nonsynonymous_SNV'     : (5, 'Missense'),
            'nonframeshift_deletion': (7, 'Nonframeshift deletion'),
            'frameshift_deletion'   : (10, 'Frameshift'),
            'stopgain'              : (11, 'Stop gained'),
            'stoploss'              : (9, 'Stop lost'),
            'startloss'             : (12, 'Start lost'),
            'startgain'             : (8, 'Start gained'),
        }

    }
    initial_length = len(gwa['header'])
    gwa['header'] += ['WBGENE', 'GENE_NAME', 'CONSEQUENCE', 'AA']
    with gzip.open(fname, 'rb') as fs:
        header = fs.readline().decode('utf8').rstrip().split(',')
        if header[5] == 'AA':
            format = 'csq'
            columns = ['WBGENE', 'GENE_NAME', 'CONSEQUENCE', 'AA']
            col_indices = [10, 11, 4, 5]
        elif header[8] == 'TRANSCRIPT':
            format = 'snpeff'
            columns = ['WBGENE', 'GENE_NAME', 'CONSEQUENCE', 'AA']
            col_indices = [9, 10, 4, 6]
        else:
            format = 'annovar_vep'
            columns = ['WBGENE', 'GENE_NAME', 'CONSEQUENCE', 'AA']
            impact_col = 5
            col_indices = [10, 11, 4, 6]
        for line in fs:
            line = line.decode('utf8').rstrip().split(',')
            marker = (line[0], int(line[1]))
            if marker not in gwa:
                continue
            WBGENE = line[col_indices[0]]
            GENE = line[col_indices[1]]
            AA = line[col_indices[3]]
            consequence = line[col_indices[2]].replace('\\x3b', '&').split('&')
            if len(consequence) == 1:
                consequence[0] = consequence[0].lstrip('*')
                if consequence[0] in consequences:
                    if consequence[0] == 'exonic':
                        consequence = consequences[consequence[0]][line[impact_col]]
                    else:
                        consequence = consequences[consequence[0]]
                    gwa[marker][-1].append((consequence, WBGENE, GENE, AA))
            else:
                for c in consequence:
                    c = c.lstrip('*')
                    if c in consequences:
                        if c == 'exonic':
                            c = consequences[c][line[impact_col]]
                        else:
                            c = consequences[c]
                        gwa[marker][-1].append((c, WBGENE, GENE, AA))
    for key in gwa:
        if key == 'header':
            continue
        if len(gwa[key][-1]) == 0:
            gwa[key] = gwa[key][:-1] + ['N/A', 'N/A', 'N/A', 'N/A']
        elif len(gwa[key][-1]) == 1:
            consequence, wbgene, gene, aa = gwa[key][-1][0]
            gwa[key] = gwa[key][:-1] + [wbgene, gene, consequence[1], aa]
        else:
            gwa[key][-1].sort()
            consequence, wbgene, gene, aa = gwa[key][-1][-1]
            gwa[key] = gwa[key][:-1] + [wbgene, gene, consequence[1], aa]
    return

main()
