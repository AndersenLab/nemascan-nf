#!/usr/bin/env python3
"""
Interactive Manhattan Plot Generator using Plotly
Creates an HTML file with an interactive Manhattan plot for GWAS results
"""

import argparse
import gzip
import math

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np


def main():
    parser = make_arg_parser()
    args = parser.parse_args()
    
    # Read GWAS data
    print(f"Reading GWAS data from {args.input}...")
    data = read_gwas_data(
        args.input,
        pos_col=args.pos_col,
        pval_col=args.pval_col,
        delimiter=args.delimiter
    )
    chrom = args.marker.split(':')[0]
    positions = {chrom: set([x for x in data['pos']])}
    
    print(f"Loaded {len(data)} SNPs")
    
    # Read variant annotation data
    annotations = load_annotations(args.annotation, positions=positions)

    # Read gene data
    genes = load_genes(args.genes)

    # Create Manhattan plot
    print("Creating Interval variant plot...")
    fig = create_interval_variant_plot(
        data,
        annotations[np.where(annotations['chr'] == chrom)],
        genes,
        args.marker,
        title=args.title,
        point_size=args.point_size,
    )
    
    # Save to HTML
    fig.write_html(
        args.output,
        include_plotlyjs=True,
        full_html=True,
        config={
            'displayModeBar': True,
            'modeBarButtonsToAdd': ['drawline', 'drawopenpath', 'eraseshape'],
            'toImageButtonOptions': {
                'format': 'png',
                'filename': 'interval_variant_plot',
                'height': args.dimensions[1],
                'width': args.dimensions[0],
                'scale': 2
            }
        }
    )
    
    print(f"Interval variant plot saved to: {args.output}")
    
def read_gwas_data(filename, chrom_names=None, pos_col='POS', pval_col='P', delimiter='\t'):
    """
    Read GWAS summary statistics from a file
    
    Parameters:
    -----------
    filename : str
        Path to the GWAS results file
    pos_col : str
        Column name for position
    pval_col : str
        Column name for p-value
    delimiter : str
        File delimiter
    
    Returns:
    --------
    list of dicts containing GWAS data
    """
    data = []
    with open(filename, 'r') as f:
        header = [x.upper() for x in f.readline().rstrip().split(delimiter)]
        pos_col = header.index(pos_col)
        pval_col = header.index(pval_col)
        
        for line in f:
            row = line.rstrip().split("\t")
            
            pos = int(row[pos_col])
            pval = float(row[pval_col])
            
            # Skip invalid p-values
            if pval <= 0 or pval > 1:
                continue
            
            data.append((pos, pval, -math.log10(pval)))
                
    dtype = np.dtype([('pos', np.int32), ('pval', np.float32), ('log_pval', np.float32)])
    data = np.array(data, dtype=dtype)
    data = data[np.argsort(data['pos'])]
    return data

def load_annotations(fname, positions):
    data = {}
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
            if line[0] not in positions:
                continue
            line[1] = int(line[1])
            if line[1] not in positions[line[0]]:
                continue
            data.setdefault(line[0], {})
            data[line[0]].setdefault(line[1], [])
            WBGENE = line[col_indices[0]]
            if WBGENE == 'N/A':
                WBGENE = ""
            GENE = line[col_indices[1]]
            if GENE == 'N/A':
                GENE = ""
            AA = line[col_indices[2]]
            if AA == 'N/A':
                AA = ""
            consequence = line[col_indices[2]].replace('\\x3b', '&').split('&')
            if len(consequence) == 1:
                consequence[0] = consequence[0].lstrip('*')
                if consequence[0] in consequences:
                    if consequence[0] == 'exonic':
                        consequence = consequences[consequence[0]][line[impact_col]]
                    else:
                        consequence = consequences[consequence[0]]
                    data[line[0]][line[1]].append((consequence, WBGENE, GENE, AA))
                elif not consequence[0].startswith("@"):
                    print(consequence[0])
            else:
                for c in consequence:
                    c = c.lstrip('*')
                    if c in consequences:
                        if c == 'exonic':
                            c = consequences[c][line[impact_col]]
                        else:
                            c = consequences[c]
                        data[line[0]][line[1]].append((c, WBGENE, GENE, AA))
                    elif not c.startswith("@"):
                        print(c)
    new_data = []
    str_len = [0, 0, 0, 0, 0]
    for chrom in data.keys():
        str_len[0] = max(str_len[0], len(chrom))
        for pos, entries in data[chrom].items():
            if len(entries) == 0:
                continue
            if len(entries) > 1:
                entries.sort()
            new_data.append(tuple([chrom, pos, entries[-1][1], entries[-1][2], entries[-1][0][1], entries[-1][3]]))
            str_len[1] = max(str_len[1], len(new_data[-1][2]))
            str_len[2] = max(str_len[2], len(new_data[-1][3]))
            str_len[3] = max(str_len[3], len(new_data[-1][4]))
    data = np.array(new_data, dtype=np.dtype([('chr', f'U{str_len[0]}'), ('pos', np.int32)] + [(name, f'U{str_len[i + 1]}') for i, name in enumerate(columns)]))
    data = data[np.lexsort((data['pos'], data['chr']))]
    return data

def load_genes(fname):
    data = {}
    with gzip.open(fname, 'rb') as fs:
        for line in fs:
            line = line.decode('utf8').rstrip().split("\t")
            chrom = line[0]
            gene = line[3].split(':')[0].split('.')
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
    return data

def create_interval_variant_plot(data,
                                 variants,
                                 genes,
                                 marker,
                                 title="QTL Variant Plot", 
                                 point_size=4,
                                ):
    """
    Create an interactive Manhattan plot using Plotly
    
    Parameters:
    -----------
    data : list
        List of dictionaries containing GWAS data
    variants : dict
        Dictionary of variants keyed by chromosome
    genes : dict
        Dictionary of genes keyed by chromosome
    marker : str
        Anchor marker
    title : str
        Plot title
    output_file : str
        Output HTML filename
    point_size : int
        Size of the data points
    """
    
    marker_index = np.searchsorted(data['pos'], int(marker.split(':')[-1]))
    marker_pos = int(marker.split(':')[-1]) / 1e6
    start = data['pos'][0]
    end = data['pos'][-1]
    chrom = marker.split(':')[0]
    colors = {
        'Intergenic': "#777777",
        'Genic Non-coding': '#ffb00c',
        'Coding-change': '#cc2181',
        'Synonymous': '#1d8838',
        'RNA': '#4477aa'
    }
    groups = {
        'Intergenic': 'Intergenic',
        'Non-coding': 'Genic Non-coding',
        'Intronic': 'Genic Non-coding',
        "3'-UTR": 'Genic Non-coding',
        "5'-UTR": 'Genic Non-coding',
        'Splice region': 'Genic Non-coding',
        'Missense': 'Coding-change',
        'Synonymous': 'Synonymous',
        'Frameshift': 'Coding-change',
        'Nonframeshift deletion': 'Coding-change',
        'Start lost': 'Coding-change',
        'Start gained': 'Coding-change',
        'Stop retained': 'Coding-change',
        'Stop list': 'Coding-change',
        'Stop gained': 'Coding-change',
        'ncRNA Splicing': 'RNA',
        'ncRNA Exonic': 'RNA',
        'ncRNA Intronic': 'RNA',
    }

    startX = np.amin(variants['pos']) / 1e6
    stopX = np.amax(variants['pos']) / 1e6
    spanX = stopX - startX
    startX -= 0.005 * spanX
    stopX += 0.005 * spanX
    startY = np.inf
    stopY = -np.inf
    grouped_variants = {}
    for name, group in groups.items():
        where = np.where(variants['CONSEQUENCE'] == name)[0]
        grouped_variants.setdefault(group, [])
        grouped_variants[group].append(variants[where])
    grouped_Y = {}
    for name in grouped_variants.keys():
        grouped_variants[name] = np.concatenate(grouped_variants[name], axis=0)
        grouped_variants[name] = grouped_variants[name][np.argsort(grouped_variants[name]['pos'])]
        indices = np.searchsorted(data['pos'], grouped_variants[name]['pos'])
        grouped_Y[name] = data['log_pval'][indices]
        if indices.shape[0] > 0:
            startY = min(startY, np.amin(data['log_pval'][indices]))
            stopY = max(stopY, np.amax(data['log_pval'][indices]))
    spanY = stopY - startY
    startY -= 0.005 * spanY
    stopY += 0.005 * spanY

    # Create figure
    if genes is None:
        rows = 1
        row_heights = [1]
    else:
        rows = 2
        row_heights = [6, 1]
    fig = make_subplots(rows=rows, cols=1, row_heights=row_heights, shared_xaxes=True, vertical_spacing=0)
    fig.add_trace(go.Scatter(
        x=[marker_pos, marker_pos],
        y=[startY, stopY],
        mode="lines",
        line=dict(
            color='black',
            dash='dash'
        ),
        showlegend=False,
    ), row=1, col=1)
    fig.add_trace(go.Scatter(
        x=[startX, stopX],
        y=[stopY, stopY-0.001 * spanY],
        mode="lines",
        line=dict(
            color='black',
            width=2,
        ),
        showlegend=False,
    ), row=1, col=1)

    if genes is not None:
        fig.add_trace(go.Scatter(
            x=[marker_pos, marker_pos],
            y=[startY, stopY],
            mode="lines",
            line=dict(
                color='black',
                dash='dash'
            ),
            showlegend=False,
        ), row=2, col=1)

        level = [[-1, -1]]
        gene_colors = {'+': '#000000', '-': '#666666'}
        gene_names = {'+': 'Positive Strand Gene', '-': 'Negative Strand Gene'}
        shown_genes = {'+': 0, '-': 0}
        for gene in genes[chrom]:
            if gene[1] < start or gene[0] > end:
                continue
            current_level = 0
            while gene[0] < level[current_level][1] + 2e3:
                current_level += 1
                if len(level) <= current_level:
                    level.append([-1, -1])
            level[current_level] = [gene[0], gene[1]]
            Y0 = current_level
            s = gene[0] / 1e6
            e = gene[1] / 1e6
            offset = 0.2
            if gene[3] == gene[4]:
                s1 = s
                e1 = e
                offset = 0.05
            else:
                s1 = gene[3] / 1e6
                e1 = gene[4] / 1e6
            X = [s]
            Y = [offset]
            exons = np.array(gene[5]) / 1e6
            utr_index = np.searchsorted(exons[:, 0], s1)
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
            Y = np.r_[Y, -Y[::-1], Y[0]] + Y0
            fig.add_trace(go.Scatter(
                x=X,
                y=Y,
                mode='lines',
                fill='toself',
                fillcolor=gene_colors[gene[2]],
                line=dict(color=gene_colors[gene[2]], width=0),
                showlegend=shown_genes[gene[2]] == 0,
                name = gene_names[gene[2]],
                text=f"{gene[6]} {gene[2]}",
                hoverinfo='text',
                hoveron='fills',
                ), row=2, col=1)
            shown_genes[gene[2]] += 1
        fig.update_yaxes(
            fixedrange=True,
            row=2,
            col=1,
            linewidth=2,
            linecolor='black',
            mirror=True,
            tickvals=[],
            ticktext=[],
            title='Genes')
        
    fig.update_xaxes(
        row=rows,
        col=1,
        linewidth=2,
        linecolor='black',
        ticks='outside',
        title='Genomic Position (Mb)',
        range=[startX, stopX]
    )

    for group in grouped_variants.keys():
        labels = []
        V = grouped_variants[group]
        for i in range(V.shape[0]):
            if V['WBGENE'][i] != "":
                gene = f" {V['GENE_NAME'][i]} ({V['WBGENE'][i]})"
            elif V['GENE_NAME'][i] != "":
                gene = f" {V['GENE_NAME'][i]}"
            else:
                gene = ""
            if V['AA'][i] != "":
                aa = f" {V['AA'][i]}"
            else:
                aa = ""
            labels.append(f"{chrom}:{V['pos'][i]} {V['CONSEQUENCE'][i]}{aa}{gene}")
        fig.add_trace(go.Scatter(
            x=V['pos'] / 1e6,
            y=grouped_Y[group],
            name=group,
            mode='markers',
            marker=dict(
                color=colors[group],
                size=point_size,
                symbol=23,
            ),
            text=labels,
            hoverinfo='text',
            showlegend=True,
        ))
    fig.update_yaxes(
        fixedrange=True,
        range=[startY, stopY],
        row=1,
        col=1,
        linewidth=2,
        linecolor='black',
        mirror=True,
        ticks='outside',
        title='-log<sub>10</sub>(P-value)')

    if title is not None:
        fig.update_layout(
            title=dict(
                text=title,
                font=dict(size=20),
                x=0.5
            ))
        upper_margin = 40
    else:
        upper_margin = 30
    fig.update_layout(
        xaxis=dict(
            range=[startX, stopX],
        ),
        margin=dict(l=30, r=0, t=upper_margin, b=0)
    )
    return fig

def make_arg_parser():
    parser = argparse.ArgumentParser(
        description='Create an interactive QTL variant plot from GWAS results'
    )
    
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input GWAS summary statistics file'
    )
    
    parser.add_argument(
        '-a', '--annotation',
        required=True,
        help='Input variant annotation file'
    )
    
    parser.add_argument(
        '-g', '--genes',
        required=True,
        help='Input gene file'
    )
    
    parser.add_argument(
        '-m', '--marker',
        required=True,
        help='Input anchor marker'
    )

    parser.add_argument(
        '-o', '--output',
        default='interval_variant_plot.html',
        help='Output HTML file (default: interval_variant_plot.html)'
    )
    
    parser.add_argument(
        '--title',
        default='QTL Variant Plot',
        help='Plot title'
    )
    
    parser.add_argument(
        '--pos-col',
        default='POS',
        help='Column name for position (default: POS)'
    )
    
    parser.add_argument(
        '--pval-col',
        default='P',
        help='Column name for p-value (default: P)'
    )
    
    parser.add_argument(
        '--delimiter',
        default='\t',
        help='File delimiter (default: tab)'
    )
    
    parser.add_argument(
        '--chrom-names',
        default=None,
        help='File with chromosome name mappings'
    )
    
    parser.add_argument(
        '--point-size',
        type=int,
        default=8,
        help='Point size (default: 8)'
    )
    
    parser.add_argument(
        '--dimensions',
        nargs=2,
        default=[1200, 400],
        help='Dimensions of plot (width, height)'
    )
    
    return parser


if __name__ == "__main__":
    main()