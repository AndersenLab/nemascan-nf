#!/usr/bin/env python3
"""
Interactive LD Plot Generator using Plotly
Creates an HTML file with an interactive LD plot for GWAS results
"""

import argparse
import math

import plotly.graph_objects as go
import plotly.express as pe
from plotly.subplots import make_subplots
import numpy as np

def main():
    parser = make_arg_parser()
    args = parser.parse_args()
    
    markers = args.markers.split(',')
    if len(markers) <= 1:
        print("Not enough SNPs (need more than 1) to calculate LD")
        return

    genotype_matrix, strains = load_genotype_matrix(args.genotypes, set(markers))

    LD_matrix = find_LD(genotype_matrix, markers)

    fig = plot_LD_matrix(LD_matrix, markers, title=args.title)

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
                'filename': 'pxg_plot',
                'height': args.dimensions[1],
                'width': args.dimensions[0],
                'scale': 2
            }
        }
    )
    print(f"LD plot saved to: {args.output}")
    
def load_genotype_matrix(gm_fname, markers=None, chrom_names=None):
    if chrom_names is not None:
        chrom_names = {line.rstrip().split("\t")[1]: line.split("\t")[0] for line in open(chrom_names)}
    else:
        chrom_names = {}

    geno = {'-1': 0, '1': 2, 'NA': 1}
    with open(gm_fname, 'r') as fs:
        header = fs.readline().rstrip().split("\t")
        strains = header[4:]
        data = []
        for line in fs:
            line = line.rstrip().split("\t")
            if line[0] in chrom_names:
                line[0] = chrom_names[line[0]]
            marker = f"{line[0]}:{line[1]}"
            if markers is not None and marker not in markers:
                continue
            data.append((marker, line[1], tuple([geno[x] for x in line[4:]])))
    markerlen = max([len(line[0]) for line in data])
    data = np.array(data, dtype=np.dtype([('marker', f"U{markerlen}"), ('pos', np.int32), ("genotype", np.uint8, (len(strains),))]))
    strains = np.array(strains, 'U8')
    order = np.argsort(strains)
    strains = strains[order]
    data['genotype'] = data['genotype'][:, order]
    return data, strains

def find_LD(genotype_matrix, markers):
    LD_matrix = np.ones((len(markers), len(markers)), np.float32)
    for i in range(len(markers) - 1):
        index1 = np.where(genotype_matrix['marker'] == markers[i])[0][0]
        valid1 = genotype_matrix['genotype'][index1, :] != 1
        for j in range(i + 1, len(markers)):
            index2 = np.where(genotype_matrix['marker'] == markers[j])[0][0]
            valid2 = genotype_matrix['genotype'][index2, :] != 1
            valid = np.where(np.logical_and(valid1, valid2))[0]
            LD_matrix[i, j] = np.corrcoef(genotype_matrix['genotype'][index1, valid],
                                          genotype_matrix['genotype'][index2, valid])[0, 1]
            LD_matrix[j, i] = LD_matrix[i, j]
    return LD_matrix

def plot_LD_matrix(LD_matrix,
                   markers,
                   title='Phenotype-Genotype Plot'):
    fig = go.Figure()
    fig.add_trace(go.Heatmap(
        z=LD_matrix,
        zmax=1,
        x=[str(x) for x in markers],
        y=[str(x) for x in markers],
        text=[[float(x) for x in row] for row in LD_matrix],
        texttemplate="%{text:0.3f}"
    ))
    fig.update_xaxes(showline=True,
                     linewidth=2,
                     linecolor='black',
                     mirror=True,)
    fig.update_yaxes(showline=True,
                     linewidth=2,
                     linecolor='black',
                     mirror=True,)
    if title is not None:
        fig.update_layout(
            title=dict(
                text=title,
                font=dict(size=20),
                x=0.5
            ),
        )
        upper_margin = 50
    else:
        upper_margin = 0
    fig.update_layout(
        margin=dict(l=0, r=0, t=upper_margin, b=0))
    return fig


def make_arg_parser():
    parser = argparse.ArgumentParser(
        description='Create an interactive LD plot from GWAS results'
    )
    
    parser.add_argument(
        '-g', '--genotypes',
        required=True,
        help='Input genotype matrix file'
    )
    
    parser.add_argument(
        '-m', '--markers',
        required=True,
        help='Comma-separated list of markers to plot'
    )
    
    parser.add_argument(
        '-o', '--output',
        default='ld_plot.html',
        help='Output HTML file (default: ld_plot.html)'
    )
    
    parser.add_argument(
        '--title',
        default='LD Plot',
        help='Plot title'
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