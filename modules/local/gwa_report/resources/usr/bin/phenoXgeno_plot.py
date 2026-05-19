#!/usr/bin/env python3
"""
Interactive Phenotype-Genotype Plot Generator using Plotly
Creates an HTML file with an interactive Phenotype-Genotype plot for GWAS results
"""

import argparse
import csv
import math

import plotly.graph_objects as go
import plotly.express as pe
from plotly.subplots import make_subplots
import numpy as np

def main():
    parser = make_arg_parser()
    args = parser.parse_args()
    
    markers = args.marker.split(',')
    if args.strains is not None:
        highlight_strains = set(args.strains.split(','))
    else:
        highlight_strains = None

    genotype_matrix, strains = load_genotype_matrix(args.genotypes, set(markers))
    traits = load_traits(args.traits)

    fig = plot_phenotype_genotype(genotype_matrix, traits, strains, markers, highlight_strains=highlight_strains,
                                  title=args.title, colors=args.colors, highlight_color=args.highlight_color)

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
    print(f"Phenotype-genotype plot saved to: {args.output}")
    
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

def load_traits(traits_fname):
    with open(traits_fname, 'r') as fs:
        header = fs.readline()
        data = np.array([tuple(line.rstrip().split("\t")[1:]) for line in fs],
                        dtype=np.dtype([('strain', 'U8'), ('value', np.float32)]))
    data = data[np.argsort(data['strain'])]
    return data

def plot_phenotype_genotype(genotype_matrix,
                            traits,
                            strains,
                            markers,
                            highlight_strains=None,
                            title='Phenotype-Genotype Plot',
                            colors=['#005AB5', '#1A85FF'],
                            highlight_color='#D41159'):

    fig = make_subplots(rows=(len(markers) - 1) // 6 + 1,
                        cols=min(len(markers), 6),
                        subplot_titles=markers,
                        vertical_spacing=0.075,
                        horizontal_spacing=0.03)
    
    for i, marker in enumerate(markers):
        index = np.where(genotype_matrix['marker'] == marker)[0][0]
        valid_strains = set([str(strains[x]) for x in range(len(strains))
                             if genotype_matrix['genotype'][index, x] != 1]).intersection(set([str(strain) for strain in traits['strain']]))
        valid_genotypes = genotype_matrix['genotype'][index, [x for x in range(len(strains)) if str(strains[x]) in valid_strains]]
        valid_traits = traits[[x for x in range(traits.shape[0]) if str(traits['strain'][x]) in valid_strains]]

        if highlight_strains is not None:
            highlight_indices = [i for i in range(valid_traits.shape[0]) if str(valid_traits['strain'][i]) in highlight_strains]
            h_traits = valid_traits[highlight_indices]
            h_genotypes = valid_genotypes[highlight_indices]
            normal_indices = [i for i in range(valid_traits.shape[0]) if str(valid_traits['strain'][i]) not in highlight_strains]
            valid_traits = valid_traits[normal_indices]
            valid_genotypes = valid_genotypes[normal_indices]

        ref = np.where(valid_genotypes == 0)[0]
        alt = np.where(valid_genotypes == 2)[0]

        if highlight_strains is not None:
            highlight_X = [int(str(h_genotypes[x]) == 0) + np.random.random()*0.1 - 0.05 + 2 * i for x in range(h_genotypes.shape[0])]
            fig.add_trace(go.Scatter(
                    x=highlight_X,
                    y=h_traits['value'],
                    text=[str(strain) for strain in h_traits['strain']],
                    textposition='top left',
                    textfont=dict(color=highlight_color),
                    mode='markers+text',
                    hoverinfo='text',
                    marker=dict(color=highlight_color),
                    showlegend=False,),
                row=i // 6 + 1,
                col=i % 6 + 1,)
            
        fig.add_trace(go.Box(y=valid_traits['value'][ref],
                             name='REF',
                             text=valid_traits['strain'][ref],
                             fillcolor='#FFFFFF',
                             line=dict(
                                 color='#000000',
                             ),
                             marker=dict(
                                 color='#000000',
                                opacity=0.5,
                             ),
                             hoverinfo='text',
                             hoveron='points',
                             boxpoints='all',
                             showlegend=False,
                             pointpos=0,),
                      row=i // 6 + 1,
                      col=i % 6 + 1,)
        fig.add_trace(go.Box(y=valid_traits['value'][alt],
                             name='ALT',
                             text=valid_traits['strain'][alt],
                             fillcolor='#FFFFFF',
                             line=dict(
                                 color='#000000',
                             ),
                             marker=dict(
                                 color='#000000',
                                opacity=0.5,
                             ),
                             hoverinfo='text',
                             hoveron='points',
                             boxpoints='all',
                             showlegend=False,
                             pointpos=0,),
                      row=i // 6 + 1,
                      col=i % 6 + 1,)

        if (i % 6 + 1) == 1:
            fig.update_yaxes(title_text="Phenotype Value", row=i // 6 + 1, col=i % 6 + 1)
        if i // 6 == len(markers) // 6 or (i // 6 + 1) * 6 + (i % 6) >= len(markers):
            fig.update_xaxes(title_text="Genotype", row=i // 6 + 1, col=i % 6 + 1)
        fig.update_xaxes(showline=True,
                         linewidth=2,
                         linecolor='black',
                         mirror=True,
                         row=i // 6 + 1,
                         col=i % 6 + 1,
                         ticks='outside',
                         tickvals=np.arange(2) + 2 * i,
                         ticktext=['REF', 'ALT'],
                         tickmode='array',)
        fig.update_yaxes(showline=True,
                         linewidth=2,
                         linecolor='black',
                         ticks='outside',
                         mirror=True,
                         row=i // 6 + 1,
                         col=i % 6 + 1)
    if title is not None:
        fig.update_layout(
            title=dict(
                text=title,
                font=dict(size=20),
                x=0.5
            ),
        )
        upper_margin=80
    else:
        upper_margin=40
    fig.update_layout(
        margin=dict(l=0, r=0, t=upper_margin, b=0),
    )
    return fig


def make_arg_parser():
    parser = argparse.ArgumentParser(
        description='Create an interactive Q-Q plot from GWAS results'
    )
    
    parser.add_argument(
        '-t', '--traits',
        required=True,
        help='Input phenotype traits file'
    )
    
    parser.add_argument(
        '-g', '--genotypes',
        required=True,
        help='Input genotype matrix file'
    )
    
    parser.add_argument(
        '-m', '--marker',
        required=True,
        help='Comma-separated list of markers to plot'
    )
    
    parser.add_argument(
        '-o', '--output',
        default='pxg_plot.html',
        help='Output HTML file (default: pxg_plot.html)'
    )
    
    parser.add_argument(
        '--title',
        default='Phenotype-Genotype Plot',
        help='Plot title'
    )
    
    parser.add_argument(
        '--strains',
        default=None,
        help='Comma-separated list of strains to highlight'
    )
    
    parser.add_argument(
        '--colors',
        nargs=2,
        default=['#005AB5', '#1A85FF'],
        help='List of colors (ref, alt)'
    )
    
    parser.add_argument(
        '--highlight-color',
        default='#D41159',
        help='Strain highlight color'
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