#!/usr/bin/env python3
"""
Interactive Manhattan Plot Generator using Plotly
Creates an HTML file with an interactive Manhattan plot for GWAS results
"""

import argparse
import math

import plotly.graph_objects as go
import numpy as np


def main():
    parser = make_arg_parser()
    args = parser.parse_args()
    
    # Read GWAS data
    print(f"Reading GWAS data from {args.input}...")
    data = read_gwas_data(
        args.input,
        chrom_names=args.chrom_names,
        chrom_col=args.chrom_col,
        pos_col=args.pos_col,
        pval_col=args.pval_col,
        delimiter=args.delimiter
    )
    
    print(f"Loaded {len(data)} SNPs")
    
    # Create Manhattan plot
    print("Creating Manhattan plot...")
    fig = create_manhattan_plot(
        data,
        title=args.title,
        bf_sig=args.bf_sig,
        eigen_sig=args.eigen_sig,
        user_sig=args.user_sig,
        point_size=args.point_size,
        sig_colors=args.sig_colors,
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
                'filename': 'manhattan_plot',
                'height': args.dimensions[1],
                'width': args.dimensions[0],
                'scale': 2
            }
        }
    )
    
    print(f"Manhattan plot saved to: {args.output}")
    
def read_gwas_data(filename, chrom_names=None, chrom_col='CHR', pos_col='POS',
                   pval_col='P', snp_col='SNP', delimiter='\t'):
    """
    Read GWAS summary statistics from a file
    
    Parameters:
    -----------
    filename : str
        Path to the GWAS results file
    chrom_col : str
        Column name for chromosome
    pos_col : str
        Column name for position
    pval_col : str
        Column name for p-value
    snp_col : str
        Column name for SNP identifier
    delimiter : str
        File delimiter
    
    Returns:
    --------
    list of dicts containing GWAS data
    """
    data = []
    if chrom_names is not None:
        chrom_names = {line.rstrip().split("\t")[1]: line.split("\t")[0] for line in open(chrom_names)}
    else:
        chrom_names = {}

    chrom_map = {
        'I'    : 1,
        'II'   : 2,
        'III'  : 3,
        'IV'   : 4,
        'V'    : 5,
        'X'    : 23,
        'Y'    : 24,
        'MT'   : 25,
        'M'    : 25,
        'MTDNA': 25
    }
    min_chromname = 0
    min_snp = 0
    with open(filename, 'r') as f:
        header = [x.upper() for x in f.readline().rstrip().split(delimiter)]
        chrom_col = header.index(chrom_col)
        pos_col = header.index(pos_col)
        pval_col = header.index(pval_col)
        
        for line in f:
            row = line.rstrip().split("\t")
            # try:
            chrom = row[chrom_col].replace('chr', '').replace('Chr', '')
            if chrom in chrom_names:
                chrom = chrom_names[chrom]
            
            # Handle X, Y, MT chromosomes
            if chrom.upper() in chrom_map:
                chrom_num = chrom_map[chrom.upper()]
            else:
                chrom_num = int(chrom)
            
            pos = int(row[pos_col])
            ref = row[pos_col + 1]
            alt = row[pos_col + 2]
            pval = float(row[pval_col])
            
            # Skip invalid p-values
            if pval <= 0 or pval > 1:
                continue
            
            snp = f"{chrom}:{pos}"
            
            min_chromname = max(min_chromname, len(chrom))
            min_snp = max(min_snp, len(snp))
            data.append((chrom_num, chrom, pos, ref, alt, pval, -math.log10(pval), snp))
                
    dtype = np.dtype([('chrom', np.int32), ('chrom_label', f'U{min_chromname}'), ('pos', np.int32),
                      ('ref', 'U1'), ('alt', 'U1'), ('pval', np.float32), ('log_pval', np.float32),
                      ('marker', f'U{min_snp}')])
    data = np.array(data, dtype=dtype)
    data = data[np.lexsort((data['pos'], data['chrom']))]
    return data

def create_manhattan_plot(data, title="Manhattan Plot", 
                          bf_sig=5e-8, 
                          eigen_sig=5e-6,
                          user_sig=None,
                          output_file="manhattan_plot.html",
                          point_size=4,
                          sig_colors=None,
                          dimensions=[1200, 400],
                          ):
    """
    Create an interactive Manhattan plot using Plotly
    
    Parameters:
    -----------
    data : list
        List of dictionaries containing GWAS data
    title : str
        Plot title
    bf_sig : float
        Bonferonni significance threshold (default: 5e-8)
    eigen_sig : float
        Eigen-adjusted significance threshold (default: 5e-6)
    user_sig : float
        User-specified significance threshold (default: None)
    output_file : str
        Output HTML filename
    point_size : int
        Size of the data points
    colors : list
        List of two colors for alternating chromosomes
    sig_colors : list
        List of colors to user for different significance threholds
    """
    
    # Calculate cumulative positions
    chr_indices = np.r_[0, np.where(data['chrom'][1:] != data['chrom'][:-1])[0] + 1, data.shape[0]]
    chrom_lengths = data['pos'][chr_indices[1:] - 1]

    chromosomes = data['chrom'][chr_indices[1:] - 1]
    chrN = chromosomes.shape[0]
    chrom_labels = data['chrom_label'][chr_indices[1:] - 1]
    cumulative_lengths = np.r_[0, np.cumsum(chrom_lengths)]
    spacer = cumulative_lengths[-1] * 0.1 / chromosomes.shape[0]
    chrom_ranges = np.zeros((chrN, 2), np.int32)
    chrom_ranges[:, 0] = cumulative_lengths[:-1] + np.arange(chrN) * spacer
    chrom_ranges[:, 1] = cumulative_lengths[1:] + np.arange(chrN) * spacer

    sig_order = [(-math.log10(bf_sig), 'Bonferonni'), (-math.log10(eigen_sig), 'Eigen')] 
    if user_sig is not None:
        sig_order.append((-math.log10(user_sig), 'User'))
    sig_order.sort()
    sig_order.reverse()

    maxY = max(sig_order[0][0], np.amax(data['log_pval'])) * 1.01
    # Create figure
    fig = go.Figure()
    
    # Add traces for each chromosome
    for i, chrom in enumerate(chromosomes):
        s, e = chr_indices[i:i+2]
        chrom_data = data[s:e]

        if chrom != 25:
            fig.add_trace(go.Scatter(
                x=[chrom_ranges[i, 0], chrom_ranges[i, 0], chrom_ranges[i, 1], chrom_ranges[i, 1], chrom_ranges[i, 0]],
                y=[0, maxY, maxY, 0, 0],
                mode='lines',
                name=f'{chrom_labels[i]} box',
                showlegend=False,
                line=dict(color="black"),
            ))
        
        # Separate highlighted and non-highlighted SNPs
        nonsig_data = chrom_data[np.where(chrom_data['log_pval'] < sig_order[-1][0])]

        # Regular points
        x_vals = nonsig_data['pos'] + chrom_ranges[i, 0]
        y_vals = nonsig_data['log_pval']
        hover_text = [
            f"SNP: {nonsig_data['marker'][d]}<br>"
            f"Chr: {nonsig_data['chrom_label'][d]}<br>"
            f"Pos: {nonsig_data['pos'][d]:,}<br>"
            f"P-value: {nonsig_data['pval'][d]:.2e}<br>"
            f"-log10(P): {nonsig_data['log_pval'][d]:.2f}"
            for d in range(nonsig_data.shape[0])
        ]
        
        fig.add_trace(go.Scattergl(
            x=x_vals,
            y=y_vals,
            mode='markers',
            marker=dict(
                size=point_size,
                color='black',
                opacity=0.3
            ),
            text=hover_text,
            hoverinfo='text',
            name=f'Chr {chrom_labels[i]}',
            showlegend=False
        ))
        
        # Significant points
        for j in range(len(sig_order)):
            if i == 0:
                name = f'{sig_order[j][1]} (P={bf_sig:.0e})'
            else:
                name = f'{sig_order[j][1]} {chrom_labels[i]}'
            fig.add_trace(go.Scatter(
                x=[chrom_ranges[i, 0], chrom_ranges[i, 1]],
                y=[sig_order[j][0], sig_order[j][0]],
                mode='lines',
                line=dict(color='black', width=1, dash=['solid', 'dash', 'dotted'][j]),
                name=name,
                hoverinfo='name',
                showlegend=bool(i == 0)
            ))

            if j == 0:
                sig_data = chrom_data[np.where(chrom_data['log_pval'] >= sig_order[j][0])]
            else:
                sig_data = chrom_data[np.where(np.logical_and(chrom_data['log_pval'] < sig_order[j - 1][0], chrom_data['log_pval'] >= sig_order[j][0]))]
            if sig_data.shape[0] == 0:
                continue
            x_vals = sig_data['pos'] + cumulative_lengths[i]
            y_vals = sig_data['log_pval']
            hover_text = [
                f"SNP: {sig_data['marker'][d]}<br>"
                f"Chr: {sig_data['chrom_label'][d]}<br>"
                f"Pos: {sig_data['pos'][d]:,}<br>"
                f"P-value: {sig_data['pval'][d]:.2e}<br>"
                f"-log10(P): {sig_data['log_pval'][d]:.2f}"
                for d in range(sig_data.shape[0])
            ]
            
            fig.add_trace(go.Scattergl(
                x=x_vals,
                y=y_vals,
                mode='markers',
                marker=dict(
                    size=point_size,
                    color=sig_colors[j],
                    opacity=0.8
                ),
                text=hover_text,
                hoverinfo='text',
                name=f'{sig_order[j][1]} {chrom_labels}',
                showlegend=False
            ))
    
    for i, chrom in enumerate(chrom_labels):
        fig.add_annotation(
            x=np.mean(chrom_ranges[i, :]),
            y=maxY,
            text=chrom,
            showarrow=False,
            xanchor='center',
            yanchor='bottom',
            font=dict(
                size=14,
                weight='bold',
            )
        )

    # Update layout
    xtickvals = []
    xticktext = []
    for i in range(chrN):
        xtickvals += list(np.arange(1, int(chrom_lengths[i] / 5e6) + 1) * 5e6 + chrom_ranges[i, 0])
        xticktext += [f"{x}" for x in (np.arange(1, int(chrom_lengths[i] / 5e6) + 1) * 5)]
    if title is not None:
        fig.update_layout(
            title=dict(
                text=title,
                font=dict(size=20),
                x=0.5
            ))
        upper_margin = 80
    else:
        upper_margin = 40
    fig.update_layout(
        xaxis=dict(
            title='Genomic Position (Mb)',
            tickmode='array',
            ticks='outside',
            tickvals=xtickvals,
            ticktext=xticktext,
            showgrid=False,
            zeroline=False,
            rangeslider=None,
            range=[0, chrom_ranges[-1, 1] * 1.035],
        ),
        yaxis=dict(
            title='-log<sub>10</sub>(P-value)',
            showgrid=False,
            ticks='outside',
            zeroline=False,
            range=[0, maxY],
        ),
        plot_bgcolor='white',
        hovermode='closest',
        legend=dict(
            orientation='h',
            yanchor='bottom',
            y=1.08,
            xanchor='right',
            x=1
        ),
        margin=dict(l=30, r=0, t=upper_margin, b=0)
    )
    
    # Add range slider for zooming
    fig.update_layout(
        xaxis=dict(
            rangeslider=dict(visible=False),
            type='linear'
        )
    )
    
    return fig

def make_arg_parser():
    parser = argparse.ArgumentParser(
        description='Create an interactive Manhattan plot from GWAS results'
    )
    
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input GWAS summary statistics file'
    )
    
    parser.add_argument(
        '-o', '--output',
        default='manhattan_plot.html',
        help='Output HTML file (default: manhattan_plot.html)'
    )
    
    parser.add_argument(
        '--title',
        default='Manhattan Plot',
        help='Plot title'
    )
    
    parser.add_argument(
        '--chrom-col',
        default='CHR',
        help='Column name for chromosome (default: CHR)'
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
        '--bf-sig',
        type=float,
        default=5e-8,
        help='Bonferonni significance threshold (default: 5e-8)'
    )
    
    parser.add_argument(
        '--eigen-sig',
        type=float,
        default=5e-6,
        help='Eigen-corrected significance threshold (default: 5e-6)'
    )
    
    parser.add_argument(
        '--user-sig',
        type=float,
        default=None,
        help='User-specified significance threshold (default: None)'
    )
    
    parser.add_argument(
        '--point-size',
        type=int,
        default=6,
        help='Point size (default: 6)'
    )
    
    parser.add_argument(
        '--sig-colors',
        nargs=3,
        default=['#D41159', '#DC3220', '#D35FB7'],
        help='List of significance colors'
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