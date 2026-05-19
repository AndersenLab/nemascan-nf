#!/usr/bin/env python3
"""
Interactive Q-Q Plot Generator using Plotly
Creates an HTML file with an interactive Q-Q plot for GWAS results
"""

import argparse
import csv
import math

import plotly.graph_objects as go
import numpy as np



def main():
    parser = make_arg_parser()
    args = parser.parse_args()
    
    # Read GWAS data
    print(f"Reading GWAS data from {args.loco_input}...")
    loco_data = read_gwas_data(
        args.loco_input,
        chrom_names=args.chrom_names,
        chrom_col=args.chrom_col,
        pval_col=args.pval_col,
        delimiter=args.delimiter
    )

    print(f"Reading GWAS data from {args.inbred_input}...")
    inbred_data = read_gwas_data(
        args.inbred_input,
        chrom_names=args.chrom_names,
        chrom_col=args.chrom_col,
        pval_col=args.pval_col,
        delimiter=args.delimiter
    )

    # Create Manhattan plot
    print("Creating Q-Q plot...")
    fig, LOCO_lambda_gc, inbred_lambda_gc = create_qq_plot(
        loco_data,
        inbred_data,
        sig_cutoff=args.sig,
        colors=args.colors,
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
                'filename': 'qq_plot',
                'height': args.dimensions[1],
                'width': args.dimensions[0],
                'scale': 2
            }
        }
    )
    print(f"Q-Q plot saved to: {args.output}")
    
def read_gwas_data(filename,
                   chrom_col='CHR',
                   pval_col='P',
                   delimiter='\t',
                   chrom_names=None):
    """
    Read GWAS summary statistics from a file
    
    Parameters:
    -----------
    filename : str
        Path to the GWAS results file
    pval_col : str
        Column name for p-value
    delimiter : str
        File delimiter
    
    Returns:
    --------
    list of dicts containing GWAS data
    """
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

    data = []
    with open(filename, 'r') as f:
        header = [x.upper() for x in f.readline().rstrip().split(delimiter)]
        chrom_col = header.index(chrom_col)
        pval_col = header.index(pval_col)
        
        for line in f:
            row = line.rstrip().split("\t")

            chrom = row[chrom_col].replace('chr', '').replace('Chr', '')
            if chrom in chrom_names:
                chrom = chrom_names[chrom]
            
            # Handle X, Y, MT chromosomes
            if chrom.upper() in chrom_map:
                chrom_num = chrom_map[chrom.upper()]
            else:
                chrom_num = int(chrom)

            pval = float(row[pval_col])
            
            # Skip invalid p-values
            if pval <= 0 or pval > 1:
                continue
            data.append((chrom_num, chrom, pval))
            min_chromname = max(min_chromname, len(chrom))

    dtype = np.dtype([('chrom', np.int32), ('chrom_label', f'U{min_chromname}'), ('pval', np.float32)])
    data = np.array(data, dtype=dtype)
    data = data[np.argsort(data['chrom'])]
    return data

def create_qq_plot(LOCO_data,
                   inbred_data,
                   sig_cutoff,
                   title="Q-Q Plot",                   
                   colors=['#005AB5', '#1A85FF'],
                   sig_colors=['#D41159', '#D35FB7']):
    """
    Create an interactive Q-Q plot using Plotly
    
    Parameters:
    -----------
    observed_pvals : np.array
        Array containing GWAS p-values
    output_file : str
        Output HTML filename
    title : str
        Plot title
    dimnesions : list
        List of width and height for the plot
    """
    # Convert to -log10 scale
    n = LOCO_data.shape[0]
    LOCO_observed_pvals = np.copy(LOCO_data['pval'])
    LOCO_observed_log = -np.log10(LOCO_observed_pvals)

    n = LOCO_data.shape[0]
    inbred_observed_pvals = np.copy(inbred_data['pval'])
    inbred_expected_pvals = (np.arange(n) + 0.5) / n
    inbred_observed_log = -np.log10(inbred_observed_pvals)
    inbred_expected_log = -np.log10(inbred_expected_pvals)
    
    # Calculate genomic inflation factor (lambda)
    median_chi2_LOCO_observed = np.median(-2 * np.log(10) * LOCO_observed_log)
    median_chi2_inbred_observed = np.median(-2 * np.log(10) * inbred_observed_log)
    median_chi2_expected = 0.4549364  # median of chi-squared with df=1
    LOCO_lambda_gc = median_chi2_LOCO_observed / median_chi2_expected
    inbred_lambda_gc = median_chi2_inbred_observed / median_chi2_expected
    
    LOCO_chr_indices = np.r_[0, np.where(LOCO_data['chrom'][1:] != LOCO_data['chrom'][:-1])[0] + 1, LOCO_data.shape[0]]
    inbred_chr_indices = np.r_[0, np.where(inbred_data['chrom'][1:] != inbred_data['chrom'][:-1])[0] + 1, inbred_data.shape[0]]
    maxN = max(np.amax(LOCO_chr_indices[1:] - LOCO_chr_indices[:-1]),
               np.amax(inbred_chr_indices[1:] - inbred_chr_indices[:-1]))
    maxX = -np.log10(0.5 / maxN)

    chrom_labels = LOCO_data['chrom_label'][LOCO_chr_indices[1:] - 1]
    chromosomes = LOCO_data['chrom'][LOCO_chr_indices[1:] - 1]
    chrN = chrom_labels.shape[0]
    spacer = maxX * 0.1
    chrom_ranges = np.zeros((chrN, 2), np.float32)
    chrom_ranges[:, 1] = np.arange(1, chrN + 1) * maxX
    chrom_ranges[1:, 0] = chrom_ranges[:-1, 1]
    chrom_ranges += np.arange(chrN).reshape(-1, 1) * spacer
    
    # Create figure
    fig = go.Figure()  
    maxY = max(np.amax(inbred_observed_log), np.amax(LOCO_observed_log)) * 1.01
    X = {'LOCO_nonsig': [], 'LOCO_sig': [], 'Inbred_nonsig': [], 'Inbred_sig': []}
    Y = {'LOCO_nonsig': [], 'LOCO_sig': [], 'Inbred_nonsig': [], 'Inbred_sig': []}

    for i, chrom in enumerate(chromosomes):
        fig.add_trace(go.Scatter(
            x=[chrom_ranges[i, 0], chrom_ranges[i, 0], chrom_ranges[i, 1], chrom_ranges[i, 1], chrom_ranges[i, 0]],
            y=[0, maxY, maxY, 0, 0],
            mode='lines',
            name=f'{chrom_labels[i]} box',
            showlegend=False,
            line=dict(color="black"),
        ))
        
        s, e = LOCO_chr_indices[i:i+2]
            
        # Add diagonal line (y = x)
        fig.add_trace(go.Scatter(
            x=[chrom_ranges[i, 0], chrom_ranges[i, 1]],
            y=[0, chrom_ranges[i, 1] - chrom_ranges[i, 0]],
            mode='lines',
            line=dict(color='black', width=1, dash='solid'),
            name='Expected (y=x)',
            hoverinfo='name',
            showlegend=False,
        ))
    
        n = e - s
        LOCO_observed_chr_pvals = LOCO_observed_pvals[s:e]
        LOCO_observed_chr_pvals.sort()
        LOCO_observed_chr_log = LOCO_observed_log[s:e]
        LOCO_observed_chr_log.sort()
        LOCO_observed_chr_log = LOCO_observed_chr_log[::-1]
        LOCO_expected_pvals = (np.arange(n) + 0.5) / n
        LOCO_expected_chr_log = -np.log10(LOCO_expected_pvals) + chrom_ranges[i, 0]

        s, e = inbred_chr_indices[i:i+2]
        n = e - s
        inbred_observed_chr_pvals = inbred_observed_pvals[s:e]
        inbred_observed_chr_pvals.sort()
        inbred_observed_chr_log = inbred_observed_log[s:e]
        inbred_observed_chr_log.sort()
        inbred_observed_chr_log = inbred_observed_chr_log[::-1]
        inbred_expected_pvals = (np.arange(n) + 0.5) / n
        inbred_expected_chr_log = -np.log10(inbred_expected_pvals) + chrom_ranges[i, 0]

        nonsig = np.where(LOCO_observed_chr_pvals > sig_cutoff)
        X['LOCO_nonsig'] += list(LOCO_expected_chr_log[nonsig])
        Y['LOCO_nonsig'] += list(LOCO_observed_chr_log[nonsig])
        sig = np.where(LOCO_observed_chr_pvals <= sig_cutoff)
        X['LOCO_sig'] += list(LOCO_expected_chr_log[sig])
        Y['LOCO_sig'] += list(LOCO_observed_chr_log[sig])

        nonsig = np.where(inbred_observed_chr_pvals > sig_cutoff)
        X['Inbred_nonsig'] += list(inbred_expected_chr_log[nonsig])
        Y['Inbred_nonsig'] += list(inbred_observed_chr_log[nonsig])
        sig = np.where(inbred_observed_chr_pvals <= sig_cutoff)
        X['Inbred_sig'] += list(inbred_expected_chr_log[sig])
        Y['Inbred_sig'] += list(inbred_observed_chr_log[sig])

    # Add observed vs expected points
    fig.add_trace(go.Scattergl(
        x=X['LOCO_nonsig'],
        y=Y['LOCO_nonsig'],
        mode='markers',
        marker=dict(
            size=8,
            color=colors[0],
            opacity=0.6
        ),
        name='LOCO',
        hovertemplate=(
            'Expected -log10(P): %{x:.2f}<br>'
            'Observed -log10(P): %{y:.2f}<br>'
            '<extra></extra>'
        ),
        showlegend=True,
    ))

    fig.add_trace(go.Scattergl(
        x=X['Inbred_nonsig'],
        y=Y['Inbred_nonsig'],
        mode='markers',
        marker=dict(
            size=8,
            color=colors[1],
            opacity=0.6
        ),
        name='Inbred',
        hovertemplate=(
            'Expected -log10(P): %{x:.2f}<br>'
            'Observed -log10(P): %{y:.2f}<br>'
            '<extra></extra>'
        ),
        showlegend=True,
    ))

    fig.add_trace(go.Scattergl(
        x=X['LOCO_sig'],
        y=Y['LOCO_sig'],
        mode='markers',
        marker=dict(
            size=8,
            color=sig_colors[0],
            opacity=0.6
        ),
        name='LOCO-significant',
        hovertemplate=(
            'Expected -log10(P): %{x:.2f}<br>'
            'Observed -log10(P): %{y:.2f}<br>'
            '<extra></extra>'
        ),
        showlegend=True,
    ))

    fig.add_trace(go.Scattergl(
        x=X['Inbred_sig'],
        y=Y['Inbred_sig'],
        mode='markers',
        marker=dict(
            size=8,
            color=sig_colors[1],
            opacity=0.6
        ),
        name='Inbred-significant',
        hovertemplate=(
            'Expected -log10(P): %{x:.2f}<br>'
            'Observed -log10(P): %{y:.2f}<br>'
            '<extra></extra>'
        ),
        showlegend=True,
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
    chrom_length = int(chrom_ranges[0, 1] - chrom_ranges[0, 0])
    for i in range(chrN):
        xtickvals += list(np.arange(chrom_length + 1) + chrom_ranges[i, 0])
        xticktext += [f"{x}" for x in (np.arange(chrom_length + 1))]

    if title is not None:
        fig.update_layout(
            title=dict(
                text=f"{title}",
                font=dict(size=18),
                x=0.5
            ))
        upper_margin=80
    else:
        upper_margin=40
    fig.update_layout(
        xaxis=dict(
            title='Expected -log<sub>10</sub>(P-value)',
            showgrid=False,
            tickmode='array',
            ticks='outside',
            tickvals=xtickvals,
            ticktext=xticktext,
            zeroline=False,
            range=[0, chrom_ranges[-1, 1]],
        ),
        yaxis=dict(
            title='Observed -log<sub>10</sub>(P-value)',
            showgrid=False,
            zeroline=False,
            ticks='outside',
            range=[0, maxY],
        ),
        plot_bgcolor='white',
        hovermode='closest',
        showlegend=True,
        margin=dict(l=60, r=40, t=upper_margin, b=80)
    )
    
    return fig, LOCO_lambda_gc, inbred_lambda_gc
        
def make_arg_parser():
    parser = argparse.ArgumentParser(
        description='Create an interactive Q-Q plot from GWAS results'
    )
    
    parser.add_argument(
        '-l', '--loco-input',
        required=True,
        help='Input LOCO GWAS summary statistics file'
    )
    
    parser.add_argument(
        '-i', '--inbred-input',
        required=True,
        help='Input inbred GWAS summary statistics file'
    )
    
    parser.add_argument(
        '-o', '--output',
        default='qq_plot.html',
        help='Output HTML file (default: qq_plot.html)'
    )
    
    parser.add_argument(
        '--title',
        default='Q-Q Plot',
        help='Plot title'
    )
    
    parser.add_argument(
        '--chrom-names',
        default=None,
        help='File with chromosome name mappings'
    )
    
    parser.add_argument(
        '--chrom-col',
        default='CHR',
        help='Column name for chromosome (default: CHR)'
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
        '--colors',
        nargs=2,
        default=['#005AB5', '#1A85FF'],
        help='List of colors (LOCO, Inbred)'
    )
    
    parser.add_argument(
        '--sig-colors',
        nargs=2,
        default=['#D41159', '#DC3220'],
        help='List of significance colors (LOCO, Inbred)'
    )

    parser.add_argument(
        '--dimensions',
        nargs=2,
        default=[1200, 400],
        help='Dimensions of plot (width, height)'
    )
    
    parser.add_argument(
        '--sig',
        type=float,
        default=5e-8,
        help='Significance threshold (default: 5e-8)'
    )
    
    return parser


if __name__ == "__main__":
    main()