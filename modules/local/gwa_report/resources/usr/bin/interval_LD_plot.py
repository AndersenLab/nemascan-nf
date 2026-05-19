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
        snp_col=args.snp_col,
        delimiter=args.delimiter
    )
    markers = set([str(x) for x in data['marker']])
    
    print(f"Loaded {len(data)} SNPs")
    
    # Read genotype matrix data
    genotype_matrix = load_genotype_matrix(
        args.genotype_matrix,
        markers=markers,
        chrom_names=args.chrom_names
    )

    # Create Manhattan plot
    print("Creating Interval LD plot...")
    fig = create_interval_ld_plot(
        data,
        genotype_matrix,
        args.marker,
        title=args.title,
        point_size=args.point_size,
        gradient=args.gradient,
        peak_color=args.peak_color,
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
                'filename': 'interval_ld_plot',
                'height': args.dimensions[1],
                'width': args.dimensions[0],
                'scale': 2
            }
        }
    )
    
    print(f"Interval LD plot saved to: {args.output}")
    
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
    max_chromname = 0
    max_snp = 0
    with open(filename, 'r') as f:
        header = [x.upper() for x in f.readline().rstrip().split(delimiter)]
        chrom_col = header.index(chrom_col)
        pos_col = header.index(pos_col)
        pval_col = header.index(pval_col)
        snp_col = header.index(snp_col)
        
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
            
            pos = int(row[pos_col])
            pval = float(row[pval_col])
            
            # Skip invalid p-values
            if pval <= 0 or pval > 1:
                continue
            
            snp = row[snp_col]
            
            max_chromname = max(max_chromname, len(chrom))
            max_snp = max(max_snp, len(snp))
            data.append((chrom_num, chrom, pos, pval, -math.log10(pval), f"{chrom}:{pos}"))
                
    dtype = np.dtype([('chrom', np.int32), ('chrom_label', f'U{max_chromname}'), ('pos', np.int32),
                      ('pval', np.float32), ('log_pval', np.float32), ('marker', f'U{max_snp}')])
    data = np.array(data, dtype=dtype)
    data = data[np.lexsort((data['pos'], data['chrom']))]
    return data

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
    return data

def find_correlations(genotype_matrix, index):
    correlations = np.zeros(genotype_matrix.shape[0], np.float32)
    marker_genotypes = genotype_matrix['genotype'][index, :]
    for i in range(genotype_matrix.shape[0]):
        if i == index:
            continue
        correlations[i] = np.corrcoef(genotype_matrix['genotype'][i, :], marker_genotypes)[0, 1] ** 2
    return correlations

def create_interval_ld_plot(data,
                            genotype_matrix,
                            marker,
                            title="Interval LD Plot", 
                            point_size=4,
                            gradient='YlGnBu',
                            peak_color='#000000',
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
    
    marker_index = np.searchsorted(genotype_matrix['pos'], int(marker.split(':')[-1]))
    correlations = find_correlations(genotype_matrix, marker_index)

    start = np.searchsorted(data['pos'], genotype_matrix['pos'][0])
    indices = np.r_[np.arange(marker_index), np.arange(marker_index + 1, genotype_matrix.shape[0])]

    # Create figure
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=genotype_matrix['pos'][indices] / 1e6,
        y=data['log_pval'][indices + start],
        mode='markers',
        marker=dict(
            color=correlations[indices],
            colorscale=gradient,
            line=dict(
                color='black',
                width=1,
            ),
            showscale=True,
            size=point_size,
        ),
        text=[f"{genotype_matrix['marker'][x]} {correlations[x]:0.4f}" for x in indices],
        hoverinfo='text',
        showlegend=False,
    ))
    fig.add_trace(go.Scatter(
        x=[genotype_matrix['pos'][marker_index] / 1e6],
        y=[data['log_pval'][marker_index + start]],
        mode='markers',
        marker=dict(
            color=peak_color,
            size=point_size,
            line=dict(
                color='black',
                width=1,
            ),
        ),
        text=[f"{genotype_matrix['marker'][x]} 1.0" for x in [marker_index]],
        hoverinfo='text',
        showlegend=False,
    ))

    if title is not None:
        fig.update_layout(
            title=dict(
                text=title,
                font=dict(size=20),
                x=0.5
            ))
        upper_margin=40
    else:
        upper_margin=30
    fig.update_layout(
        xaxis=dict(
            title='Genomic Position (Mb)',
            showline=True,
            ticks='outside',
            linewidth=2,
            linecolor='black',
            mirror=True,
            showgrid=False,
            zeroline=False,
        ),
        yaxis=dict(
            title='-log<sub>10</sub>(P-value)',
            showline=True,
            ticks='outside',
            linewidth=2,
            linecolor='black',
            mirror=True,
            showgrid=False,
            zeroline=False,
        ),
        plot_bgcolor='white',
        hovermode='closest',
        margin=dict(l=30, r=0, t=upper_margin, b=0)
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
        '-g', '--genotype-matrix',
        required=True,
        help='Input genotype matrix file'
    )
    
    parser.add_argument(
        '-m', '--marker',
        required=True,
        help='Input anchor marker'
    )

    parser.add_argument(
        '-o', '--output',
        default='interval_ld_plot.html',
        help='Output HTML file (default: interval_ld_plot.html)'
    )
    
    parser.add_argument(
        '--title',
        default='QTL Interval LD Plot',
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
        '--snp-col',
        default='SNP',
        help='Column name for SNP identifier (default: SNP)'
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
        default=10,
        help='Point size (default: 10)'
    )
    
    parser.add_argument(
        '--gradient',
        default='YlGnBu',
        help='Color gradient for alternating chromosomes (default: )'
    )
    
    parser.add_argument(
        '--peak-color',
        nargs=3,
        default='#DC3220',
        help='Anchor peak color'
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