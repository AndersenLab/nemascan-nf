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
    
    # Read genotype matrix data
    print(f"Reading genotype matrix data from {args.input}...")
    data, strains = load_genotype_matrix(
        args.input,
        markers = set([args.marker]),
        chrom_names=args.chrom_names
    )
    
    print(f"Loaded genotype matrix")
    
    # Read HDR data
    hdrs = load_HDRs(args.hdr)

    # Create HDR plot
    print("Creating HDR plot...")
    fig = create_HDR_plot(
        data['genotype'][np.where(data['marker'] == args.marker)[0][0], :],
        HDRs=hdrs,
        strains=strains,
        marker=args.marker,
        start=args.start,
        stop=args.stop,
        title=args.title,
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
                'filename': 'HDR_plot',
                'height': args.dimensions[1],
                'width': args.dimensions[0],
                'scale': 2
            }
        }
    )
    
    print(f"HDR plot saved to: {args.output}")
   
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
    return data, strains

def load_HDRs(fname):
    data = {}
    with gzip.open(fname, 'rb') as fs:
        for line in fs:
            line = line.decode('utf8').rstrip().split("\t")
            data.setdefault(line[0], {})
            data[line[0]].setdefault(line[3], [])
            data[line[0]][line[3]].append((int(line[1]), int(line[2])))
    for chrom in data.keys():
        for strain in data[chrom].keys():
            data[chrom][strain] = np.array(data[chrom][strain])
            data[chrom][strain] = data[chrom][strain][np.argsort(data[chrom][strain][:, 0])]
    return data

def create_HDR_plot(data,
                    HDRs,
                    strains,
                    marker,
                    start,
                    stop,
                    title="HDR Plot", 
                ):
    """
    Create an interactive HDR plot using Plotly
    
    Parameters:
    -----------
    data : np.array
        Array of strain genotypes
    HDRs : dict
        Dictionary of dictionaries of np.arrays keyed by chroms and strains, respectively of HDR bounds
    strains : list
        List of strain names
    marker : str
        Anchor marker
    start: int
        Start coordinate of region to plot
    stop: int
        Stop coordinate of region to plot
    title : str
        Plot title
    """
    
    marker_pos = int(marker.split(':')[-1]) / 1e6
    chrom = marker.split(':')[0]
    start /= 1e6
    stop /= 1e6
    ref_strains = []
    alt_strains = []
    for i in range(len(strains)):
        if data[i] == 2:
            ref_strains.append(strains[i])
        else:
            alt_strains.append(strains[i])

    # Create figure
    fig = make_subplots()
    divergent = 0
    fig.add_trace(go.Scatter(
        x=[start, stop, stop, start, start],
        y=[0, 0, len(ref_strains), len(ref_strains), 0],
        fillcolor='#bbbbbb',
        mode='lines',
        fill='toself',
        line=dict(color='#bbbbbb', width=0),
        name='Not Divergent',
        showlegend=True,
        hoverinfo='skip',
    ))
    fig.add_trace(go.Scatter(
        x=[start, stop, stop, start, start],
        y=[len(ref_strains) + 3, len(ref_strains) + 3, len(strains) + 3, len(strains) + 3, len(ref_strains) + 3],
        fillcolor='#bbbbbb',
        mode='lines',
        fill='toself',
        line=dict(color='#bbbbbb', width=0),
        showlegend=False,
        hoverinfo='skip',
    ))
    fig.add_trace(go.Scatter(
        x=[start, stop, stop, start, start],
        y=[0, 0, len(ref_strains), len(ref_strains), 0],
        mode='lines',
        line=dict(color='black', width=1),
        showlegend=False,
        hoverinfo='skip'))
    fig.add_trace(go.Scatter(
        x=[start, stop, stop, start, start],
        y=[len(ref_strains) + 3, len(ref_strains) + 3, len(strains) + 3, len(strains) + 3, len(ref_strains) + 3],
        mode='lines',
        line=dict(color='black', width=1),
        showlegend=False,
        hoverinfo='skip'))
    for i in range(len(ref_strains)):
        if ref_strains[i] in HDRs[chrom]:
            for j in range(HDRs[chrom][ref_strains[i]].shape[0]):
                rstart, rstop = HDRs[chrom][ref_strains[i]][j, :]
                text = f"{ref_strains[i]} {chrom}:{rstart}-{rstop}"
                rstart = rstart / 1e6
                rstop = rstop / 1e6
                if rstart > stop or rstop < start:
                    continue
                rstart = max(rstart, start)
                rstop = min(rstop, stop)
                fig.add_trace(go.Scatter(
                    x=[rstart, rstop, rstop, rstart, rstart],
                    y=[i, i, i + 1, i + 1, i],
                    mode='lines',
                    fill='toself',
                    fillcolor='#cc2181',
                    line=dict(color='#cc2181', width=0),
                    name='Divergent',
                    showlegend=divergent == 0,
                    text=text,
                    hoverinfo='text',
                    hoveron='fills',))
                divergent += 1
    for i in range(1, len(ref_strains)):
        fig.add_trace(go.Scatter(
            x=[start, stop],
            y=[i, i],
            mode='lines',
            line=dict(
                color='black',
                width=0.5
            ),
            showlegend=False,
            hoverinfo='skip'))

    for i in range(len(alt_strains)):
        y = i + len(ref_strains) + 3
        if alt_strains[i] in HDRs[chrom]:
            for j in range(HDRs[chrom][alt_strains[i]].shape[0]):
                rstart, rstop = HDRs[chrom][alt_strains[i]][j, :]
                text = f"{alt_strains[i]} {chrom}:{rstart}-{rstop}"
                rstart = rstart / 1e6
                rstop = rstop / 1e6
                if rstart > stop or rstop < start:
                    continue
                rstart = max(rstart, start)
                rstop = min(rstop, stop)
                fig.add_trace(go.Scatter(
                    x=[rstart, rstop, rstop, rstart, rstart],
                    y=[y, y, y + 1, y + 1, y],
                    mode='lines',
                    fill='toself',
                    fillcolor='#cc2181',
                    line=dict(color='#cc2181', width=0),
                    showlegend=False,
                    text=text,
                    hoverinfo='text',
                    hoveron='fills'))
    for i in range(1, len(alt_strains)):
        y = i + len(ref_strains) + 3
        fig.add_trace(go.Scatter(
            x=[start, stop],
            y=[y, y],
            mode='lines',
            line=dict(
                color='black',
                width=0.5
            ),
            showlegend=False,
            hoverinfo='skip'))
        
    fig.add_trace(go.Scatter(
        x=[marker_pos, marker_pos],
        y=[0, len(ref_strains)],
        mode='lines',
        line=dict(
            color='black',
            width=2,
            dash='dash',
        ),
        showlegend=False,
        hoverinfo='skip'))
        
    fig.add_trace(go.Scatter(
        x=[marker_pos, marker_pos],
        y=[len(ref_strains) + 3, len(strains) + 3],
        mode='lines',
        line=dict(
            color='black',
            width=2,
            dash='dash',
        ),
        showlegend=False,
        hoverinfo='skip'))

    fig.update_yaxes(
        tickvals=np.arange(len(strains) + 3) + 0.5,
        ticktext=(ref_strains + ([""] * 3) + alt_strains),
        title='Strain',
        secondary_y=False,
        range=[0, len(strains) + 3],
        tickfont=dict(
            size=8,
        ),)
    fig.update_xaxes(
        ticks='outside',
        title='Genomic Position (Mb)',
        range=[start, stop],
        )
    
    fig.add_annotation(
        x=stop,
        y=len(ref_strains) / 2,
        text=f"REF ({len(ref_strains)})",
        showarrow=False,
        textangle=90,
        xanchor='left',
        )
    fig.add_annotation(
        x=stop,
        y=len(ref_strains) + 3 + len(alt_strains) / 2,
        text=f"ALT ({len(alt_strains)})",
        showarrow=False,
        textangle=90,
        xanchor='left',
        )

    fig.update_layout(
        title=dict(
            text=title,
            font=dict(size=20),
            x=0.5
        ),
        legend=dict(
            yanchor='bottom',
            xanchor='right',
            y=1.005,
            x=0.995
        ),
        margin=dict(l=0, r=20, t=80, b=0)
    )
    return fig

def make_arg_parser():
    parser = argparse.ArgumentParser(
        description='Create an interactive hyper-divergent region plot from genotype matrix'
    )
    
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input genotype matrix file'
    )
    
    parser.add_argument(
        '--hdr',
        required=True,
        help='Input HDR file'
    )
    
    parser.add_argument(
        '-m', '--marker',
        required=True,
        help='Input anchor marker'
    )

    parser.add_argument(
        '-o', '--output',
        default='HDR_plot.html',
        help='Output HTML file (default: HDR_plot.html)'
    )
    
    parser.add_argument(
        '--title',
        default='Hyper-divergent Regions Plot',
        help='Plot title'
    )
    
    parser.add_argument(
        '--chrom-names',
        default=None,
        help='File with chromosome name mappings'
    )

    parser.add_argument(
        '--start',
        required=True,
        type=int,
        help='Region start coordinate'
    )

    parser.add_argument(
        '--stop',
        required=True,
        type=int,
        help='Region stop coordinate'
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