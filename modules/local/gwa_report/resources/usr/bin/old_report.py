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

import numpy as np

from manhattan_plot import read_gwas_data, create_manhattan_plot
from qq_plot import create_qq_plot
from phenoXgeno_plot import load_genotype_matrix, load_traits, plot_phenotype_genotype
from LD_plot import find_LD, plot_LD_matrix
from interval_LD_plot import find_correlations, create_interval_ld_plot
from interval_variant_plot import load_annotations, load_genes, create_interval_variant_plot
from HDR_plot import load_HDRs, create_HDR_plot
from haplotype_plot import load_haplotypes, create_haplotype_plot


def main():
    parser = make_arg_parser()
    args = parser.parse_args()

    broad_gwa = {}
    for method, gwa_fname, qtl_fname in args.gwas:
        if method == 'loco':
            pos_col = "BP"
        else:
            pos_col = "POS"
        gwa = read_gwas_data(gwa_fname)
        qtl = [line.rstrip().split("\t") for line in open(qtl_fname)]
        broad_gwa[method] = [gwa, qtl]
        
    


    # loco_bf_sig = args.alpha / np.where(np.logical_not(np.isnan(loco_gwa['log_pval'])))[0].shape[0]
    # inbred_bf_sig = args.alpha / np.where(np.logical_not(np.isnan(inbred_gwa['log_pval'])))[0].shape[0]
    # eigen_sig = args.alpha / args.eigen
    # if args.sig_method == "BF":
    #     sig_cutoff = loco_bf_sig
    # elif args.sig_method == "EIGEN":
    #     sig_cutoff = eigen_sig
    # else:
    #     sig_cutoff = args.user_sig

    # with open(args.loco_qtl) as fs:
    #     QTL_header = [fs.readline().rstrip().split("\t")]
    # loco_QTLs = [[item for i, item in enumerate(line.rstrip().split("\t")) if QTL_header[0][i] != 'peak_id'] for line in open(args.loco_qtl)]
    # inbred_QTLs = [[item for i, item in enumerate(line.rstrip().split("\t")) if QTL_header[0][i] != 'peak_id'] for line in open(args.inbred_qtl)]
    # QTL_header = loco_QTLs[0]
    # QTLs = list(loco_QTLs[1:]) + list(inbred_QTLs[1:])
    # QTLs.sort(key=lambda row: (row[0], int(row[3]), row[6]))
    # all_markers = set([row[1] for row in QTLs])

    # genotype_matrix, strains = load_genotype_matrix(args.genotypes, all_markers)

    # traits = load_traits(args.traits)

    # if args.genes is not None:
    #     genes = load_genes(args.genes)
    # else:
    #     genes = None

    # if args.hdr is not None:
    #     hdrs = load_HDRs(args.hdr)
    # else:
    #     hdrs = None

    # if args.haplotype is not None:
    #     haplotypes = load_haplotypes(args.haplotype)
    # else:
    #     haplotypes = None

    # marker_data = {}
    # positions = {}
    # for method, marker, marker_gwa_fname, marker_gt_fname in args.marker:
    #     if method == 'loco':
    #         pos_col = 'BP'
    #     else:
    #         pos_col = 'POS'
    #     marker_gwa = read_gwas_data(
    #         filename=marker_gwa_fname,
    #         chrom_names=args.chrom_names,
    #         pos_col=pos_col
    #     )
    #     interval_markers = set([str(x) for x in data['marker']])
    #     marker_gt, _ = load_genotype_matrix(
    #         marker_gt_fname,
    #         markers=interval_markers,
    #         chrom_names=args.chrom_names)
    #     chrom = marker.split(':')[0]
    #     positions.setdefault(chrom, [])
    #     positions[chrom] += [x for x in marker_gwa['pos']]
    
    # for chrom in positions:
    #     positions[chrom] = set(positions[chrom])

    # annotations = {}
    # if args.annotation is not None:
    #     annotations['nonmt'] = load_annotation(args.annotation, positions=positions)
    # if mt_marker:
    #     annotations['mt'] = load_annotation(args.mt_annotation, positions=positions)


# def generate_report():
#     # Initialize report generator
#     report = HTMLReportGenerator(
#         title=f"NemaScan Report for {args.trait_name}",
#     )
    
#     if args.strain_issues is not None:
#         issue_text = "<br>".join(open(args.strain_issues).read().rstrip().split("\n"))
#         report.add_section(Section(
#             id="strain_issues",
#             title="Strain Issues",
#             level=1,
#             content=f"""
#             <p>All strain names were converted to the corresponding isotype name, which can be looked up here:
#             <a href "https://caendr.org/" class="uri">https://caendr.org/</a>
#             . If you submitted replicate data, replicates for a given isotype were averaged to one mean value.</p>
#             <div class="note-box">
#             {issue_text}
#             </div>
#             """
#         ))

#     loco_manhattan_fig = create_manhattan_plot(
#             loco_gwa,
#             title="LOCO",
#             bf_sig=loco_bf_sig,
#             eigen_sig=eigen_sig,
#             user_sig=args.user_sig,
#             point_size=6,
#             sig_colors=['#D41159', '#DC3220', '#D35FB7'],
#         )
#     inbred_manhattan_fig = create_manhattan_plot(
#             inbred_gwa,
#             title="Inbred",
#             bf_sig=inbred_bf_sig,
#             eigen_sig=eigen_sig,
#             user_sig=args.user_sig,
#             point_size=6,
#             sig_colors=['#D41159', '#DC3220', '#D35FB7'],
#         )
#     manhattan_section = report.add_section(Section(
#         id="manhattan_plots",
#         title="Manhattan Plots",
#         level=1,
#         content=f"""
#         <p>A genome-wide association study (GWAS) was performed by testing whether marker genotype differences can
#         explain phenotypic variation. These tests correct for relatedness among individuals in the population using
#         a genomic relatedness matrix (or “kinship matrix”). This anlaysis was performed with GCTA using two different
#         kinship matrices: one constructed specifically with inbred model organisms in mind (INBRED) and one which is
#         constructed from all markers except those on the chromosome of the tested marker (“leave-one-chromosome-out”;
#         LOCO). The INBRED kinship matrix more heavily corrects for genetic stratification at the tested marker, while
#         the LOCO kinship matrix does not, and may therefore increase power in certain contexts. <br>
#         <ul>
#         <li>Every dot is a SNV marker.</li>
#         <li>SNVs are colored if they pass the genome-wide corrected significance threshold:
#             <ul>
#             <li>The horizontal solid line corresponds to stricter Bonferroni (BF) threshold which is based on the
#             number of markers in the analysis.</li>
#             <li>The horizontal dash line corresponds to more permissive EIGEN threshold, which corrects for the number
#             of independent markers in your data set. This threshold takes advantage of the extensive LD in C. elegans
#             to limit the number of “unique” markers. (See Zdraljevic et al. 2019 (PMID: 30958264) for more)</li>
#             <li>If you selected a custom threshold, only this threshold is shown as a dotted line.</li>
#             </ul>
#         </li>
#         </ul>
#         </p>
#         """ + loco_manhattan_fig.to_html(full_html=False) + "<br>" + inbred_manhattan_fig.to_html(full_html=False) + "<br>" + \
#         report.create_table(headers=QTL_header, rows=QTLs, id='all_qtls')
#     ))


#     qq_fig, loco_lambda, inbred_lambda = create_qq_plot(
#         LOCO_data=loco_gwa,
#         inbred_data=inbred_gwa,
#         sig_cutoff=sig_cutoff,
#         title=None
#     )
#     qq_section = report.add_section(Section(
#         id="qq_plot",
#         title="Genomic Inflation",
#         level=1,
#         content=f"""
#         <p>The p-values calculated from each marker association test were compared to the theoretical distribution of p-values
#         under the null hypothesis. This comparison is displayed for each chromosome in the quantile-quantile plots
#         (Q-Q plots) below. The genomic inflation factor (λ_GC) estimates the inflation of observed p-values compared to a
#         theoretical χ^2 [0.5,1]. Mappings producing genomic inflation factors greater than 1.25 may indicate some systematic
#         bias, such as strong population stratification of phenotype values. <br>
#         """ + qq_fig.to_html(full_html=False) + "<br>" +\
#         f"""
#         <p><b>The genomic inflation factor is {inbred_lambda} for the INBRED mapping and {loco_lambda} for the LOCO mapping</b><br>
#         <i>The following sections of the report are shown for mappings performed using both the INBRED and LOCO kinship matrix
#         construction approaches. It is recommended you choose one set of results based on the previous diagnostic plots. These
#         results may vary between different traits.</i></p>
#         """
#     ))

#     qtl_datasets = {}
#     positions = {}
#     mt_marker = False
#     nonmt_marker = False
#     for method, marker, marker_gwa_fname, marker_gt_fname in args.marker:
#         qtl_datasets.setdefault(method, {})
#         qtl_datasets[method][marker] = []
#         qtl_datasets[method][marker].append(
#             read_gwas_data(
#                 filename=marker_gwa_fname,
#                 chrom_names=args.chrom_names,
#                 pos_col='POS'))
#         qtl_markers = set([str(x) for x in qtl_datasets[method][marker][0]['marker']])
#         qtl_datasets[method][marker].append(
#             load_genotype_matrix(
#                 marker_gt_fname,
#                 markers=qtl_markers,
#                 chrom_names=args.chrom_names)[0])
#         chrom = marker.split(':')[0]
#         if re.fullmatch("M|Mt|MtDNA|MT", chrom):
#             mt_marker = True
#         else:
#             nonmt_marker = True
#         positions.setdefault(chrom, set([]))
#         positions[chrom] = positions[chrom].union(set([x for x in qtl_datasets[method][marker][0]['pos']]))

#     annotations = {}
#     if args.annotation is not None and nonmt_marker:
#         annotations['nonmt'] = load_annotations(args.annotation, positions=positions)
#     else:
#         annotations['nonmt'] = None
#     if args.mt_annotation is not None and mt_marker:
#         annotations['mt'] = load_annotations(args.mt_annotation, positions=positions)
#     else:
#         annotations['mt'] = None

#     if args.strains is not None:
#         highlight_strains = args.strains.split(',')
#     else:
#         highlight_strains = None
    
#     for QTLs, method in [(inbred_QTLs, 'inbred'), (loco_QTLs, 'loco')]:
#         if method == 'inbred':
#             content = ["""
#                 <p>This is the default kinship matrix construction approach, designed for inbred model organisms (See 
#                 <a href="https://yanglab.westlake.edu.cn/software/gcta/#MakingaGRM">https://yanglab.westlake.edu.cn/software/gcta/#MakingaGRM</a>
#                 for more info).</p><br>
#                 """]
#         else:
#             content = ["""
#                 <p>LOCO may provide increased power to detect QTL because it does not correct for relatedness (or stratification) on the
#                 chromosome of each tested marker, sometimes providing higher power to detect linked QTL or QTL within divergent regions.
#                 However, this higher power also comes with higher false discovery rates. For more info, check out
#                 <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC9258552/">Widmayer et al. 2022.</a></p><br>
#                 """]
#         if len(QTLs) > 1:
#             content.append(report.create_table(headers=QTLs[0], rows=QTLs[1:], id=f'{method}_qtls'))
#         else:
#             content.append(f"""<div class="note-box">No significant QTLs were found using the {method} approach.</div>""")
#         report.add_section(Section(
#             id=method,
#             title=f'{method.upper()}',
#             level=1,
#             content="\n".join(content)
#         ))
#         if len(QTLs) > 1:
#             marker_names = [row[1] for row in QTLs[1:]]

#             pXg_fig = plot_phenotype_genotype(
#                 genotype_matrix,
#                 traits,
#                 strains,
#                 marker_names,
#                 highlight_strains=highlight_strains,
#                 title=None)
#             report.add_section(Section(
#                 id=f'{method}_pXg',
#                 title='Phenotype by Genotype Split',
#                 level=2,
#                 content="""
#                 <p>For each detected QTL, we can observe the phenotypes of the strains with the reference (REF) allele (i.e. same
#                 genotype as N2) compared to the phenotypes of the strains with the alternative (ALT) allele (i.e. genotype different
#                 than N2). A QTL is defined as a region where genetic variation is correlated with phenotypic variation, so we expect
#                 to see a difference in phenotype between the REF and ALT groups. In a best-case scenario, we like to see a large split
#                 between REF and ALT and a good number of strains in both groups. It is also important to ensure that the mean phenotype
#                 of neither group is driven by a small number of outlier strains.</p><br>""" + pXg_fig.to_html(full_html=False)
#             ))

#             content = """
#                 <p>If your trait has multiple QTL, we calculate linkage disequilibrium (LD) between them. This is useful because
#                 sometimes we find that one strong QTL might be in linkage disequilibrium to a secondary QTL (even if it exists on
#                 another chromosome). If this is the case, the secondary QTL might not contain a true causal variant, thus it is
#                 important to check this before narrowing the QTL experimentally.</p><br>"""
#             if len(QTLs) > 2:
#                 LD_matrix = find_LD(genotype_matrix, marker_names)
#                 fig = plot_LD_matrix(LD_matrix, marker_names, title=None)
#                 content += fig.to_html(full_html=False)
#             else:
#                 content += """<div class="note-box">Not enough QTL regions to calculate linkage disequilibrium.</div>"""
#             report.add_section(Section(
#                 id=f'{method}_LD',
#                 title='Linkage Disequilibrium',
#                 level=2,
#                 content=content,
#             ))

#             for marker in marker_names:
#                 report.add_section(Section(
#                     id=f'{marker}_{method}',
#                     title=f'{marker}',
#                     level=2,
#                     content=""
#                 ))

#                 content = ["""
#                 <p>Fine mapping was performed by evaluating the genotype-phenotype relationship for variants nearby the QTL
#                 identified from GWA mapping using a vcf containing imputed variants to avoid removing variants with missing
#                 genotype information for one or a few strains. Only SNVs were considered in this mapping.</p>"""]
#                 chrom = marker.split(':')[0]
#                 if re.fullmatch("M|Mt|MT|MtDNA", chrom):
#                     chrom_type = 'mt'
#                 else:
#                     chrom_type = 'nonmt'
#                 if (chrom_type == 'mt' and annotations['mt'] is not None) or (chrom_type == 'nonmt' and annotations['nonmt'] is not None):
#                     content.append("""
#                     <p>Each variant is represented by a vertical line, colored by the predicted variant impact (i.e. HIGH impact
#                     variants could be variants that introduce a change in the amino acid sequence or a stop-gain). Genes are
#                     represented by horizontal lines with an arrow showing the direction of the gene.</p><br>""")
#                     content.append(create_interval_variant_plot(
#                         qtl_datasets[method][marker][0],
#                         annotations[chrom_type][np.where(annotations[chrom_type]['chr'] == chrom)],
#                         genes,
#                         marker,
#                         point_size=8,
#                         title=None).to_html(full_html=False))
#                     content.append("""
#                     <div class="note-box">This second plot is very similar to the first. Here, each variant is represented by a
#                     circle colored by the linkage to the peak marker (colored in red). This plot can be useful to determine what
#                     the structure of your region looks like. If you have many variants with high linkage to your peak marker, it
#                     is important to remember that any of those variants could be causal.</div><br>""")
#                 else:
#                     content.append("""
#                     Each variant is represented by a circle colored by the linkage to the peak marker (colored in red). This plot
#                     can be useful to determine what the structure of your region looks like. If you have many variants with high
#                     linkage to your peak marker, it is important to remember that any of those variants could be causal.<br>""")
#                 content.append(create_interval_ld_plot(
#                     qtl_datasets[method][marker][0],
#                     qtl_datasets[method][marker][1],
#                     marker,
#                     point_size=8,
#                     title=None).to_html(full_html=False))

#                 report.add_section(Section(
#                     id=f'{marker}_{method}_finemapping',
#                     title=f'Fine Mapping',
#                     level=3,
#                     content="\n".join(content)
#                 ))

#                 if (chrom_type == 'mt' and annotations['mt'] is not None) or (chrom_type == 'nonmt' and annotations['nonmt'] is not None):
#                     header, variant_list = compile_variant_list(
#                         qtl_datasets[method][marker][0],
#                         qtl_datasets[method][marker][1],
#                         annotations[chrom_type],
#                         strains)
#                     content = report.create_table(headers=header, rows=variant_list, id=f'{method}_variant_list')

#                     report.add_section(Section(
#                         id=f'{marker}_{method}_variants',
#                         title=f'All variants in interval',
#                         level=3,
#                         content=content
#                     ))

#                 content = ["""
#                     <p>Mediation analysis was performed to analyze if gene expression variation is significantly correlated with the
#                     phenotype (overlap of phenotype QTL with expression QTL). Top candidates whose expression might mediate the phenotype
#                     QTL are shown below. (Note: expression data currently unpublished). For more information about mediation analysis,
#                     check out <a href="https://www.micropublication.org/journals/biology/micropub-biology-000305">Evans and Andersen 2020
#                     (PMID: 32385045)</a>.</p>
#                     """]
#                 if True:
#                     content.append("""
#                         <div class="note-box">Mediation was not performed for this analysis</div>
#                         """)
#                 report.add_section(Section(
#                     id=f'{marker}_{method}_mediation',
#                     title=f'Mediation analysis',
#                     level=3,
#                     content="\n".join(content)
#                 ))

#                 if hdrs is not None:
#                     content = ["""
#                         <p>We recently published about punctuated hyper-divergent regions in C. elegans (<a href="https://www.nature.com/articles/s41559-021-01435-x">
#                         Lee et al. 2021 (PMID: 32385045)</a>). Within these divergent regions, we are less confident about the variant calls and
#                         even the gene content between strains. For these reasons, if your QTL falls within a divergent region it may complicate
#                         your analyses and requires extra careful interpretation of fine-mapping results.<br>
#                         The following plot shows divergent regions for each strain across the QTL region. Strains are split by genotype at the peak
#                         marker. You should be careful if many strains are divergent, especially if most of the strains in the ALT group are divergent,
#                         for example.</p>
#                         """]
#                     content.append(create_HDR_plot(
#                         qtl_datasets[method][marker][1]['genotype'][np.where(qtl_datasets[method][marker][1]['marker'] == marker)[0][0], :],
#                         HDRs=hdrs,
#                         strains=strains,
#                         marker=marker,
#                         start=qtl_datasets[method][marker][0]['pos'][0],
#                         stop=qtl_datasets[method][marker][0]['pos'][-1],
#                         title=''
#                     ).to_html(full_html=False, default_height=8*(len(strains) + 3) + 80))
#                     report.add_section(Section(
#                         id=f'{marker}_{method}_hdr',
#                         title=f'Hyper-divergent Regions',
#                         level=3,
#                         content="\n".join(content)
#                     ))
                
#                 if haplotypes is not None:
#                     content = ["""
#                         <p>The following plot shows the genome-wide haplotype (genetic relatedness) of mapped strains split by REF or ALT genotype.
#                         This plot can be useful to help identify how many unique haplotypes are present in the REF or ALT groups. If you want to
#                         choose parent strains for a NIL cross to validate this QTL, you might want to choose strains in the major haplotype of the
#                         REF/ALT groups that also have distinct phenotypes.</p>
#                         """]
#                     content.append(create_haplotype_plot(
#                         qtl_datasets[method][marker][1]['genotype'][np.where(qtl_datasets[method][marker][1]['marker'] == marker)[0][0], :],
#                         haplotypes=haplotypes,
#                         strains=strains,
#                         marker=marker,
#                         start=qtl_datasets[method][marker][0]['pos'][0],
#                         stop=qtl_datasets[method][marker][0]['pos'][-1],
#                         title=''
#                     ).to_html(full_html=False, default_height=8*(len(strains) + 3) + 40))
#                     report.add_section(Section(
#                         id=f'{marker}_{method}_haplotype',
#                         title=f'Haplotypes',
#                         level=3,
#                         content="\n".join(content)
#                     ))
#     report.save(args.output)

def make_arg_parser():
    parser = argparse.ArgumentParser(
        description='Create an interactive Manhattan plot from GWAS results'
    )
     
    parser.add_argument(
        '--output',
        required=True,
        help='Output report filename'
    )
    
    parser.add_argument(
        '--gwas',
        required=True,
        action='append',
        nargs=3, # method, gwa, sigQTL
        help='GWAS configuration file'
    )
    
    parser.add_argument(
        '--finemapping',
        action='append',
        nargs=3, # method, marker, annotated_gwa
        help='GWAS configuration file'
    )
    
    parser.add_argument(
        '--mediation',
        action='append',
        nargs=4, # method, marker, multimed, singlemed
        help='GWAS configuration file'
    )
    
    parser.add_argument(
        '-g', '--genotype-matrix',
        required=True,
        help='Genotype matrix file'
    )
    
    parser.add_argument(
        '-t', '--trait',
        required=True,
        help='Trait name'
    )
    
    parser.add_argument(
        '-p', '--phenotype',
        required=True,
        help='Trait phenotype file'
    )
    
    parser.add_argument(
        '--hdr',
        help='Input HDR file'
    )
    
    parser.add_argument(
        '--haplotype',
        help='Input haplotype file'
    )
    
    parser.add_argument(
        '--strain-issues',
        default=None,
        help='File with strain issues'
    )
    
    parser.add_argument(
        '--strains',
        default=None,
        help='Comma-separated list of strains to highlight'
    )
    
    parser.add_argument(
        '--alpha',
        type=float,
        default=0.05,
        help='Significance threshold (default: 0.05)'
    )
    
    parser.add_argument(
        '--independent-tests',
        type=float,
        default=None,
        help='Eigen number of independent tests'
    )
    
    parser.add_argument(
        '--significance',
        default="BF",
        help='Significance method (BF or EIGEN), or a user-specified significance threshold (default: BF)'
    )
    
    return parser

def read_gwas_data(filename, pos_col='POS', extended=False):
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
    data = []
    with open(filename, 'r') as f:
        header = [x.upper() for x in f.readline().rstrip().split("\t")]
        chrom_col = header.index("CHR")
        pos_col = header.index(pos_col)
        pval_col = header.index('P')
        if extended:
            ld_col = header.index('LD')
            wbgene_col = header.index('WBGENE')
            gene_col = header.index('GENE_NAME')
            consequence_col = header.index('CONSEQUENCE')
            aa_col = header.index('AA')
        
        for line in f:
            row = line.rstrip().split("\t")
            pval = float(row[pval_col])
            # Skip invalid p-values
            if pval <= 0 or pval > 1:
                continue
            
            chrom = row[chrom_col].replace('chr', '').replace('Chr', '')
            pos = int(row[pos_col])
            snp = f"{chrom}:{pos}"
            ref = row[pos_col + 1]
            alt = row[pos_col + 2]
            if extended:
                ld = float(row[ld_col])
                wbgene = row[wbgene_col]            
                gene = row[gene_col]            
                consequence = row[consequence_col]
                aa = row[aa_col]
                data.append((chrom, pos, ref, alt, snp, pval, -math.log10(pval), ld, wbgene, gene, consequence, aa))
            else:
                data.append((chrom, pos, ref, alt, snp, pval, -math.log10(pval)))

    dtype = [
        ('chrom', f'U{min([len(row[0]) for row in data])}'),
        ('pos', np.int32),
        ('ref', 'U1'),
        ('alt', 'U1'),
        ('marker', f'U{min([len(row[4]) for row in data])}'),
        ('pval', np.float32),
        ('log_pval', np.float32)]
    if extended:
        dtype += [
            ('ld', np.float32),
            ('wbgene', f'U{min([len(row[8]) for row in data])}'),
            ('gene', f'U{min([len(row[9]) for row in data])}'),
            ('consequence', f'U{min([len(row[10]) for row in data])}'),
            ('aa', f'U{min([len(row[11]) for row in data])}')]
    data = np.array(data, dtype=dtype)
    data = data[np.lexsort((data['pos'], data['chrom']))]
    return data

# def compile_variant_list(gwa, gt, variants, strains):
#     valid = np.where(np.logical_and(np.logical_and(variants['chr'] == gwa['chrom_label'][0], variants['pos'] >= gwa['pos'][0]),
#                                     variants['pos'] <= gwa['pos'][-1]))[0]
#     gwa_indices = np.searchsorted(gwa['pos'], variants['pos'][valid])
#     gt_indices = np.searchsorted(gt['pos'], variants['pos'][valid])
#     variant_list = []
#     for i, v in enumerate(valid):
#         j = gwa_indices[i]
#         k = gt_indices[i]
#         if variants['pos'][v] != gwa['pos'][j] or variants['pos'][v] != gt['pos'][k]:
#             continue
#         alt = np.where(gt['genotype'][k, :] == 0)[0]
#         variant_list.append([
#             str(gwa['marker'][j]),
#             str(gwa['ref'][j]),
#             str(gwa['alt'][j]),
#             str(variants['WBGENE'][v]),
#             str(variants['GENE_NAME'][v]),
#             str(variants['CONSEQUENCE'][v]),
#             f"{gwa['log_pval'][j]:0.4f}",
#             ','.join([strains[x] for x in alt])])
#     return ['MARKER', 'REF', 'ALT', 'WBGeneID', 'GENE NAME', 'CONSEQUENCE', 'VARIANT LOG10p', 'ALT STRAINS'], variant_list



# @dataclass
# class Section:
#     """Represents a section in the report"""

#     def __init__(self, id: str, title: str = "section", level: int = 1, content: str = "", subsections: list = []):
#         self.id = id
#         self.title = title
#         self.level = level
#         self.content = content
#         self.subsections = subsections

#     def add_subsection(self, subsection):
#         """Add a section to the report"""
#         self.subsections.append(subsection)
    

# class HTMLReportGenerator:
#     """Main class for generating HTML reports"""
    
#     def __init__(self, title: str = "Report"):
#         self.title = title
#         self.sections: List[Section] = []
#         self.custom_css = ""
#         self.body_script = ["let tables = {};"]
        
#     def _get_css(self) -> str:
#         """Return the CSS styles for the report"""
#         return """
#         :root {
#             --primary-color: #2c3e50;
#             --secondary-color: #3498db;
#             --accent-color: #e74c3c;
#             --background-color: #f8f9fa;
#             --text-color: #333;
#             --border-color: #dee2e6;
#         }
        
#         * {
#             box-sizing: border-box;
#             margin: 0;
#             padding: 0;
#         }
        
#         body {
#             font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
#             line-height: 1.6;
#             color: var(--text-color);
#             background-color: var(--background-color);
#         }
        
#         .container {
#             max-width: 1200px;
#             margin: 0 auto;
#             padding: 20px;
#         }
        
#         /* Header Styles */
#         .report-header {
#             background: linear-gradient(135deg, var(--primary-color), var(--secondary-color));
#             color: white;
#             padding: 40px;
#             text-align: center;
#             margin-bottom: 30px;
#             border-radius: 8px;
#             box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
#         }
        
#         .report-header h1 {
#             font-size: 2.5em;
#             margin-bottom: 10px;
#         }
        
#         .report-header .meta {
#             font-size: 0.9em;
#             opacity: 0.9;
#         }
        
#         .report-header .logo {
#             max-height: 80px;
#             margin-bottom: 20px;
#         }
        
#         /* Table of Contents */
#         .toc {
#             background: white;
#             padding: 30px;
#             border-radius: 8px;
#             box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
#             margin-bottom: 0px;
#         }
        
#         .toc h2 {
#             color: var(--primary-color);
#             border-bottom: 3px solid var(--secondary-color);
#             padding-bottom: 10px;
#             margin-bottom: 20px;
#         }
        
#         .toc ul {
#             list-style: none;
#         }
        
#         .toc li {
#             padding: 5px 0;
#         }
        
#         .toc a {
#             color: var(--text-color);
#             text-decoration: none;
#             transition: color 0.3s;
#         }
        
#         .toc a:hover {
#             color: var(--secondary-color);
#         }
        
#         .toc .toc-level-1 { margin-left: 0; font-weight: bold; }
#         .toc .toc-level-2 { margin-left: 20px; }
#         .toc .toc-level-3 { margin-left: 40px; font-size: 0.9em; }
        
#         /* Section Styles */
#         .section {
#             background: white;
#             padding-left: 30px;
#             padding-right: 30px;
#             padding-bottom: 10px;
#             border-radius: 8px;
#             box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
#             margin-bottom: 10px;
#         }
        
#         .section h1, .section h2, .section h3 {
#             color: var(--primary-color);
#             margin-bottom: 20px;
#         }
        
#         .section h1 {
#             font-size: 2em;
#             border-bottom: 3px solid var(--secondary-color);
#             padding-bottom: 10px;
#         }
        
#         .section h2 {
#             font-size: 1.5em;
#             border-bottom: 2px solid var(--border-color);
#             padding-bottom: 8px;
#         }
        
#         .section h3 {
#             font-size: 1.2em;
#         }
        
#         .section p {
#             margin-bottom: 15px;
#             text-align: justify;
#         }
        
#         .section ul {
#             padding-inline-start: 40px;
#         }
        
#         /* Table Styles */
#         .data-table {
#             width: 100%;
#             border-collapse: collapse;
#             margin: 20px 0;
#             font-size: 0.9em;
#         }
        
#         .data-table thead {
#             background: var(--primary-color);
#             color: white;
#         }
        
#         .data-table th, .data-table td {
#             padding: 12px 15px;
#             text-align: left;
#             border-bottom: 1px solid var(--border-color);
#         }
        
#         .data-table tbody tr:hover {
#             background-color: #f5f5f5;
#         }
        
#         .data-table tbody tr:nth-child(even) {
#             background-color: #fafafa;
#         }
        
#         /* Chart Styles */
#         .chart-container {
#             text-align: center;
#             margin: 20px 0;
#         }
        
#         .chart-container img {
#             max-width: 100%;
#             height: auto;
#             border-radius: 4px;
#         }
        
#         .chart-caption {
#             font-style: italic;
#             color: #666;
#             margin-top: 10px;
#             font-size: 0.9em;
#         }
        
#         /* Alert/Info Boxes */
#         .alert {
#             padding: 15px 20px;
#             border-radius: 4px;
#             margin: 15px 0;
#         }
        
#         .alert-info {
#             background-color: #d1ecf1;
#             border-left: 4px solid #17a2b8;
#             color: #0c5460;
#         }
        
#         .alert-warning {
#             background-color: #fff3cd;
#             border-left: 4px solid #ffc107;
#             color: #856404;
#         }
        
#         .alert-success {
#             background-color: #d4edda;
#             border-left: 4px solid #28a745;
#             color: #155724;
#         }
        
#         .alert-danger {
#             background-color: #f8d7da;
#             border-left: 4px solid #dc3545;
#             color: #721c24;
#         }
        
#         /* Code Block */
#         .code-block {
#             background-color: #2d2d2d;
#             color: #f8f8f2;
#             padding: 15px;
#             border-radius: 4px;
#             overflow-x: auto;
#             font-family: 'Consolas', 'Monaco', monospace;
#             font-size: 0.9em;
#             margin: 15px 0;
#         }
        
#         /* Footer */
#         .report-footer {
#             text-align: center;
#             padding: 20px;
#             color: #666;
#             font-size: 0.9em;
#             border-top: 1px solid var(--border-color);
#             margin-top: 30px;
#         }
        
#         /* Print Styles */
#         @media print {
#             body {
#                 background-color: white;
#             }
            
#             .section {
#                 box-shadow: none;
#                 break-inside: avoid;
#             }
            
#             .toc {
#                 break-after: page;
#             }
#         }
        
#         /* Responsive */
#         @media (max-width: 768px) {
#             .container {
#                 padding: 10px;
#             }
            
#             .report-header {
#                 padding: 20px;
#             }
            
#             .report-header h1 {
#                 font-size: 1.8em;
#             }
            
#             .section {
#                 padding: 20px;
#             }
#         }
#         .note-box {
#             padding: 15px; /* Adds space inside the box */
#             margin: 20px 0; /* Adds space above and below the box */
#             background-color: #ccdded; /* A light, soft blue background */
#             border-left: 5px solid #007bff; /* A prominent left border for emphasis */
#             border-radius: 4px; /* Slightly rounded corners */
#             font-family: sans-serif; /* Readable font */
#         }
#         """ + self.custom_css
    
#     def add_custom_css(self, css: str):
#         """Add custom CSS to the report"""
#         self.custom_css += css
    
#     def add_section(self, section: Section):
#         """Add a section to the report"""
#         self.sections.append(section)
#         return section
    
#     def create_section(self, id: str, title: str, content: str = "", 
#                        level: int = 1) -> Section:
#         """Create and add a new section"""
#         section = Section(id=id, title=title, content=content, level=level)
#         self.sections.append(section)
#         return section
    
#     def create_table(self, headers: List[str], rows: List[List[Any]], 
#                      caption: str = "", id: str = "table") -> str:
#         """Create an HTML table"""
#         html = f'<table class="display" id="{id}">'
        
#         if caption:
#             html += f'<caption>{caption}</caption>'
        
#         # Header
#         html += '<thead><tr>'
#         for header in headers:
#             html += f'<th>{header}</th>'
#         html += '</tr></thead>'
        
#         # Body
#         html += '<tbody>'
#         for row in rows:
#             html += '<tr>'
#             for cell in row:
#                 html += f'<td>{cell}</td>'
#             html += '</tr>'
#         html += '</tbody></table>'
        
#         self.body_script.append(f"tables[{id}] = new DataTable('#{id}');")
#         return html
    
#     @staticmethod
#     def create_chart_html(chart_base64: str, caption: str = "") -> str:
#         """Create HTML for embedding a chart"""

#         html = f'''
#         <div class="chart-container">
#             <img src="{chart_base64}" alt="Chart">
#             {f'<p class="chart-caption">{caption}</p>' if caption else ''}
#         </div>
#         '''
#         return html
    
#     @staticmethod
#     def create_alert(message: str, alert_type: str = "info") -> str:
#         """Create an alert/info box"""
#         return f'<div class="alert alert-{alert_type}">{message}</div>'
    
#     @staticmethod
#     def create_code_block(code: str) -> str:
#         """Create a code block"""
#         return f'<pre class="code-block">{code}</pre>'
    
#     def _generate_toc(self) -> str:
#         """Generate table of contents HTML"""
#         html = '<nav class="toc"><h2>Table of Contents</h2><ul>'
        
#         def process_sections(sections: List[Section]) -> str:
#             toc_html = ""
#             for section in sections:
#                 toc_html += f'''
#                 <li class="toc-level-{section.level}">
#                     <a href="#{section.id}">{section.title}</a>
#                 </li>
#                 '''
#                 if section.subsections:
#                     toc_html += f'<ul>{process_sections(section.subsections)}</ul>'
#             return toc_html
        
#         html += process_sections(self.sections)
#         html += '</ul></nav>'
#         return html
    
#     def _generate_sections_html(self) -> str:
#         """Generate HTML for all sections"""
#         def render_section(section: Section) -> str:
#             tag = f'h{min(section.level, 6)}'
#             html = f'''
#             <div class="section" id="{section.id}">
#                 <{tag}>{section.title}</{tag}>
#                 {section.content}
#             '''
            
#             for subsection in section.subsections:
#                 html += render_section(subsection)
            
#             html += '</div>'
#             return html
        
#         return ''.join(render_section(s) for s in self.sections)
    
#     def generate(self, include_toc: bool = True) -> str:
#         """Generate the complete HTML report"""
#         generated_date = datetime.now().strftime("%B %d, %Y at %H:%M")
#         body_script = "\n".join(self.body_script)

#         html = f'''<!DOCTYPE html>
# <html lang="en">
# <head>
#     <meta charset="UTF-8">
#     <meta name="viewport" content="width=device-width, initial-scale=1.0">
#     <meta name="generated" content="{generated_date}">
#     <title>{self.title}</title>
#     <style>
#     {self._get_css()}
#     </style>
#     <script src="https://code.jquery.com/jquery-4.0.0.slim.min.js"></script>
#     <link href="https://cdn.datatables.net/2.3.7/css/dataTables.dataTables.min.css" rel="stylesheet" crossorigin="anonymous"> 
#     <script src="https://cdn.datatables.net/2.3.7/js/dataTables.min.js" crossorigin="anonymous"></script>
#     </head>
# <body>
#     <div class="container">
#         <header class="report-header">
#             <h1>{self.title}</h1>
#             <p class="meta">
#                 Generated: {generated_date}
#             </p>
#         </header>
        
#         {self._generate_toc() if include_toc else ''}
        
#         <main>
#             {self._generate_sections_html()}
#         </main>
        
#         <footer class="report-footer">
#             <p>Report generated on {generated_date}</p>
#         </footer>
#     </div>
# <script>
# {body_script}
# </script>
# </body>
# </html>'''
        
#         return html
    
#     def save(self, filepath: str, include_toc: bool = True):
#         """Save the report to a file"""
#         html = self.generate(include_toc)
#         with open(filepath, 'w', encoding='utf-8') as f:
#             f.write(html)
#         print(f"Report saved to: {filepath}")


# ============================================================================
# EXAMPLE USAGE
# ============================================================================

def create_sample_report():
    """Create a sample report demonstrating all features"""
    
    # Initialize report generator
    report = HTMLReportGenerator(
        title="Monthly Sales Analysis Report",
        author="Data Analytics Team"
    )
    
    # Section 1: Executive Summary
    section1 = report.create_section(
        id="executive-summary",
        title="Executive Summary",
        level=1,
        content="""
        <p>This report provides a comprehensive analysis of the monthly sales performance
        across all regions. The data covers the period from January to December 2024.</p>
        
        <p>Key findings indicate a <strong>15% increase</strong> in overall sales compared
        to the previous year, with the Western region showing the most significant growth.</p>
        """ + report.create_alert(
            "Highlight: Q4 sales exceeded targets by 23%!", 
            "success"
        )
    )
    
    # Section 2: Sales Overview with Charts
    sales_data = {
        'January': 45000,
        'February': 52000,
        'March': 48000,
        'April': 61000,
        'May': 55000,
        'June': 67000
    }
    
    bar_chart = ChartGenerator.bar_chart(
        data=sales_data,
        title="Monthly Sales (First Half 2024)",
        xlabel="Month",
        ylabel="Sales ($)",
        color="#3498db"
    )
    
    section2 = report.create_section(
        id="sales-overview",
        title="Sales Overview",
        level=1,
        content="""
        <p>The following chart illustrates the monthly sales trends for the first half of 2024.</p>
        """ + report.create_chart_html(bar_chart, "Figure 1: Monthly Sales Performance")
    )
    
    # Section 3: Regional Analysis with Pie Chart
    regional_data = {
        'North': 125000,
        'South': 98000,
        'East': 110000,
        'West': 145000
    }
    
    pie_chart = ChartGenerator.pie_chart(
        data=regional_data,
        title="Sales Distribution by Region"
    )
    
    section3 = report.create_section(
        id="regional-analysis",
        title="Regional Analysis",
        level=1,
        content="""
        <p>Sales distribution across regions shows the West leading in total revenue.</p>
        """ + report.create_chart_html(pie_chart, "Figure 2: Regional Sales Distribution")
    )
    
    # Section 4: Detailed Data Table
    table_headers = ["Region", "Q1 Sales", "Q2 Sales", "Q3 Sales", "Q4 Sales", "Total"]
    table_rows = [
        ["North", "$28,500", "$31,200", "$35,800", "$29,500", "$125,000"],
        ["South", "$22,000", "$25,500", "$24,000", "$26,500", "$98,000"],
        ["East", "$26,000", "$28,000", "$27,500", "$28,500", "$110,000"],
        ["West", "$32,000", "$36,500", "$38,000", "$38,500", "$145,000"],
        ["<strong>Total</strong>", "<strong>$108,500</strong>", "<strong>$121,200</strong>", 
         "<strong>$125,300</strong>", "<strong>$123,000</strong>", "<strong>$478,000</strong>"]
    ]
    
    section4 = report.create_section(
        id="detailed-data",
        title="Detailed Sales Data",
        level=1,
        content="""
        <p>The table below provides a quarterly breakdown of sales by region.</p>
        """ + report.create_table(
            headers=table_headers,
            rows=table_rows,
            caption="Table 1: Quarterly Sales by Region"
        ) + report.create_alert(
            "Note: All figures are in USD and represent gross sales before returns.",
            "info"
        )
    )
    
    # Section 5: Trend Analysis with Line Chart
    trend_data = {
        '2023': [40000, 42000, 45000, 48000, 50000, 55000],
        '2024': [45000, 52000, 48000, 61000, 55000, 67000]
    }
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun']
    
    line_chart = ChartGenerator.line_chart(
        data=trend_data,
        x_labels=months,
        title="Year-over-Year Sales Comparison",
        xlabel="Month",
        ylabel="Sales ($)"
    )
    
    section5 = report.create_section(
        id="trend-analysis",
        title="Trend Analysis",
        level=1,
        content="""
        <p>Comparing 2023 and 2024 performance shows consistent growth across all months.</p>
        """ + report.create_chart_html(line_chart, "Figure 3: Year-over-Year Comparison")
    )
    
    # Section 6: Performance Metrics with Horizontal Bar
    performance_data = {
        'Customer Satisfaction': 4.5,
        'Delivery Speed': 4.2,
        'Product Quality': 4.7,
        'Value for Money': 4.0,
        'Support Response': 3.8
    }
    
    hbar_chart = ChartGenerator.horizontal_bar_chart(
        data=performance_data,
        title="Customer Satisfaction Metrics (out of 5)",
        xlabel="Rating",
        ylabel="Category"
    )
    
    section6 = report.create_section(
        id="performance-metrics",
        title="Performance Metrics",
        level=1,
        content="""
        <p>Customer feedback metrics show strong performance across all categories.</p>
        """ + report.create_chart_html(hbar_chart, "Figure 4: Customer Satisfaction Ratings")
    )
    
    # Section 7: Conclusions
    section7 = report.create_section(
        id="conclusions",
        title="Conclusions and Recommendations",
        level=1,
        content="""
        <p>Based on our analysis, we recommend the following actions:</p>
        <ul>
            <li>Increase marketing spend in the South region to boost underperforming sales</li>
            <li>Maintain current strategies in the West region</li>
            <li>Focus on improving support response times based on customer feedback</li>
            <li>Plan for continued growth in Q1 2025</li>
        </ul>
        """ + report.create_alert(
            "Action Required: Schedule strategy meeting for Q1 2025 planning.",
            "warning"
        ) + """
        <h3>Technical Notes</h3>
        <p>Data processing methodology:</p>
        """ + report.create_code_block("""
# Data aggregation script
import pandas as pd

sales_df = pd.read_csv('sales_data.csv')
quarterly_sales = sales_df.groupby(['region', 'quarter']).sum()
        """)
    )
    
    # Generate and save the report
    report.save("sample_report.html")
    
    return report


if __name__ == "__main__":
    main()
