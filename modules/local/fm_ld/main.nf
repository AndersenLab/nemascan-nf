nextflow.enable.types = true

process FM_LD {
    tag "${meta.id}"
    label 'process_medium_progressive'

    conda "${moduleDir}/environment.yml"
    container "docker://quay.io/biocontainers/plink:1.90b6.21--h779adbc_1"

    input:
    record(
        meta: Map,
        vcf: Path,
        vcf_index: Path,
        chrom_numbers: Path
    )

    output:
    record(
        meta: meta,
        ld: file("${meta.id}.LD")
    )

    topic:
    record(tool:"plink", version:eval("plink --version 2>&1 | sed 's/^PLINK v//' | sed 's/..-bit.*//'")) >> 'versions'

    script:
    """
    nsnps=`zcat ${vcf} | wc -l`
    chrom_num=`grep -w ${meta.chrom} ${chrom_numbers} | cut -f 2`

    plink --r2 with-freqs \\
        --threads ${task.cpus} \\
        --allow-extra-chr \\
        --snps-only \\
        --ld-window-r2 0 \\
        --ld-snp \$chrom_num:${meta.peak} \\
        --ld-window \$nsnps \\
        --ld-window-kb 6000 \\
        --chr \$chrom_num \\
        --out ${meta.id}.QTL \\
        --set-missing-var-ids @:# \\
        --vcf ${vcf}
    tail -n +2 ${meta.id}.QTL.ld | \\
        sed 's/^[[:space:]]*//;s/[[:space:]]*\$//' | \\
        tr -s ' ' | \\
        sed 's/ /\\t/g' | \\
        sort -k6,6n -k7,7 | \\
        cut -f 7,9 > ${meta.id}.LD
    """

    stub:
    """
    touch ${meta.id}.LD
    """
}
