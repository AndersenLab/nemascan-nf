nextflow.enable.types = true

process MAKE_FM_GENOTYPE_MATRIX {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f"

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
        matrix: file("${meta.id}_genotype_matrix.tsv")
    )

    topic:
    record(tool:"bcftools", version:eval("bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//'")) >> 'versions'

    script:
    """
    bcftools query \\
        --print-header \\
        --format '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' \\
        ${vcf} | \\
    sed 's/[[# 0-9]*]//g' | \\
    sed 's/:GT//g' | \\
    sed 's/0|0/-1/g' | \\
    sed 's/1|1/1/g' | \\
    sed 's/0|1/NA/g' | \\
    sed 's/1|0/NA/g' | \\
    sed 's/.|./NA/g'  | \\
    sed 's/0\\/0/-1/g' | \\
    sed 's/1\\/1/1/g'  | \\
    sed 's/0\\/1/NA/g' | \\
    sed 's/1\\/0/NA/g' | \\
    sed 's/.\\/./NA/g' | \\
    awk '{
        if (NR == FNR) CHROMS[\$2] = \$1;
        else \$1 = CHROMS[\$1];
    }' ${chrom_numbers} - > ${meta.id}_genotype_matrix.tsv
    """

    stub:
    """
    touch ${meta.id}_genotype_matrix.tsv
    """
}
