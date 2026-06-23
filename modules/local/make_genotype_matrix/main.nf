nextflow.enable.types = true

process MAKE_GENOTYPE_MATRIX {
    label 'process_medium_small_progressive'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f"

    input:
    record(
        vcf: Path,
        vcf_index: Path,
        markers: Path
    )

    output:
    record(
        matrix: file("Genotype_Matrix.tsv"),
        vcf: record(
            vcf: file("pruned.vcf.gz"),
            vcf_index: file("pruned.vcf.gz.tbi")
        )
    )

    topic:
    record(tool:"bcftools", version:eval("bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//'")) >> 'versions'

    script:
    """
    awk -F":" '\$1=\$1' OFS="\\t" ${markers} | \\
    sort -k1,1d -k2,2n > markers.txt

    bcftools view \\
        --types snps \\
        --regions-file markers.txt \\
        --output-type z \\
        --output pruned.vcf.gz \\
        ${vcf}
    bcftools index -t pruned.vcf.gz

    bcftools query \\
        --print-header \\
        --format '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' \\
        pruned.vcf.gz | \\
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
    sed 's/.\\/./NA/g' > Genotype_Matrix.tsv
    """

    stub:
    """
    touch Genotype_Matrix.tsv
    touch pruned.vcf.gz
    touch pruned.vcf.gz.tbi
    """
}
