nextflow.enable.types = true

process ANNOTATE_CHROMS {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f"

    input:
    record(
        vcf: Path,
        vcf_index: Path,
        chrom_numbering: Path
    )

    output:
    record(
        vcf: file("annotated.vcf.gz"),
        vcf_index: file("annotated.vcf.gz.tbi")
    )

    topic:
    record(tool:"bcftools", version:eval("bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//'")) >> 'versions'

    script:
    """
    tail -n +2 ${chrom_numbering} | cut -f 1,2 > chrom_numbering_noheader.txt
    bcftools annotate \\
        --rename-chrs chrom_numbering_noheader.txt \\
        --output-type z \\
        --threads ${task.cpus} \\
        --output annotated.vcf.gz \\
        ${vcf}
    bcftools index -t annotated.vcf.gz
    """

    stub:
    """
    touch annotated.vcf.gz
    touch annotated.vcf.gz.tbi
    """
}
