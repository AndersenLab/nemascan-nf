nextflow.enable.types = true

process VCF_FILTER {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f"

    input:
    record(
        vcf: Path,
        vcf_index: Path,
        samples: Path
    )

    output:
    record(
        vcf: file("filtered.vcf.gz"),
        vcf_index: file("filtered.vcf.gz.tbi")
    )

    topic:
    record(tool:"bcftools", version:eval("bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//'")) >> 'versions'

    script:
    """
    bcftools view \\
        --output-type v \\
        --samples-file ${samples} \\
        --threads ${task.cpus} \\
        ${vcf} | \\
    bcftools filter \\
        --include 'N_MISSING=0 & TYPE="snp"' \\
        --output filtered.vcf.gz \\
        --threads ${task.cpus} \\
        --output-type z
        
    bcftools index -t filtered.vcf.gz
    """

    stub:
    """
    touch filtered.vcf.gz
    touch filtered.vcf.gz.tbi
    """
}
