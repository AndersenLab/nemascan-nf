nextflow.enable.types = true

process VCF_FM_FILTER {
    tag "${meta.id}"
    label 'process_medium_small'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f"

    input:
    record(
        meta: Map,
        vcf: Path,
        vcf_index: Path,
        strains: Path,
        chrom_numbers: Path
    )

    output:
    record(
        meta: meta,
        vcf: file("${meta.id}.vcf.gz"),
        vcf_index: file("${meta.id}.vcf.gz.tbi")
    )

    topic:
    record(tool:"bcftools", version:eval("bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//'")) >> 'versions'

    script:
    def mt_option = meta.chrom == 'MtDNA' ? "& GT!=\"het\"" : ""
    """
    echo -e "${meta.chrom}\\t${meta.start}\\t${meta.end}" > roi.txt

    bcftools view \\
        --output-type v \\
        --samples-file ${strains} \\
        --threads ${task.cpus} \\
        --regions-file roi.txt \\
        ${vcf} | \\
    bcftools filter \\
        --include 'N_MISSING=0 & TYPE="snp" ${mt_option}' \\
        --threads ${task.cpus} | \\
    bcftools annotate \\
        --rename-chrs ${chrom_numbers} \\
        --output-type z \\
        --output ${meta.id}.vcf.gz
        
    bcftools index -t ${meta.id}.vcf.gz
    """

    stub:
    """
    touch ${meta.id}.vcf.gz
    touch ${meta.id}.vcf.gz.tbi
    """
}
