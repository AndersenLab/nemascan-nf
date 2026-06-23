nextflow.enable.types = true

process VCF_TO_FM_PLINK {
    tag "${meta.id}"
    label 'process_medium_small'

    conda "${moduleDir}/environment.yml"
    container "docker://quay.io/biocontainers/plink:1.90b6.21--h779adbc_1"

    input:
    record(
        meta: Map,
        vcf: Path,
        vcf_index: Path
    )
    maf: BigDecimal

    output:
    record(
        meta: meta + [plink_prefix: "${meta.id}"],
        bed: file("${meta.id}.bed"),
        bim: file("${meta.id}.bim"),
        fam: file("${meta.id}.fam")
    )

    topic:
    record(tool:"plink", version:eval("plink --version 2>&1 | sed 's/^PLINK v//' | sed 's/..-bit.*//'")) >> 'versions'

    script:
    """
    plink \\
        --vcf ${vcf} \\
        --threads ${task.cpus} \\
        --make-bed \\
        --snps-only \\
        --maf ${maf} \\
        --biallelic-only \\
        --set-missing-var-ids @:# \\
        --geno \\
        --out ${meta.id} \\
        --allow-extra-chr
    """

    stub:
    """
    touch ${meta.id}.bed
    touch ${meta.id}.bim
    touch ${meta.id}.fam
    """
}


