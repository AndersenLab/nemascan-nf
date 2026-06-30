nextflow.enable.types = true

process VCF_TO_PLINK {
    label 'process_medium_small'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/plink:1.90b6.21--h779adbc_1"

    input:
    record(
        vcf: Path,
        vcf_index: Path
    )

    output:
    record(
        meta: [plink_prefix: "gwa"],
        bed: file("gwa.bed"),
        bim: file("gwa.bim"),
        fam: file("gwa.fam")
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
        --biallelic-only \\
        --set-missing-var-ids @:# \\
        --geno \\
        --recode \\
        --out gwa \\
        --allow-extra-chr
    """

    stub:
    """
    touch gwa.bed
    touch gwa.bim
    touch gwa.fam
    """
}
