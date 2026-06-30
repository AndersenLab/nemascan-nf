nextflow.enable.types = true

process LD_PRUNED_MARKERS {
    label 'process_medium_progressive'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/plink:1.90b6.21--h779adbc_1"

    input:
    record(
        vcf: Path,
        vcf_index: Path
    )
    maf: BigDecimal

    output:
    file("plink.prune.in")

    topic:
    record(tool:"plink", version:eval("plink --version 2>&1 | sed 's/^PLINK v//' | sed 's/..-bit.*//'")) >> 'versions'

    script:
    """
    plink \\
        --vcf ${vcf} \\
        --threads ${task.cpus} \\
        --snps-only \\
        --biallelic-only \\
        --maf ${maf} \\
        --set-missing-var-ids @:# \\
        --indep-pairwise 50 10 0.8 \\
        --geno \\
        --allow-extra-chr
    """

    stub:
    """
    touch plink.prune.in
    """
}
