nextflow.enable.types = true

process SINGLE_MEDIATION {
    tag "${meta.id}"
    label "single_mediation"
    maxRetries 0

    conda "${moduleDir}/environment.yml"
    container "docker://andersenlab/mediation:20220407173046db3227"

    input:
    record (
        meta: Map,
        trait: Path,
        eqtl: Path,
        expression: Path,
        genotype_matrix: Path,
        medmulti: Path,
        genes: Path
    )

    output:
    record (
        meta: meta,
        trait: file("${trait}"),
        eqtl: file("${eqtl}"),
        expression: file("${expression}"),
        genotype_matrix: file("${genotype_matrix}"),
        medmulti: file("${medmulti}"),
        medsingle: file("${meta.id}_med.tsv"),
        genes: file("${genes}"),
    )

    topic:
    record(tool:"R", version:eval("R --version |& grep 'R version' | cut -f 3 -d' '")) >> 'versions'
    record(tool:"R-MultiMed", version:eval("Rscript -e 'packageVersion(\"MultiMed\")' |& cut -f 2 -d\" \" |& sed 's/’//g' |& sed 's/‘//g'")) >> 'versions'

    script:
    """
    touch mediation_results.tsv
    single_mediation.R \\
        ${genes} \\
        ${genotype_matrix} \\
        ${expression} \\
        ${eqtl} \\
        ${trait} \\
        ${meta.chrom} \\
        ${meta.peak} \\
        ${meta.trait} \\
        ${meta.id}_med.tsv
    """

    stub:
    """
    touch ${meta.id}_med.tsv
    """
}