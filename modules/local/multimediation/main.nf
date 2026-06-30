nextflow.enable.types = true

process MULTIMEDIATION {
    tag "${meta.id}"
    label "multimediation"
    label 'process_small'

    conda "${moduleDir}/environment.yml"
    container "andersenlab/mediation:20220407173046db3227"

    input:
    record (
        meta: Map,
        trait: Path,
        eqtl: Path,
        expression: Path,
        genotype_matrix: Path
    )

    output:
    record (
        meta: meta,
        trait: file("${trait}"),
        eqtl: file("${eqtl}"),
        expression: file("${expression}"),
        genotype_matrix: file("${genotype_matrix}"),
        medmulti: file("${meta.id}_medmulti.tsv"),
        genes: file("${meta.id}_genes.tsv")
    )

    topic:
    record(tool:"R", version:eval("R --version |& grep 'R version' | cut -f 3 -d' '")) >> 'versions'
    record(tool:"R-MultiMed", version:eval("Rscript -e 'packageVersion(\"MultiMed\")' |& cut -f 2 -d\" \" |& sed 's/’//g' |& sed 's/‘//g'")) >> 'versions'

    script:
    """
    touch ${meta.id}_genes.tsv
    multimediation.R \\
        ${genotype_matrix} \\
        ${expression} \\
        ${eqtl} \\
        ${trait} \\
        ${meta.chrom} \\
        ${meta.peak} \\
        ${meta.id}
    """

    stub:
    """
    touch ${meta.id}_genes.tsv
    touch ${meta.id}_medmulti.tsv
    """
}