nextflow.enable.types = true

process NARROW_H2 {
    tag "${meta.trait_name}"
    maxRetries 0

    conda "${moduleDir}/environment.yml"
    container "docker://andersenlab/r-sommer:20260121"

    input:

    record(
        meta: Map,
        trait: Path,
        genotype_matrix: Path
    )

    output:
    record(
        meta: meta,
        trait: trait,
        matrix: genotype_matrix,
        h2: file("${meta.trait_name}_narrow_h2.txt")
    )

    topic:
    record(tool:"R", version:eval("Rscript  Rscript --version | cut -d' ' -f4")) >> 'versions'
    record(tool:"R-sommer", version:eval("""Rscript -e "library(sommer); cat(as.character(packageVersion('sommer')))" """)) >> 'versions'

    script:
    """
    find_narrow_h2.r ${genotype_matrix} ${trait} ${meta.trait_name} ${meta.trait_name}_narrow_h2.txt
    """

    stub:
    """
    touch "${meta.trait_name}_narrow_h2.txt"
    """
}