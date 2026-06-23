nextflow.enable.types = true

process GCTA64_PCA {
    tag "${meta.method}"
    label 'process_small_progressive'

    conda "${moduleDir}/environment.yml"
    container "docker://quay.io/biocontainers/gcta:1.94.1--h9ee0642_0"

    input:
    record(
        meta: Map,
        bed: Path,
        bim: Path,
        fam: Path,
        grm_bin: Path,
        grm_id: Path,
        grm_N: Path?,
        sparse_sp: Path?,
        sparse_id: Path?,
    )

    output:
    record(
        meta: meta,
        bed: file("${bed}"),
        bim: file("${bim}"),
        fam: file("${fam}"),
        grm_bin: file("${grm_bin}"),
        grm_id: file("${grm_id}"),
        grm_N: file("${grm_N}"),
        sparse_sp: file("${sparse_sp}", optional: true),
        sparse_id: file("${sparse_id}", optional: true),
        pca: file("${meta.grm_prefix}.eigenvec")
    )

    topic:
    record(tool:"gcta64", version:eval("gcta64 | grep version | cut -f 3 -d' '")) >> 'versions'

    script:
    """
    gcta64 --grm ${meta.grm_prefix} \\
           --pca 1 \\
           --out ${meta.grm_prefix} \\
           --thread-num ${task.cpus}
    """

    stub:
    """
    touch ${meta.grm_prefix}.eigenvec
    """
}