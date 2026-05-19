nextflow.enable.types = true

process GCTA64_MAKE_GRM {
    tag "${meta.id}"
    maxRetries 0
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "docker://quay.io/biocontainers/gcta:1.94.1--h9ee0642_0"

    input:
    record(
        meta: Map,
        bed: Path,
        bim: Path,
        fam: Path
    )
    maf: BigDecimal
    sparse_cut: BigDecimal

    output:
    record(
        meta: meta + [
            grm_prefix: "${meta.plink_prefix}_gcta_grm_${meta.method}",
            sparse_prefix: "${meta.plink_prefix}_sparse_grm_${meta.method}"
        ],
        bed: file("${bed}"),
        bim: file("${bim}"),
        fam: file("${fam}"),
        grm_bin: file("${meta.plink_prefix}_gcta_grm_${meta.method}.grm.bin"),
        grm_id: file("${meta.plink_prefix}_gcta_grm_${meta.method}.grm.id"),
        grm_N: file("${meta.plink_prefix}_gcta_grm_${meta.method}.grm.N.bin", optional: true),
        sparse_sp: file("${meta.plink_prefix}_sparse_grm_${meta.method}.grm.sp", optional: true),
        sparse_id: file("${meta.plink_prefix}_sparse_grm_${meta.method}.grm.id", optional: true),
    )

    topic:
    record(tool:"gcta64", version:eval("gcta64 | grep version | cut -f 3 -d' '")) >> 'versions'

    script:
    def command = meta.method == "inbred" ? "--make-grm-inbred" : "--make-grm" 
    """
    gcta64 --bfile ${meta.plink_prefix} \\
           --autosome \\
           --maf ${maf} \\
           ${command} \\
           --out ${meta.plink_prefix}_gcta_grm_${meta.method} \\
           --thread-num ${task.cpus}

    if [[ "${meta.method}" == "inbred" ]]; then
        gcta64 --grm ${meta.plink_prefix}_gcta_grm_${meta.method} \\
               --make-bK-sparse ${sparse_cut} \\
               --out ${meta.plink_prefix}_sparse_grm_${meta.method} \\
               --thread-num ${task.cpus}
    fi
    """

    stub:
    """
    touch ${meta.plink_prefix}_gcta_grm_${meta.method}.grm.bin
    touch ${meta.plink_prefix}_gcta_grm_${meta.method}.grm.id
    if [[ ${meta.method} == "inbred" ]]; then
        touch ${meta.plink_prefix}_sparse_grm_${meta.method}.grm.bin
        touch ${meta.plink_prefix}_sparse_grm_${meta.method}.grm.id
    fi
    if [[ ${meta.method} == "loco" ]]; then
        touch ${meta.plink_prefix}_gcta_grm_${meta.method}.N.grm.bin
    fi
    """
}