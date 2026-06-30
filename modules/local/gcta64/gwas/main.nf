nextflow.enable.types = true

process GCTA64_GWAS {
    tag "${meta.id}"
    label 'process_small_progressive'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/gcta:1.94.1--h9ee0642_0"

    input:
    record(
        meta: Map,
        bed: Path,
        bim: Path,
        fam: Path,
        grm_bin: Path,
        grm_id: Path,
        grm_N: Path,
        sparse_sp: Path?,
        sparse_id: Path?,
        pca: Path?,
        trait: Path
    )

    output:
    record(
        meta: meta,
        gwa: file("${meta.id}_lmm-exact.gwa"),
        trait: file("${trait}")
    )

    topic:
    record(tool:"gcta64", version:eval("gcta64 | grep version | cut -f 3 -d' '")) >> 'versions'

    script:
    def command = meta.method == "loco" ? "--mlma-loco" : meta.method == "finemap" ? "--mlma" : "--fastGWA-mlm-exact"
    def grm_option = meta.method == "inbred" ? "--grm-sparse ${meta.sparse_prefix}" : "--grm ${meta.grm_prefix}"
    def pca_option = pca ? "--qcovar ${pca}" : "" 
    """
    # linear regression model
    awk '{printf "%s\\t%s\\t%s\\n", \$1, \$1, \$2}' ${trait} > reformatted_traits.tsv
    gcta64 ${command} \\
           ${grm_option} \\
           --bfile ${meta.plink_prefix} \\
           ${pca_option} \\
           --out ${meta.id}_lmm-exact \\
           --pheno reformatted_traits.tsv \\
           --thread-num ${task.cpus}
    
    if [[ ${meta.method} == "loco" ]]; then
        mv ${meta.id}_lmm-exact.loco.mlma ${meta.id}_lmm-exact.gwa
    elif [[ ${meta.method} == "finemap" ]]; then
        mv ${meta.id}_lmm-exact.mlma ${meta.id}_lmm-exact.gwa
    else
        mv ${meta.id}_lmm-exact.fastGWA ${meta.id}_lmm-exact.gwa
    fi
    """

    stub:
    """
    touch ${meta.id}_lmm-exact.gwa
    """
}