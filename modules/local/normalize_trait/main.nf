nextflow.enable.types = true

process NORMALIZE_TRAIT {
    tag "${meta.id}"

    label "local"
    maxRetries 0

    conda null
    container null

    input:
    record (
        meta: Map,
        trait: Path,
        qtl: Path
    )

    output:
    record (
        meta: meta,
        trait: file("${meta.trait_name}_norm.tsv"),
        qtl: file("${qtl}")
    )

    script:
    """
    awk '{
        if (\$2 ~ /^-?[0-9]*(.[0-9]*)?\$/) {
            TRAIT[\$1] = \$2;
            N += 1;
            MEAN += \$2;
            STD += \$2 * \$2;
        }
    }END{
        MEAN /= N;
        STD = (STD / (N - 1) - MEAN * MEAN) ^ 0.5;
        for (STRAIN in TRAIT) {
            printf "%s\\t%s\\n", STRAIN, (TRAIT[STRAIN] - MEAN) / STD;
        }
    }' ${trait} > ${meta.trait_name}_norm.tsv
    """

    stub:
    """
    touch ${meta.trait_name}_norm.tsv
    """
}