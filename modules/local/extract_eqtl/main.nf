nextflow.enable.types = true

process EXTRACT_EQTL {
    tag "${meta.id}"

    label "local"

    conda null
    container null

    input:
    record (
        meta: Map,
        trait: Path,
        eqtl: Path,
        expression: Path
    )

    output:
    record (
        meta: meta,
        trait: file("${trait}"),
        eqtl: file("${meta.id}_eqtl.tsv"),
        expression: file("${expression}")
    )

    script:
    """
    awk -v CHROM="${meta.chrom}" -v START="${meta.start - 1e6}" -v STOP="${meta.end + 1e6}" '{
        if (NR == 1) print \$0;
        else if (\$5 == CHROM && \$6 < STOP && \$7 > START) print \$0;
    }' ${eqtl} > ${meta.id}_eqtl.tsv
    """

    stub:
    """
    touch ${meta.id}_eqtl.tsv
    """
}