nextflow.enable.types = true

process AGGREGATE_MAPPINGS {
    tag "${meta.id}"
    label "process_small"

    conda "${moduleDir}/environment.yml"
    container "andersenlab/numpy:2025071813435349b371"

    input:
    record(
        meta: Map,
        trait: Path,
        gwa: Path,
        chrom_numbering: Path
    )
    independent_tests: Double
    snp_grouping: Integer
    ci_size: Integer
    significance_threshold: String

    output:
    record(
        meta: meta,
        trait: file("${trait}"),
        gwa: file("${meta.id}_chroms.gwa"),
        qtl: file("${meta.id}_qtl.tsv")
    )

    topic:
    record(tool:"python", version:eval("python --version |& sed '1!d; s/^.*Python //'")) >> 'versions'
    record(tool:"numpy", version:eval("python -c 'import numpy; print(numpy.version.version)'")) >> 'versions'

    script:
    def independent_test_option = independent_tests ? "--independent_tests ${independent_tests}" : ''
    """
    aggregate_mappings.py \\
        ${gwa} \\
        ${independent_test_option} \\
        --snp_grouping ${snp_grouping} \\
        --CI_size ${ci_size} \\
        --significance_threshold ${significance_threshold} \\
        --method ${meta.method} \\
        --trait ${meta.trait_name} \\
        --chromosome_names ${chrom_numbering} \\
        --output ${meta.id}_qtl.tsv

    awk 'BEGIN{OFS="\t"}{
        if (NR == FNR) CHROM_NAME[\$2]=\$1;
        else if (FNR == 1) print \$0;
        else {
            \$1=CHROM_NAME[\$1];
            \$2=\$1 ":" \$3;
            print \$0;
        }
    }' ${chrom_numbering} ${gwa} > ${meta.id}_chroms.gwa
    """

    stub:
    """
    touch "${meta.id}_qtl.tsv"
    touch "${meta.id}_chroms.gwa"
    """
}