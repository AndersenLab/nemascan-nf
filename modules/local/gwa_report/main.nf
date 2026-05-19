nextflow.enable.types = true

record ConfigRecord {
    trait_name: String
    config: String
    paths: List<Path>
    template: Path
}

process GWA_REPORT {
    label "gwa_report"
    maxRetries 0

    conda "${moduleDir}/environment.yml"
    container "docker://andersenlab/plotly:20260417"

    input:
    config: ConfigRecord

    output:
    record (
        trait_name: config.trait_name,
        report: file("${config.trait_name}_report.html")
    )

    topic:
    record(tool:"python", version:eval("python --version |& sed '1!d; s/^.*Python //'")) >> 'versions'

    script:
    """
    echo -e "${config.config}" > ${config.trait_name}_config.tsv


    create_report.py \\
        --trait ${config.trait_name} \\
        --config ${config.trait_name}_config.tsv \\
        --template report_template.html \\
        --output ${config.trait_name}_report.html
    """

    stub:
    """
    touch ${config.trait_name}_report.html
    """
}