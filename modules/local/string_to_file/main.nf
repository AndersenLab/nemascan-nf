nextflow.enable.types = true

process STRING_TO_FILE {
    tag "${name}"

    label "local"
    maxRetries 0

    conda null
    container null

    input:
    record (
        name: String,
        value: String
    ) 

    output:
    file("${name}")

    script:
    """
    echo -e "${value}" > ${name}
    """

    stub:
    """
    touch ${name}
    """
}