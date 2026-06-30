nextflow.enable.types = true

process STRING_TO_FILE {
    tag "${name}"

    label "local"

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