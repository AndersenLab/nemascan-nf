nextflow.enable.types = true

process COLLECT_VERSIONS {
    input:
    version_lines: Bag<String>

    output:
    file('versions.txt')

    script:
    def content = version_lines.join('\n')
    """
    echo "${content}" > versions.txt
    """

    stub:
    """
    touch versions.txt
    """
}
