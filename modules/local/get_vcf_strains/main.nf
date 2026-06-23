nextflow.enable.types = true

process GET_VCF_STRAINS {
    label "local"

    conda null
    container null

    input:
        record(vcf: Path)

    output:
    file("vcf_strains.txt")

    script:
    def command = vcf.extension == "gz" ? 'zcat' : 'cat'
    """
    # Check if on a Mac
    if [ \$(uname) == "Darwin" ] && [ "${command}" == "zcat" ]; then
        COMMAND="gunzip -c"
    else
        COMMAND="${command}"
    fi

    \${COMMAND} ${vcf} | \
        head -n 400 | \
        grep -w "#CHROM" | \
        awk '{ for (I=10;I<=NF;I++) printf "%s\\n", \$I}' > vcf_strains.txt
    """

    stub:
    """
    touch vcf_strains.txt
    """
}