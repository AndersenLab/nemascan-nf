nextflow.enable.types = true

process GET_VCF_STRAINS {
    label "local"

    input:
        record(vcf: Path)

    output:
    file("vcf_strains.txt")

    script:
    def command = vcf.extension == "gz" ? 'zcat' : 'cat'
    """
    echo "starting"

    # Check if on a Mac
    if [ \$(uname) == "Darwin" ] && [ "${command}" == "zcat" ]; then
        COMMAND="gunzip -c"
    else
        COMMAND="${command}"
    fi
    echo "command is \${COMMAND}"

    \${COMMAND} ${vcf} | \
        head -n 400 | \
        grep -w "#CHROM" | \
        awk '{ for (I=10;I<=NF;I++) printf "%s\\n", \$I}' > vcf_strains.txt
    cat vcf_strains.txt
    echo "finished"
    """

    stub:
    """
    touch vcf_strains.txt
    """
}