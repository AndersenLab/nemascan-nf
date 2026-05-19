nextflow.enable.types = true

process APPEND_ANNOTATIONS {
    tag "${meta.id}"
    label 'process_small'

    conda "${moduleDir}/environment.yml"
    container "docker://andersenlab/numpy:2025071813435349b371"

    input:
    record(
        meta: Map,
        gwa: Path,
        ld: Path,
        annotations: Path
    )

    output:
    record(
        meta: meta,
        gwa: file("${meta.id}.annotated.gwa")
    )

    topic:
    record(tool:"python", version:eval("python --version |& sed '1!d; s/^.*Python //'")) >> 'versions'

    script:
    """
    HEADER=`head -n 1 ${gwa}`
    echo -e "\${HEADER}\tLD" > tmp.LD.gwa
    tail -n +2 ${gwa} | cut -f 2-10 | sort -k1,1 > tmp.gwa
    sort -k1,1 ${ld} > tmp.ld
    join -e "NA" tmp.gwa tmp.ld | \\
        awk -v CHROM="${meta.chrom}" '{
            \$1=CHROM ":" \$2;
            printf "%s\\t%s\\n", CHROM, \$0
        }' | \\
        sort -k2,2 -k3,3n | \\
        sed 's/ /\\t/g' >> tmp.LD.gwa
        
    append_annotations.py tmp.LD.gwa ${annotations} ${meta.id}.annotated.gwa
    """

    stub:
    """
    touch ${meta.id}.annotated.gwa
    """
}
