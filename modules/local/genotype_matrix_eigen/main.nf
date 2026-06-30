nextflow.enable.types = true

process GENOTYPE_MATRIX_EIGEN {
    label 'process_bigmem_progressive'
    tag "${chromosome}"

    conda "${moduleDir}/environment.yml"
    container "andersenlab/numpy:2025071813435349b371"

    input:
    record(
        genotype_matrix: Path,
        chromosome: String
    )

    output:
    file("${chromosome}_independent_snvs.txt")

    topic:
    record(tool:"python", version:eval("python --version |& sed '1!d; s/^.*Python //'")) >> 'versions'
    record(tool:"numpy", version:eval("python -c 'import numpy; print(numpy.version.version)'")) >> 'versions'

    script:
    """
    awk -v CHROM="${chromosome}" '{if (NR == 0) print \$0; else if (\$1 == CHROM){gsub("NA", "0", \$0); print \$0;}}' ${genotype_matrix} > chr_genotype_matrix.tsv
    genotype_matrix_eigen.py chr_genotype_matrix.tsv > ${chromosome}_independent_snvs.txt
    """

    stub:
    """
    touch ${chromosome}_independent_snvs.txt
    """
}