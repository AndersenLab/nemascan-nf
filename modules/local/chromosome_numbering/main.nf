nextflow.enable.types = true

process CHROMOSOME_NUMBERING {
    label "local"

    conda null
    container null

    input:
    record(
        vcf: Path,
        vcf_index: Path
    )

    output:
    file("chrom_numbering.txt")

    script:
    """
    zcat ${vcf} | head -n 200 | grep "##contig" | \\
        awk '
        function roman_to_arabic(roman_str) {
            len = length(roman_str)
            arabic_num = 0
            for (i = 1; i <= len; i++) {
                char = substr(roman_str, i, 1)
                next_char = substr(roman_str, i+1, 1)

                val = 0
                if (char == "I") val = 1
                else if (char == "V") val = 5
                else if (char == "X") val = 10
                else if (char == "L") val = 50
                else if (char == "C") val = 100
                else if (char == "D") val = 500
                else if (char == "M") val = 1000

                # Handle the subtraction rule (e.g., IV = 4)
                if (i < len) {
                    next_val = 0
                    if (next_char == "V") next_val = 5
                    else if (next_char == "X") next_val = 10
                    else if (next_char == "L") next_val = 50
                    else if (next_char == "C") next_val = 100
                    else if (next_char == "D") next_val = 500
                    else if (next_char == "M") next_val = 1000

                    if (val < next_val) {
                        arabic_num += (next_val - val)
                        i++ # Skip the next character as it was part of the subtraction
                        continue
                    }
                }
                arabic_num += val
            }
            return arabic_num
        }{
            split(\$1,A,"=");
            split(A[3],B,",");
            CHROM=B[1];
            if (CHROM ~ /^(Chr|chr|CHR)?[IVXLCDM]*\$/){
                STRIPPED=CHROM;
                gsub("(Chr|chr|CHR)", "", STRIPPED);
                NUM=roman_to_arabic(STRIPPED);
            } else NUM=99;
            printf "%s\\t%s\\n", NUM, CHROM;
        }' | sort -k1,1n |
        awk '
        BEGIN{
            printf "chromosome\\tchrom_number\\n";
        }{
            printf "%s\\t%s\\n", \$2, NR;
        }' > chrom_numbering.txt
    """

    stub:
    """
    echo -e "chromosome\\tchrom_number\\n" > chrom_numbering.txt
    """
}