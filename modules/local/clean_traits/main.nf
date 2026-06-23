nextflow.enable.types = true

process CLEAN_TRAITS {
    label "local"

    conda null
    container null

    input:
    record(
        traits: Path,
        isogroups: Path,
        vcf_strains: Path,
        summarization_method: String,
        skip_pruning: Boolean
    )

    output:
    record(
        traits: files("split_traits/*.tsv"),
        included: file("included_strains.txt"),
        issues: file("strain_issues.txt")
    )

    script:
    """
    # Clean trait file to remove illegal characters and ensure only decimals and NAs
    # Returns cleaned_traits.tsv
    sanitize_traits.sh ${traits}

    # Create a strain -> isotype reference map and list of strain issues
    # Returns strain_mapping.tsv and strain_issues.txt
    map_strains_to_isotypes.sh ${isogroups} cleaned_traits.tsv ${vcf_strains}

    # Remove samples with missing isotypes or multiple strains per isotype and collapse replicates
    # Returns summarized_traits.tsv
    summarize_traits.sh cleaned_traits.tsv strain_mapping.tsv ${summarization_method}

    # Identify and potentially remove outliers
    # Returns filtered_traits.tsv, included_strains.txt, and omitted_strains.txt
    remove_outliers.sh summarized_traits.tsv ${skip_pruning}

    if [[ \$(wc -l omitted_strains.txt | awk '{print \$1}') -gt 0 ]]; then
        for I in \$(cat omitted_strains.txt); do
            echo "WARNING: Removing isotype reference strain \${I} because it has no valid trait data." >> strain_issues.txt
        done
    fi
    
    if [[ \$(wc -l strain_issues.txt | awk '{print \$1}') -eq 0 ]]; then
        echo "No strain issues found." >> strain_issues.txt
    fi

    # Split traits
    mkdir split_traits
    NAMES=(\$(cat ${traits} | head -n 1))
    for I in \$(seq 1 \$(expr \${#NAMES[*]} - 1)); do
        OUTNAME="\${NAMES[\${I}]}.tsv"
        cut -f1,\$(expr \${I} + 1) filtered_traits.tsv | tail -n +2 | awk '{printf "%s\\t%s\\n", \$1, \$2;}' > split_traits/\${OUTNAME}
    done
    """

    stub:
    """
    mkdir split_traits/dummy.tsv
    touch included_strains.txt
    touch strain_issues.txt
    """
}