#!/usr/bin/env bash

# map_strains_to_isotypes.sh isogroups.tsv traits.tsv vcf_strains.txt

isogroups=$1
traits=$2
vcf_strains=$3

if [[ $(wc -l ${isogroups} | awk '{print $1}') -eq 0 ]]; then
    awk '{
        if (FNR == NR) {
            VCF_STRAINS[$1] = 1;
        } else {
            if (NR > 1) {
                if ($1 in VCF_STRAINS) {
                    printf "%s\t%s\n", $1, $1 >> "strain_mapping.tsv";
                } else {
                    if (length(VCF_MISSING_ISOTYPE) == 0){
                        VCF_MISSING_ISOTYPE = STRAIN;
                    } else {
                        VCF_MISSING_ISOTYPE = VCF_MISSING_ISOTYPE + ", " + STRAIN;
                    }
                }
            }
        }
    }END{
        if (length(VCF_MISSING_ISOTYPE) > 0) {
            printf "WARNING: Removing strain(s) %s because they do not appear in the VCF samples.\n", VCF_MISSING_ISOTYPE >> "strain_issues.txt";
        } else {
            printf "No Issues" >> "strain_issues.txt";
        }
    }' ${vcf_strains} ${traits}
else
    awk 'BEGIN{
        FILE = 0;
    }{
        if (FNR == 1) FILE = FILE + 1;
        if (FILE == 1) {
            VCF_STRAINS[$1] = 1;
        } else if (FILE == 2) {
            if (FNR == 1) {
                # Determine which columns have the strain and isotype names
                for (INDEX=1; INDEX<=NF; INDEX++) {
                    if ($INDEX == "strain") {
                        STRAIN_INDEX=INDEX;
                    } else if ($INDEX == "isotype") {
                        ISOTYPE_INDEX=INDEX;
                    } else if ($INDEX == "previous_names") {
                        ALT_NAME_INDEX=INDEX;
                    }
                }
            } else {
                ISOTYPES[$STRAIN_INDEX] = $ISOTYPE_INDEX;
                if (length($ALT_NAME_INDEX) > 0) {
                    split($ALT_NAME_INDEX,ALTS,"|");
                    for (I in ALTS) {
                        ALT_NAMES[ALTS[I]]=$STRAIN_INDEX;
                    }
                }
            }
        } else {
            if (FNR > 1 && length($1) > 0) {
                if ($1 in ALT_NAMES) {
                    STRAIN=ALT_NAMES[$1];
                } else {
                    STRAIN=$1;
                }
                if (STRAIN in ISOTYPES) {
                    ISOTYPE = ISOTYPES[STRAIN]
                    MAPPING[$1] = ISOTYPE;
                } else {
                    MAPPING[$1] = "NA"
                }
            }
        }
    }END{
        ISSUES = 0
        for (STRAIN in MAPPING) {
            ISOTYPE = MAPPING[STRAIN];
            if (ISOTYPE in VCF_STRAINS || ISOTYPE == "NA"){
                ISOCOUNT[ISOTYPE] = ISOCOUNT[ISOTYPE] + 1;
                if (ISOTYPE in ISOSTRAINS){
                    ISOSTRAINS[ISOTYPE] = ISOSTRAINS[ISOTYPE] ", " STRAIN;
                }else{
                    ISOSTRAINS[ISOTYPE] = STRAIN;
                }
            } else {
                ISSUES = 1;
                if (STRAIN == ISOTYPE){
                    if (length(VCF_MISSING_ISOTYPE) == 0){
                        VCF_MISSING_ISOTYPE = STRAIN;
                    } else {
                        VCF_MISSING_ISOTYPE = VCF_MISSING_ISOTYPE ", " STRAIN;
                    }
                }else{
                    if (ISOTYPE in VCF_MISSING_STRAIN){
                        VCF_MISSING_STRAIN[ISOTYPE] = STRAIN;
                    } else {
                        VCF_MISSING_STRAIN[ISOTYPE] = VCF_MISSING_STRAIN[ISOTYPE] ", " STRAIN;
                    }
                }
            }
        }
        for (ISOTYPE in ISOCOUNT){
            if (ISOTYPE == "NA"){
                NO_ISOTYPE = ISOSTRAINS["NA"];
                ISSUES = 1;
            } else if (ISOCOUNT[ISOTYPE] == 1){
                STRAIN = ISOSTRAINS[ISOTYPE];
                if (ISOTYPE != STRAIN) {
                    ISSUES = 1;
                    RENAMED[STRAIN] = ISOTYPE;
                }
                printf "%s\t%s\n", STRAIN, ISOTYPE >> "strain_mapping.tsv";
            } else {
                ISSUES = 1;
                split(ISOSTRAINS[ISOTYPE],STRAINS,", ");
                IS_ISOTYPE = 0;
                for (IDX in STRAINS){
                    if (STRAINS[IDX] == ISOTYPE){
                        IS_ISOTYPE = 1;
                    }
                }
                if (IS_ISOTYPE > 0){
                    printf "%s\t%s\n", ISOTYPE, ISOTYPE >> "strain_mapping.tsv";
                    for (IDX in STRAINS){
                        if (STRAINS[IDX] != ISOTYPE){
                            if (ISOTYPE in EXTRA_STRAINS) {
                                EXTRA_STRAINS[ISOTYPE] = EXTRA_STRAINS[ISOTYPE] ", " STRAINS[IDX];
                            } else {
                                EXTRA_STRAINS[ISOTYPE] = STRAINS[IDX];
                            }
                        }
                    }
                } else {
                    for (IDX in STRAINS){
                        if (ISOTYPE in MULTIPLE_STRAINS) {
                            MULTIPLE_STRAINS[ISOTYPE] = MULTIPLE_STRAINS[ISOTYPE] ", " STRAINS[IDX];
                        } else {
                            MULTIPLE_STRAINS[ISOTYPE] = STRAINS[IDX];
                        }
                    }
                }
            }
        }
        if (length(NO_ISOTYPE) > 0) {
            printf "WARNING: Removing strain(s) %s because they do not fall into a defined isotype.\n", NO_ISOTYPE >> "strain_issues.txt";
        }
        if (length(VCF_MISSING_ISOTYPE) > 0) {
            printf "WARNING: Removing isotype reference strain(s) %s because they do not appear in the VCF samples.\n", VCF_MISSING_ISOTYPE >> "strain_issues.txt";
        }
        for (ISOTYPE in VCF_MISSING_STRAIN) {
            printf "WARNING: Removing non-isotype reference strain %s because its isotype %s does not appear in the VCF samples.\n", VCF_MISSING_STRAIN[ISOTYPE], ISOTYPE >> "strain_issues.txt";
        }
        for (ISOTYPE in EXTRA_STRAINS) {
            printf "WARNING: Removing non-isotype reference strain(s) %s from isotype group %s.\n", EXTRA_STRAINS[ISOTYPE], ISOTYPE >> "strain_issues.txt";
        }
        for (ISOTYPE in MULTIPLE_STRAINS) {
            printf "WARNING: Removing non-isotype reference strain(s) %s from isotype group %s. To include this isotype in the analysis, evaluate the similarity of these strains and choose one representative for the group.\n", MULTIPLE_STRAINS[ISOTYPE], ISOTYPE >> "strain_issues.txt";
        }
        for (STRAIN in RENAMED) {
            printf "NOTE: Non-isotype reference strain %s renamed to isotype %s.\n", STRAIN, RENAMED[STRAIN] >> "strain_issues.txt";
        }
    }' ${vcf_strains} ${isogroups} ${traits}
fi
