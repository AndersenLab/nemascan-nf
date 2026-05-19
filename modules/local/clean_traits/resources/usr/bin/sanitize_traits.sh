#!/usr/bin/env bash

# sanitize_traits.sh traits.tsv
traits=$1

# Remove illegal characters from trait names
head -n 1 ${traits} | tr -c '[:alnum:]\t\n' '_' > cleaned_traits.tsv

# Ensure all trait values are decimals or "NA" 
tail -n +2 ${traits} | awk '
{
    printf "%s", $1;
    for (I=2; I<=NF; I++) {
        # Check if value is already a decimal
        if ($I ~ /^-?[0-9]*(\.[0-9]*)?$/) {
            printf "\t%f", $I;
        # Check if value is in scientific notation
        } else if ($I ~ /^-?[0-9](\.[0-9]*)?*[e|E]-?[0-9]+$/) {
            printf "\t%f", $I;
        } else {
            printf "\tNA";
        }
    }
    printf "\n";
}' >> cleaned_traits.tsv
