#!/usr/bin/env bash

# summarize_traits.sh strain_mapping.tsv traits.tsv summarization_method
traits=$1
strain_mapping=$2
summarization_method=$3

# Rename/remove samples and summarize replicates

awk -v METHOD="${summarization_method}" '
function topDownMerge(B, start, middle, end, A) {
    i = start;
    j = middle;
    for (k=start;k<end;k++) {
        if (i < middle && (j >= end || A[i] < A[j])) {
            B[k] = A[i];
            i = i + 1;
        } else {
            B[k] = A[j];
            j = j + 1;
        }
    }
}

function topDownSplitMerge(B, start, stop, A) {
    if (stop - start <= 1) return;
    middle = int((stop + start) / 2);
    topDownSplitMerge(A, start, middle, B);
    middle = int((stop + start) / 2);
    topDownSplitMerge(A, middle, stop, B);
    middle = int((stop + start) / 2);
    topDownMerge(B, start, middle, stop, A);
}

function sort(ORIGINAL, SORTED) {
    N = length(ORIGINAL);
    for (K=1;K<=N;K++) {
        WORK[K] = ORIGINAL[K];
        SORTED[K] = ORIGINAL[K];
    }
    topDownSplitMerge(SORTED, 1, N+1, WORK);
}

function mean(ARRAY) {
    sum = 0;
    n = 0;
    for (K in ARRAY){
        if (K != "NA") {
            sum = sum + ARRAY[K];
            n = n + 1;
        }
    }
    if (n > 0) {
        return sum / n;
    } else {
        return "NA";
    }
}

function median(ARRAY) {
    n = 0
    for (K in ARRAY){
        if (ARRAY[K] != "NA") {
            n = n + 1;
            TEMP[n] = ARRAY[K];
        }
    }
    if (n == 0) {
        return "NA";
    } else {
        sort(TEMP, SORTED);
        mid = n / 2;
        if (mid == int(mid)) {
            return (SORTED[int(mid)] + SORTED[int(mid) + 1]) / 2;
        } else {
            return SORTED[int(mid) + 1];
        }
    }
}

{
    if (FNR == NR) {
        STRAIN_MAP[$1] = $2;
    } else {
        if (FNR == 1) {
            for (I=2; I<=NF; I++) {
                TRAITS[I-1] = $I;
            }
            N_TRAITS = length(TRAITS);
        } else {
            if ($1 in STRAIN_MAP) {
                ISOTYPE = STRAIN_MAP[$1];
                if (ISOTYPE in ISOTYPE_INDICES) {
                    I = ISOTYPE_INDICES[ISOTYPE];
                    ISOTYPE_COUNT[ISOTYPE] = ISOTYPE_COUNT[ISOTYPE] + 1;
                    for (J=2; J<=NF; J++) {
                        DATA[I " " (J - 1)] = DATA[I " " (J - 1)] " " $J;
                    }
                } else {
                    I = length(ISOTYPE_INDICES) + 1;
                    ISOTYPE_MAP[I] = ISOTYPE;
                    ISOTYPE_INDICES[ISOTYPE] = I;
                    ISOTYPE_COUNT[ISOTYPE] = 1;
                    for (J=2; J<=NF; J++) {
                        DATA[I " " (J - 1)] = $J;
                    }
                }
            }
        }
    }
}END{
    LINE = "strain";
    for (J=1; J<=N_TRAITS; J++){
        LINE = LINE "\t" TRAITS[J];
    }
    printf "%s\n", LINE;
    for (I=1; I<=length(ISOTYPE_INDICES); I++){
        ISOTYPE = ISOTYPE_MAP[I];
        if (ISOTYPE_COUNT[ISOTYPE] == 1){
            LINE = ISOTYPE;
            for (J=1; J<=N_TRAITS; J++){
                LINE = LINE "\t" DATA[I " " J];
            }
            printf "%s\n", LINE;
        } else if (ISOTYPE_COUNT[ISOTYPE] > 1){
            LINE = ISOTYPE;
            for (J=1; J<=N_TRAITS; J++){
                split(DATA[I " " J],SAMPLE_DATA," ");
                if (METHOD == "median") {
                    SUMMARY = median(SAMPLE_DATA);
                } else {
                    SUMMARY = mean(SAMPLE_DATA);
                }
                LINE = LINE "\t" SUMMARY;
            }
            printf "%s\n", LINE;
        }
    }
}' ${strain_mapping} ${traits} > summarized_traits.tsv
