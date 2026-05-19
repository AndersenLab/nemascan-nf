#!/usr/bin/env bash

# remove_outliers.sh traits.tsv skip_pruning
# Returns filtered_strains.tsv, included_strains.txt, omitted_strains.txt 

traits=$1
skip_pruning=$2

if ( ${skip_pruning} == "true" ); then
    awk '{
        if (NR == 1) print $0;
        else {
            valid = 0;
            for (I=2; I<=NF; I++) {
                if ($I != "NA") valid = valid + 1;
            }
            if (valid > 0) print $0;
        }
    }' ${traits} > filtered_traits.tsv
else
    awk '
    function quantile(X, p) {
        n = length(X);
        h = (n - 1) * p + 1;
        j = int(h);
        gamma = h - j;
        q = (1 - gamma) * X[j] + gamma * X[j+1];
        return q;
    }

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
        for (I=1;I<=N;I++) WORK[I] = ORIGINAL[I];
        for (I=1;I<=N;I++) SORTED[I] = ORIGINAL[I];
        topDownSplitMerge(SORTED, 1, N+1, WORK);
    }

    {
        if (NR == 1) {
            for (I=2;I<=NF;I++) {
                TRAITS[I-1] = $I;
            }
            N_TRAITS = length(TRAITS);
        } else {
            I = NR - 1;
            STRAINS[I] = $1;
            for (J=2;J<=NF;J++){
                ALL_DATA[I " " (J-1)] = $J;
            }
        }
    }END{
        N_STRAINS = length(STRAINS);
        for (J=1;J<=N_TRAITS;J++) {
            delete DATA;
            delete DATA_INDEX;
            delete SORTED;
            N = 0;
            for (I=1;I<=N_STRAINS;I++) {
                if (ALL_DATA[I " " J] != "NA") {
                    N = N + 1;
                    DATA[N] = ALL_DATA[I " " J];
                    DATA_INDEX[N] = I;
                } else {
                    OUTLIER[I " " J] = 1;
                }
            }
            sort(DATA, SORTED);
            Q1=quantile(SORTED, 0.25);
            Q3=quantile(SORTED, 0.75);
            IQR=Q3 - Q1;
            CUTOFFS["6L"] = Q1 - IQR * 10;
            CUTOFFS["5L"] = Q1 - IQR * 7;
            CUTOFFS["4L"] = Q1 - IQR * 5;
            CUTOFFS["3L"] = Q1 - IQR * 4;
            CUTOFFS["2L"] = Q1 - IQR * 3;
            CUTOFFS["1L"] = Q1 - IQR * 2;
            CUTOFFS["1H"] = Q3 + IQR * 2;
            CUTOFFS["2H"] = Q3 + IQR * 3;
            CUTOFFS["3H"] = Q3 + IQR * 4;
            CUTOFFS["4H"] = Q3 + IQR * 5;
            CUTOFFS["5H"] = Q3 + IQR * 7;
            CUTOFFS["6H"] = Q3 + IQR * 10;
            delete BINCOUNTS;
            for (I=1;I<=N;I++) {
                if (DATA[I] <=  CUTOFFS["6L"]) {
                    BINCOUNTS["SIXLS"] = BINCOUNTS["SIXLS"] + 1;
                    SIXLS[I] = 1;
                } else {
                    SIXLS[I] = 0;
                }
                if (DATA[I] >  CUTOFFS["6L"] && DATA[I] <= CUTOFFS["5L"]) {
                    BINCOUNTS["FIVELS"] = BINCOUNTS["FIVELS"] + 1;
                    FIVELS[I] = 1;
                } else {
                    FIVELS[I] = 0;
                }
                if (DATA[I] >  CUTOFFS["5L"] && DATA[I] <= CUTOFFS["4L"]) {
                    BINCOUNTS["FOURLS"] = BINCOUNTS["FOURLS"] + 1;
                    FOURLS[I] = 1;
                } else {
                    FOURLS[I] = 0;
                }
                if (DATA[I] >  CUTOFFS["4L"] && DATA[I] <= CUTOFFS["3L"]) {
                    BINCOUNTS["THREELS"] = BINCOUNTS["THREELS"] + 1;
                    THREELS[I] = 1;
                } else {
                    THREELS[I] = 0;
                }
                if (DATA[I] >  CUTOFFS["3L"] && DATA[I] <= CUTOFFS["2L"]) {
                    BINCOUNTS["TWOLS"] = BINCOUNTS["TWOLS"] + 1;
                    TWOLS[I] = 1;
                } else {
                    TWOLS[I] = 0;
                }
                if (DATA[I] >  CUTOFFS["2L"] && DATA[I] <= CUTOFFS["1L"]) {
                    BINCOUNTS["ONELS"] = BINCOUNTS["ONELS"] + 1;
                    ONELS[I] = 1;
                } else {
                    ONELS[I] = 0;
                }
                if (DATA[I] >= CUTOFFS["1H"] && DATA[I] < CUTOFFS["2H"]) {
                    BINCOUNTS["ONEHS"] = BINCOUNTS["ONEHS"] + 1;
                    ONEHS[I] = 1;
                } else {
                    ONEHS[I] = 0;
                }
                if (DATA[I] >= CUTOFFS["2H"] && DATA[I] < CUTOFFS["3H"]) {
                    BINCOUNTS["TWOHS"] = BINCOUNTS["TWOHS"] + 1;
                    TWOHS[I] = 1;
                } else {
                    TWOHS[I] = 0;
                }
                if (DATA[I] >= CUTOFFS["3H"] && DATA[I] < CUTOFFS["4H"]) {
                    BINCOUNTS["THREEHS"] = BINCOUNTS["THREEHS"] + 1;
                    THREEHS[I] = 1;
                } else {
                    THREEHS[I] = 0;
                }
                if (DATA[I] >= CUTOFFS["4H"] && DATA[I] < CUTOFFS["5H"]) {
                    BINCOUNTS["FOURHS"] = BINCOUNTS["FOURHS"] + 1;
                    FOURHS[I] = 1;
                } else {
                    FOURHS[I] = 0;
                }
                if (DATA[I] >= CUTOFFS["5H"] && DATA[I] < CUTOFFS["6H"]) {
                    BINCOUNTS["FIVEHS"] = BINCOUNTS["FIVEHS"] + 1;
                    FIVEHS[I] = 1;
                } else {
                    FIVEHS[I] = 0;
                }
                if (DATA[I] >= CUTOFFS["6H"]) {
                    BINCOUNTS["SIXHS"] = BINCOUNTS["SIXHS"] + 1;
                    SIXHS[I] = 1;
                } else {
                    SIXHS[I] = 0;
                }
            }
            for (I=1;I<=N;I++) {
                cuts1 = (((SIXHS[I] == 1) && (((BINCOUNTS["SIXHS"] + BINCOUNTS["FIVEHS"] + BINCOUNTS["FOURHS"]) / N_STRAINS) <= 0.05)) || \
                        ((SIXLS[I] == 1) && (((BINCOUNTS["SIXLS"] + BINCOUNTS["FIVELS"] + BINCOUNTS["FOURLS"]) / N_STRAINS) <= 0.05)));
                cuts2 = (((FIVEHS[I] == 1) && (((BINCOUNTS["SIXHS"] + BINCOUNTS["FIVEHS"] + BINCOUNTS["FOURHS"] + BINCOUNTS["THREEHS"]) / N_STRAINS) <= 0.05)) || \
                        ((FIVELS[I] == 1) && (((BINCOUNTS["SIXLS"] + BINCOUNTS["FIVELS"] + BINCOUNTS["FOURLS"] + BINCOUNTS["THREELS"]) / N_STRAINS) <= 0.05)));
                cuts3 = (((FOURHS[I] == 1) && (((BINCOUNTS["SIXHS"] + BINCOUNTS["FIVEHS"] + BINCOUNTS["FOURHS"] + BINCOUNTS["THREEHS"] + BINCOUNTS["TWOHS"]) / N_STRAINS) <= 0.05)) || \
                        ((FOURLS[I] == 1) && (((BINCOUNTS["SIXLS"] + BINCOUNTS["FIVELS"] + BINCOUNTS["FOURLS"] + BINCOUNTS["THREELS"] + BINCOUNTS["TWOLS"]) / N_STRAINS) <= 0.05)));
                if (cuts1 == 1 || cuts2 == 1 || cuts3 == 1) {
                    OUTLIER[DATA_INDEX[I] " " J] = 1;
                } else {
                    OUTLIER[DATA_INDEX[I] " " J] = 0;
                }
            }
        }
        for (I=1;I<=N_TRAITS;I++) {
            TRAIT_NAMES = TRAIT_NAMES "\t" TRAITS[I];
        }
        printf "strain%s\n", TRAIT_NAMES;
        for (I=1;I<=N_STRAINS;I++) {
            LINE = STRAINS[I];
            N_VALUES = 0
            for (J=1;J<=N_TRAITS;J++) {
                if (OUTLIER[I " " J] == 0) {
                    LINE = LINE "\t" ALL_DATA[I " " J];
                    N_VALUES = N_VALUES + 1;
                } else {
                    LINE = LINE "\tNA";
                }
            }
            if (N_VALUES > 0) {
                printf "%s\n", LINE;
            }
        }
    }' ${traits} > filtered_traits.tsv
fi

cut -f 1 filtered_traits.tsv | tail -n +2 | sort -k1,1 > included_strains.txt
if [[ $(tail -n +2 filtered_traits.tsv | wc -l | awk '{print $1}') -gt $(wc -l included_strains.txt | awk '{print $1}') ]]; then
    cut -f 1 ${traits} | tail -n +2 | grep -v -w -f included_strains.txt > omitted_strains.txt
else
    touch omitted_strains.txt
fi