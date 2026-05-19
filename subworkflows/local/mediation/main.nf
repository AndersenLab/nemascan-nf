nextflow.enable.types = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { NORMALIZE_TRAIT           } from '../../../modules/local/normalize_trait'
include { EXTRACT_EQTL              } from '../../../modules/local/extract_eqtl'
include { MULTIMEDIATION            } from '../../../modules/local/multimediation'
include { SINGLE_MEDIATION          } from '../../../modules/local/single_mediation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record EqtlRecord {
    eqtl: Path
    expression: Path
}
record QtlRecord {
    trait: String
    marker: String
    method: String
    log10p: String
    CHROM: String
    startPOS: String
    peakPOS: String
    endPOS: String
    peak_id: String
}
record SigQtlRecord {
    meta: Map
    trait: Path
    qtl: Path
}
record MediationRecord {
    meta: Map
    trait: Path
    eqtl: Path
    expression: Path
    genotype_matrix: Path
    medmulti: Path
    medsingle: Path
    genes: Path
}

workflow MEDIATION {

    take:
    ch_genotype_matrix: Channel<Path>
    ch_sig_qtl: Channel<SigQtlRecord>
    ch_eqtl: Channel<EqtlRecord>
    
    main:
    ch_normalize_trait = NORMALIZE_TRAIT (
        ch_sig_qtl
    )

    // Split significant qtl file and extract regions from eQTL file 
    ch_sig_qtl_norm = ch_sig_qtl
        .flatMap { row -> row.qtl.splitCsv ( sep:"\t", header:true ) }
        .map { row -> new QtlRecord(row) }
        .map { row -> record(key:row.trait, qtl_record:record(meta: [
            id:"${row.trait}_${row.method}_${row.CHROM}_${row.startPOS}_${row.endPOS}",
            trait_name:row.trait,
            method:row.method,
            chrom:row.CHROM,
            start:row.startPOS.toInteger(),
            peak:row.peakPOS.toInteger(),
            end:row.endPOS.toInteger()]))
        }
        .join ( ch_normalize_trait
            .map { row -> record(key: row.meta.trait_name, trait:row.trait) },
            by: 'key'
        )
        .combine ( ch_eqtl )
        .map { row -> row[0].qtl_record + record(
            trait: row[0].trait,
            eqtl: row[1].eqtl,
            expression: row[1].expression
        ) }

    ch_extract_eqtl = EXTRACT_EQTL (
        ch_sig_qtl_norm
    )

    ch_multimediation = MULTIMEDIATION (
        ch_extract_eqtl
            .combine ( ch_genotype_matrix )
            .map { row -> row[0] + record(genotype_matrix: row[1]) }
    )

    // Find mediation for each significant gene in list
    ch_single_mediation = SINGLE_MEDIATION (
        ch_multimediation
    )

    emit:
    med_results: Channel<MediationRecord> = ch_single_mediation
}

process TEST {

    label "local"
    maxRetries 0

    conda null
    container null

    input:
    record(
        trait: Path
    )

    output:
    stdout()

    script:
    """
    wc -l ${trait}
    """
}
