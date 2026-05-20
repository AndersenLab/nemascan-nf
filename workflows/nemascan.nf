nextflow.enable.types = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GENOTYPE_MATRIX           } from '../subworkflows/local/genotype_matrix'
include { GWA_MAPPING               } from '../subworkflows/local/gwa_mapping'
include { FINE_MAPPING              } from '../subworkflows/local/fine_mapping'
include { MEDIATION                 } from '../subworkflows/local/mediation'
include { GWA_REPORTING             } from '../subworkflows/local/gwa_reporting'
include { COLLECT_VERSIONS          } from '../modules/local/collect_versions'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record VcfRecord {
    vcf: Path
    vcf_index: Path
}
record EqtlRecord {
    eqtl: Path
    expression: Path
}
record BroadRecord {
    meta: Map
    trait: Path
    gwa: Path
    qtl: Path
    matrix: Path
    h2: Path
    independent_tests: Float
}
record FineMapRecord {
    meta: Map
    gwa: Path
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
record ReportRecord {
    trait_name: String
    report: Path
}
record VersionRecord {
    tool: String
    version: String
}

workflow NEMASCAN {

    take:
    ch_traits: Channel<Path>
    ch_vcf: Channel<VcfRecord>
    ch_imputed_vcf: Channel<VcfRecord>
    ch_annotation: Channel<Path>
    ch_haplotypes: Channel<Path>
    ch_isogroups: Channel<Path>
    ch_eqtl: Channel<EqtlRecord>
    ch_gwa_method: Channel<String>

    main:
    genotype_matrix_call = GENOTYPE_MATRIX (
        ch_traits,
        ch_vcf,
        ch_isogroups
    )
    
    ch_filtered_traits  = genotype_matrix_call.traits
    ch_pruned_vcf       = genotype_matrix_call.pruned_vcf
    ch_strains          = genotype_matrix_call.strains
    ch_strain_issues    = genotype_matrix_call.strain_issues
    ch_genotype_matrix  = genotype_matrix_call.genotype_matrix

    if (params.mapping) {
        gwa_mapping_call = GWA_MAPPING (
            ch_filtered_traits,
            ch_pruned_vcf,
            ch_genotype_matrix,
            ch_gwa_method
        )
        ch_broad_gwa     = gwa_mapping_call.broad_gwa
        ch_chrom_numbers = gwa_mapping_call.chromosome_numbers

        if (params.mapping && params.finemapping) {
            fine_mapping_call = FINE_MAPPING (
                ch_strains,
                ch_imputed_vcf,
                ch_broad_gwa,
                ch_chrom_numbers,
                ch_annotation,
            )
            ch_finemap_gwa = fine_mapping_call.finemap_gwa
            ch_roi_gt      = fine_mapping_call.roi_genotype_matrix
        } else {
            ch_finemap_gwa = channel.empty( )
            ch_roi_gt      = channel.empty( )
        }

        if (params.mapping && params.mediation) {
            ch_med_results = MEDIATION (
                ch_genotype_matrix,
                ch_broad_gwa,
                ch_eqtl
            )
        } else {
            ch_med_results      = channel.empty( )
        }

        if (params.mapping && params.skip_report == false) {
            ch_gwa_report = GWA_REPORTING (
                ch_strain_issues,
                ch_genotype_matrix,
                ch_broad_gwa,
                ch_finemap_gwa,
                ch_roi_gt,
                ch_med_results,
                ch_haplotypes
                )
        } else {
            ch_gwa_report = channel.empty( )
        }   
    }

    // Compile versions of tools used in the workflow
    ch_versions = channel.topic("versions")
        .map { row -> new VersionRecord(row) }
        .map { row -> "${row.tool}\t${row.version}" }
        .collect ( )

    ch_collect_versions = COLLECT_VERSIONS(
        ch_versions
    )


    emit:
    broad_gwa: Channel<BroadRecord>  = ch_broad_gwa
    fine_gwa: Channel<FineMapRecord> = ch_finemap_gwa
    med_results: Channel<MediationRecord> = ch_med_results
    gwa_report: Channel<ReportRecord> = ch_gwa_report
    strain_issues: Channel<Path> = ch_strain_issues
    versions: Value<Path> = ch_collect_versions
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
