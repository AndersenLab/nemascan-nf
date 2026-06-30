nextflow.enable.types = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GET_VCF_STRAINS                            } from '../../../modules/local/get_vcf_strains'
include { CLEAN_TRAITS                               } from '../../../modules/local/clean_traits'
include { VCF_FILTER                                 } from '../../../modules/local/vcf_filter'
include { LD_PRUNED_MARKERS                          } from '../../../modules/local/ld_pruned_markers'
include { MAKE_GENOTYPE_MATRIX                       } from '../../../modules/local/make_genotype_matrix'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record VcfRecord {
    vcf: Path
    vcf_index: Path
}
record TraitRecord {
    meta: Map
    trait: Path
}

workflow GENOTYPE_MATRIX {

    take:
    ch_traits: Channel<Path>
    ch_vcf: Channel<VcfRecord>
    ch_isogroups: Channel<Path>

    main:
    ch_get_vcf_strains = GET_VCF_STRAINS ( ch_vcf )

    // Filter traits and split into invidual trait files
    ch_clean_traits_input = ch_traits
        .combine ( ch_isogroups )
        .combine( ch_get_vcf_strains )
        .combine( channel.of ( params.summarization_method ) )
        .combine( channel.of ( params.skip_pruning ) )
        .map { row -> record(
            traits: row[0],
            isogroups: row[1],
            vcf_strains: row[2],
            summarization_method: row[3],
            skip_pruning: row[4]
        ) }

    ch_clean_traits = CLEAN_TRAITS (
        ch_clean_traits_input
    )

    ch_included_strains = ch_clean_traits
        .map { row -> row.included }

    ch_pruned_traits = ch_clean_traits
        .flatMap { row -> row.traits }
        .map { row -> record(
            meta: [trait_name:row.baseName],
            trait: row
        ) }

    // Filter VCF to only phenotyped strains, no missing GTs, SNPs-only
    ch_vcf_filter = VCF_FILTER (
        ch_vcf
            .combine ( ch_included_strains )
            .map { row -> row[0] + record(samples: row[1]) }
    )

    // Find a LD-pruned set of markers
    ch_ld_pruned_markers = LD_PRUNED_MARKERS (
        ch_vcf_filter,
        params.maf
    )

    // Prune VCF file and convert to genotype matrix
    ch_make_genotype_matrix = MAKE_GENOTYPE_MATRIX (
        ch_vcf_filter
            .combine ( ch_ld_pruned_markers )
            .map { row -> row[0] + record(markers: row[1]) }
    )

    emit:
    traits: Channel<TraitRecord>   = ch_pruned_traits
    pruned_vcf: Channel<VcfRecord> = ch_make_genotype_matrix.map { row -> row.vcf }
    strains: Channel<Path>         = ch_included_strains
    strain_issues: Channel<Path>   = ch_clean_traits.map { row -> row.issues }
    markers: Channel<Path>         = ch_ld_pruned_markers
    genotype_matrix: Channel<Path> = ch_make_genotype_matrix.map { row -> row.matrix }
}
