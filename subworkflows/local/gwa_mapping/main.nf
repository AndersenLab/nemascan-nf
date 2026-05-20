nextflow.enable.types = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CHROMOSOME_NUMBERING                          } from '../../../modules/local/chromosome_numbering'
include { ANNOTATE_CHROMS                               } from '../../../modules/local/annotate_chroms'
include { VCF_TO_PLINK                                  } from '../../../modules/local/vcf_to_plink'
include { GCTA64_MAKE_GRM                               } from '../../../modules/local/gcta64/make_grm'
include { GCTA64_PCA                                    } from '../../../modules/local/gcta64/pca'
include { GCTA64_GWAS                                   } from '../../../modules/local/gcta64/gwas'
include { GENOTYPE_MATRIX_EIGEN                         } from '../../../modules/local/genotype_matrix_eigen'
include { AGGREGATE_MAPPINGS                            } from '../../../modules/local/aggregate_mappings'
include { NARROW_H2                                     } from '../../../modules/local/narrow_h2'

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
record ChromNumRecord {
    chromosome: String
    chrom_number: Integer
}
record EigenRecord {
    eigen: String
}
record BroadRecord {
    meta: Map
    trait: Path
    gwa: Path
    qtl: Path
    matrix: Path
    h2: Path
    independent_tests: BigDecimal
}

workflow GWA_MAPPING {

    take:
    ch_filtered_traits: Channel<TraitRecord>
    ch_pruned_vcf: Channel<VcfRecord>
    ch_genotype_matrix: Channel<Path>
    ch_gwa_method: Channel<String>

    main:
    // Convert chromosome names to numbers for GCTA compatibility
    ch_chromosome_numbering = CHROMOSOME_NUMBERING (
        ch_pruned_vcf
    )

    ch_numbered_vcf = ANNOTATE_CHROMS (
        ch_pruned_vcf
            .combine ( ch_chromosome_numbering )
            .map { row -> row[0] + record(chrom_numbering: row[1]) }
    )

    // Convert VCF to plink binary format
    ch_plink = VCF_TO_PLINK (
        ch_numbered_vcf
    )
        .combine ( ch_gwa_method )
        .map { row -> row[0] + record(meta: row[0].meta + [id: row[1], method: row[1], true_method: row[1]]) }

    // Create genotype relatedness matrix and sparse GRM
    ch_gcta64_make_grm = GCTA64_MAKE_GRM (
        ch_plink,
        params.maf,
        params.sparse_cut
    )

    // If asked for, create PCA-based covariate matrix from GRM
    if (params.pca) {
        ch_gcta64_grm = GCTA64_PCA (
            ch_gcta64_make_grm
        )
    } else {
        ch_gcta64_grm = ch_gcta64_make_grm
        .map { row -> row + record(pca: null) }
    }

    // Combine with each trait from the split traits channel
    ch_plink_trait = ch_gcta64_grm
        .combine ( ch_filtered_traits )
        .map { row -> row[0] + row[1] + record(meta: row[0].meta + row[1].meta + [id: "${row[1].meta.trait_name}_${row[0].meta.method}"]) }

    ch_gcta64_gwas = GCTA64_GWAS (
        ch_plink_trait
    )

    // Extract chromosome names
    ch_genotype_matrix_eigen_input = ch_chromosome_numbering
        .flatMap { row -> row.splitCsv ( sep:"\t", header:true ) }
        .map { row -> new ChromNumRecord(row) }
        .combine (
            ch_genotype_matrix
        )
        .map { row -> record(        
            chromosome: row[0].chromosome,
            chrom_number: row[0].chrom_number,
            genotype_matrix: row[1]
        ) }

    ch_genotype_matrix_eigen = GENOTYPE_MATRIX_EIGEN (
        ch_genotype_matrix_eigen_input
    )

    ch_independent_tests = ch_genotype_matrix_eigen
        .flatMap { row -> row.splitCsv( sep:"\t", header:true ) }
        .map { row -> new EigenRecord(row) }
        .map { row -> row.eigen.toBigDecimal() }
        .reduce { a, b -> a + b }

    // Pull out significant QTL regions
    ch_aggregate_mappings = AGGREGATE_MAPPINGS (
        ch_gcta64_gwas
            .combine( ch_chromosome_numbering )
            .map { row -> record(
                meta: row[0].meta,
                trait: row[0].trait,
                gwa: row[0].gwa,
                chrom_numbering: row[1]
            ) },
        ch_independent_tests,
        params.snp_grouping,
        params.ci_size,
        params.significance_threshold
    )

    // Calculate the narrow-sense heritability of the trait
    ch_narrow_h2 = NARROW_H2 (
        ch_filtered_traits
            .combine ( ch_genotype_matrix )
            .map { row -> row[0] + record(genotype_matrix: row[1]) }
    )

    // Join narrow h2 estimates with aggregated mappings for output
    ch_broad_mappings = ch_aggregate_mappings
        .map { row -> record(key: record(trait_name: row.meta.trait_name, method: row.meta.method), gwa_record:row) }
        .join ( ch_narrow_h2
            .combine ( ch_gwa_method )
            .map { row -> record(key: record(trait_name: row[0].meta.trait_name, method: row[1]), h2_record:row[0]) },
            by: 'key'
        )
        .combine ( ch_independent_tests )
        .map { row -> row[0].h2_record + row[0].gwa_record + record(meta: row[0].h2_record.meta + row[0].gwa_record.meta + [independent_tests: row[1]]) }

    emit:
    broad_gwa: Channel<BroadRecord>   = ch_broad_mappings
    chromosome_numbers: Channel<Path> = ch_chromosome_numbering
}
