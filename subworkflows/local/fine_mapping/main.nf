nextflow.enable.types = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { VCF_FM_FILTER             } from '../../../modules/local/vcf_fm_filter'
include { VCF_TO_FM_PLINK           } from '../../../modules/local/vcf_to_fm_plink'
include { GCTA64_MAKE_GRM           } from '../../../modules/local/gcta64/make_grm'
include { GCTA64_PCA                } from '../../../modules/local/gcta64/pca'
include { GCTA64_GWAS               } from '../../../modules/local/gcta64/gwas'
include { MAKE_FM_GENOTYPE_MATRIX   } from '../../../modules/local/make_fm_genotype_matrix'
include { FM_LD                     } from '../../../modules/local/fm_ld'
include { APPEND_ANNOTATIONS        } from '../../../modules/local/append_annotations'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


record BroadRecord {
    meta: Map
    trait: Path
    gwa: Path
    qtl: Path
    matrix: Path
    h2: Path
    independent_tests: BigDecimal
}
record VcfRecord {
    meta: Map
    vcf: Path
    vcf_index: Path
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
record FineMapRecord {
    meta: Map
    gwa: Path
}

workflow FINE_MAPPING {

    take:
    ch_strains: Channel<Path>
    ch_imputed_vcf: Channel<VcfRecord>
    ch_broad_gwa: Channel<BroadRecord>
    ch_chrom_numbers: Channel<Path>
    ch_annotations: Channel<Path>

    main:
    // Update metadata values and create region of interest files 
    ch_qtl = ch_broad_gwa
        .map { row -> row.qtl }
        .flatMap { csv -> csv.splitCsv ( sep:"\t", header:true ) }
        .map{ row -> new QtlRecord(row) }
        .map { row -> record(
            meta: [
                id: "${row.trait}_${row.method}_${row.CHROM}_${row.startPOS}_${row.endPOS}",
                trait_name:row.trait,
                marker: row.marker,
                log10p: row.log10p.toDouble(),
                method:"finemap",
                true_method:row.method,
                chrom:row.CHROM,
                start:row.startPOS.toInteger(),
                peak:row.peakPOS.toInteger(),
                end:row.endPOS.toInteger()
            ])
        }
        .combine ( ch_imputed_vcf )
        .combine ( ch_strains )
        .combine ( ch_chrom_numbers)
        .map { row -> row[1] + row[0] + record(strains: row[2], chrom_numbers: row[3]) }


    // Filter imputed VCF by strains and ROIs, followed by missing GTs and SNPs and converting chrom names to numbers
    ch_vcf_fm_filter = VCF_FM_FILTER(
        ch_qtl
    )
        .combine ( ch_chrom_numbers)
        .map { row -> row[0] + record(chrom_numbers: row[1]) }

    // Convert VCF to plink binary format
    ch_vcf_to_fm_plink = VCF_TO_FM_PLINK (
        ch_vcf_fm_filter,
        params.maf
    )

    // Create genotype relatedness matrix and sparse GRM
    ch_gcta64_make_grm = GCTA64_MAKE_GRM (
        ch_vcf_to_fm_plink,
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

    // Combine with correct trait from the split traits channel
    ch_plink_trait = ch_gcta64_grm
        .map { row -> record(key: row.meta.trait_name, grm_record: row) }
        .join (
            ch_broad_gwa
                .map { row -> record(key: row.meta.trait_name, gwa_record: record(trait: row.trait)) },
            by: "key"
        )
        .map { row ->
            row.grm_record + record( meta: row.grm_record.meta + [trait_name: row.key]) + row.gwa_record
        }

    ch_gcta64_gwas = GCTA64_GWAS (
        ch_plink_trait
    )
        .map { row -> record(meta: [
                id: row.meta.id,
                trait_name:row.meta.trait_name,
                method:row.meta.true_method,
                chrom:row.meta.chrom,
                start:row.meta.start,
                peak:row.meta.peak,
                end:row.meta.end
            ],
            gwa: row.gwa) }

    // Create genotype matrix for ROI
    ch_make_fm_genotype_matrix = MAKE_FM_GENOTYPE_MATRIX (
        ch_vcf_fm_filter
    )

    // Create LD table for ROI
    ch_fm_ld = FM_LD (
        ch_vcf_fm_filter
    )

    // Add variant annotations
    ch_ld_gwa = ch_fm_ld
        .map { row -> record(key: row.meta.id, ld_record: row) }
        .join (
            ch_gcta64_gwas
                .map { row -> record(key: row.meta.id, gwa_record: row) },
            by: "key"
        )
        .combine ( ch_annotations )
        .map { row -> row[0].gwa_record + record(ld: row[0].ld_record.ld, annotations: row[1]) }

    ch_append_annotations = APPEND_ANNOTATIONS (
        ch_ld_gwa
    )

    emit:
    finemap_gwa: Channel<FineMapRecord> = ch_append_annotations
    roi_genotype_matrix: Channel<Path>  = ch_make_fm_genotype_matrix
}
