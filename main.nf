#!/usr/bin/env nextflow

nextflow.enable.types = true

// Global default params, used in configs
params {

    // TODO nf-core: Specify your pipeline's command line flags
    // Input options
    traits: String                 = "${projectDir}/tests/data/debug/ben1_pheno.tsv"
    species: String                = null
    vcf: String                    = null
    imputed: String                = null
    annotation: String             = null
    haplotypes: String             = null
    isogroups: String              = null
    mapping: Boolean               = true
    finemapping: Boolean           = true
    mediation: Boolean             = false
    transcript_eqtl: String        = null
    transcript_expression: String  = null
    genes: String                  = null
    skip_pruning: Boolean          = false
    summarization_method: String   = "median"
    maf: BigDecimal                = 0.05
    pca: Boolean                   = false
    snp_grouping: Integer          = 1000
    ci_size: Integer               = 150
    significance_threshold: String = "BF"
    sparse_cut: BigDecimal         = 0.05
    alpha: BigDecimal              = 0.05
    gwa_method: String             = "both"
    highlight_strains: String      = ""
    skip_report: Boolean           = false

    // Boilerplate options
    outdir: Path                   = null
    publish_dir_mode: String       = 'copy'
    monochrome_logs: Boolean       = false
    version: Boolean               = false

    // Config options
    // config_profile_name: String    = null
    // config_profile_description: String = null

    // custom_config_version: String  = 'master'
    // custom_config_base: Path       = "https://raw.githubusercontent.com/nf-core/configs/${params.custom_config_version}"
    // config_profile_contact: String = null
    // config_profile_url: Path       = null
    validate_params: Boolean       = true
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    andersenlab/nemascan-nf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/andersenlab/nemascan-nf
----------------------------------------------------------------------------------------
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_nemascan-nf_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_nemascan-nf_pipeline'
include { NEMASCAN                } from './workflows/nemascan'


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

workflow {

    main:
    if (!nextflow.version.matches('>=26.04.0')) {
        error("This pipeline requires Nextflow version 26.0 or higher, but you are running version ${nextflow.version}")
    }

    pipeline_initialisation_call = PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        workflow.outputDir
    )

    //
    // WORKFLOW: Run main workflow
    //
    nemascan_call = NEMASCAN (
        pipeline_initialisation_call.traits,
        pipeline_initialisation_call.vcf,
        pipeline_initialisation_call.imputed_vcf,
        pipeline_initialisation_call.annotations,
        pipeline_initialisation_call.haplotypes,
        pipeline_initialisation_call.isogroups,
        pipeline_initialisation_call.eqtl,
        pipeline_initialisation_call.gwa_method
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        workflow.outputDir,
        params.monochrome_logs,
    )

    publish:
    // broad_gwa     = channel.empty()
    // fine_gwa      = channel.empty()
    // med_results   = channel.empty()
    // gwa_report    = channel.empty()
    // strain_issues = channel.empty()
    // versions      = channel.empty()
    broad_gwa     = nemascan_call.broad_gwa
    fine_gwa      = nemascan_call.fine_gwa
    med_results   = nemascan_call.med_results
    gwa_report    = nemascan_call.gwa_report
    strain_issues = nemascan_call.strain_issues
    versions      = nemascan_call.versions
}

output {
    broad_gwa: Channel<BroadRecord> {
        path { sample ->
            sample.gwa >> "${workflow.outputDir}/${sample.meta.trait_name}/gwa/${sample.meta.method}/"
            sample.qtl >> "${workflow.outputDir}/${sample.meta.trait_name}/gwa/${sample.meta.method}/"
            sample.h2  >> "${workflow.outputDir}/${sample.meta.trait_name}/gwa/${sample.meta.method}/"
        }
        mode params.publish_dir_mode
        overwrite true
    }
    fine_gwa: Channel<FineMapRecord> {
        path { sample -> 
            sample.gwa >> "${workflow.outputDir}/${sample.meta.trait_name}/fine_mapping/${sample.meta.method}/${sample.meta.chrom}_${sample.meta.start}_${sample.meta.end}/"
        }
        mode params.publish_dir_mode
        overwrite true
    }
    med_results: Channel<MediationRecord> {
        path { sample ->
            sample.medmulti  >> "${workflow.outputDir}/${sample.meta.trait_name}/mediation/${sample.meta.method}/${sample.meta.chrom}_${sample.meta.start}_${sample.meta.end}/"
            sample.medsingle >> "${workflow.outputDir}/${sample.meta.trait_name}/mediation/${sample.meta.method}/${sample.meta.chrom}_${sample.meta.start}_${sample.meta.end}/"
        }
        mode params.publish_dir_mode
        overwrite true
    }
    gwa_report: Channel<ReportRecord> {
        path { sample ->
            sample.report >> "${workflow.outputDir}/${sample.trait_name}/report/"
        }
        mode params.publish_dir_mode
        overwrite true
    }
    strain_issues: Channel<Path> {
        path { sample ->
            sample >> "${workflow.outputDir}/"
        }
        mode params.publish_dir_mode
        overwrite true
    }
    versions {
        path '.'
        mode params.publish_dir_mode
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
