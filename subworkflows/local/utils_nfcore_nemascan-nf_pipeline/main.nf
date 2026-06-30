nextflow.enable.types = true

//
// Subworkflow with functionality specific to the andersenlab/fq-processing-nf pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
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

workflow PIPELINE_INITIALISATION {

    take:
    version: Boolean                // boolean: Display version and exit
    validate_params: Boolean        // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs: Boolean        // boolean: Do not use coloured log outputs
    nextflow_cli_args: List<String> // array: List of positional nextflow CLI args
    outdir: Path                    // path: The output directory where the results will be saved
    

    main:

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        false,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Create channel from input file provided through params.input
    //

    if (params.mapping) {
        finemapping = params.finemapping
        mediation = params.mediation
        matrix_only = false
    } else {
        finemapping = false
        mediation = false
        matrix_only = true
    }

    //
    // Check that required file names have been defined
    //
    error_message = ''
    if (params.traits  == null) {
        error_message = error_message + "A phenotype trait file path must be passed with --traits\n"
    }
    if (params.vcf == null) {
        error_message = error_message + "A vcf file path must be passed with --vcf\n"
    }
    if (finemapping && params.imputed == null) {
        error_message = error_message + "For fine-mapping an imputed vcf file path must be passed with --imputed\n"
    }
    if (mediation && params.transcript_eqtl == null) {
        error_message = error_message + "For mediation a transcript eqtl file path must be passed with --transcript_eqtl\n"
    }
    if (mediation && params.transcript_expression == null) {
        error_message = error_message + "For mediation a transcript expression file path must be passed with --transcript_expression\n"
    }
    if (error_message.length() > 0) {
        error(error_message)
    }
    ch_traits = channel.fromPath ( params.traits, checkIfExists: true )

    ch_vcf = channel.fromPath ( params.vcf, checkIfExists: true )
        .combine (
            channel.fromPath ( params.vcf_index, checkIfExists: true )
        )
        .map { row -> record(vcf: row[0], vcf_index: row[1]) }

    if (finemapping) {
        ch_imputed_vcf = channel.fromPath ( params.imputed, checkIfExists: true )
        .combine (
            channel.fromPath ( params.imputed_index, checkIfExists: true )
        )
            .map { row -> record(vcf: row[0], vcf_index: row[1]) }
    } else {
        ch_imputed_vcf = channel.fromPath ( params.vcf, checkIfExists: true )
        .combine (
            channel.fromPath ( params.vcf_index, checkIfExists: true )
        )
            .map { row -> record(vcf: row[0], vcf_index: row[1]) }
    }

    if (matrix_only == false && params.annotation != null) {
        ch_annotations = channel.fromPath ( params.annotation, checkIfExists: true )
    } else {
        ch_annotations = channel.fromPath ( "${projectDir}/assets/NO_FILE", checkIfExists:true )
    }

    if (matrix_only == false && params.skip_report == false && params.haplotypes != null) {
        ch_haplotypes = channel.fromPath ( params.haplotypes, checkIfExists: true )
    } else {
        ch_haplotypes = channel.empty ( )
    }
    
    if (params.isogroups != null) {
        ch_isogroups = channel.fromPath ( params.isogroups, checkIfExists: true )
    } else {
        ch_isogroups = channel.fromPath ( "${projectDir}/assets/NO_FILE", checkIfExists:true )
    }

    if (mediation) {
        ch_eqtl = channel.fromPath ( params.transcript_eqtl, checkIfExists: true )
            .combine (
                channel.fromPath ( params.transcript_expression, checkIfExists: true )
            )
            .map { row -> record(eqtl: row[0], expression: row[1]) }
    } else {
        ch_eqtl = channel.empty ( )
    }

    if (params.gwa_method == "both" || params.gwa_method == "loco" || params.gwa_method == "inbred") {
        if (params.gwa_method == "both") {
            ch_gwa_method = channel.of ( "loco", "inbred" )
        } else {
            ch_gwa_method = channel.of ( params.gwa_method )
        }
    } else {
        error("The gwa_method parameter must be 'both', 'loco', or 'inbred', not '${params.gwa_method}'.")
    }

    emit:
    traits: Channel<Path>           = ch_traits
    vcf: Channel<VcfRecord>         = ch_vcf
    imputed_vcf: Channel<VcfRecord> = ch_imputed_vcf
    annotations: Channel<Path>      = ch_annotations
    haplotypes: Channel<Path>       = ch_haplotypes
    isogroups: Channel<Path>        = ch_isogroups
    eqtl: Channel<EqtlRecord>       = ch_eqtl
    gwa_method: Channel<String>     = ch_gwa_method
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    outdir: Path          //    path: Path to output directory where results will be published
    monochrome_logs: Boolean // boolean: Disable ANSI colour codes in log output

    main:

    //
    // Completion email and summary
    //
    workflow.onComplete {

        completionSummary(monochrome_logs)
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
