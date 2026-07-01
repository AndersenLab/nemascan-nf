nextflow.enable.types = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GWA_REPORT                       } from '../../../modules/local/gwa_report'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record BroadGwaRecord {
    meta: Map
    trait: Path
    gwa: Path
    qtl: Path
    h2: Path
}
record FineMapGwaRecord {
    meta: Map
    gwa: Path
}
record RoiRecord {
    meta: Map
    matrix: Path
}
record MediationRecord {
    meta: Map
    medsingle: Path
}
record ConfigRecord {
    trait_name: String
    config: String
    paths: List<Path>
    template: Path
}
record GitRecord {
    repo: String
    branch: String
    commit: String
}
record ReportRecord {
    trait_name: String
    report: Path
}

workflow GWA_REPORTING {

    take:
    ch_strain_issues: Channel<Path>
    ch_genotype_matrix: Channel<Path>
    ch_broad_gwa: Channel<BroadGwaRecord>
    ch_finemap_gwa: Channel<FineMapGwaRecord>
    ch_roi_gt: Channel<RoiRecord>
    ch_med_results: Channel<MediationRecord>
    ch_haplotypes: Channel<Path>
    ch_git: Channel<GitRecord>

    main:
    ch_trait_names = ch_broad_gwa
        .map { row -> row.meta.trait_name }
        .unique ( )

    if (params.genes != null) {
        ch_genes = channel.fromPath( params.genes, checkIfExists:true )
    } else {
        ch_genes = channel.empty()
    }

    if (params.git_head != null) {
        ch_git_head = channel.fromPath( params.git_head, checkIfExists:true )
    } else if (workflow.repository != null ) {
        ch_git_params = channel.of( record(
            repository: workflow.repository,
            revision: workflow.revision,
            commitId: workflow.commitId
         ))
        //  ch_git_head = GIT2FILE( ch_git_params ).out
    } else {
        ch_git_head = channel.empty()
    }

    // Assemble config file for gwa report
    ch_config_collection = channel.empty ( )
        .mix (
            ch_broad_gwa
                .map { row -> tuple(row.meta.trait_name, record(
                    entry: "trait\t${row.trait.getName()}",
                    files: [row.trait]
                )) }
                .unique ( )
        )
        .mix (
            ch_broad_gwa
                .map { row -> tuple(row.meta.trait_name, record(
                    entry: "broad_gwa\t${row.meta.method}\t${row.gwa.getName()}",
                    files: [row.gwa]
                )) }
        )
        .mix (
            ch_broad_gwa
                .map { row -> tuple(row.meta.trait_name, record(
                    entry: "broad_qtl\t${row.meta.method}\t${row.qtl.getName()}",
                    files: [row.qtl]
                )) }
        )
        .mix (
            ch_broad_gwa
                .map { row -> tuple(row.meta.trait_name, record(
                    entry: "narrow_h2\t${row.h2.getName()}",
                    files: [row.h2]
                )) }
                .unique ( )
        )
        .mix (
            ch_broad_gwa
                .map { row -> tuple(row.meta.trait_name, record(
                    entry: "independent_tests\t${row.meta.independent_tests}",
                    files: []
                )) }
                .unique ( )
        )
        .mix (
            ch_finemap_gwa
                .map { row -> tuple(row.meta.trait_name, record(
                    entry: "finemap_gwa\t${row.meta.method}\t${row.meta.chrom}\t${row.meta.start}\t${row.meta.end}\t${row.meta.peak}\t${row.gwa.getName()}",
                    files: [row.gwa]
                )) }
        )
        .mix (
            ch_roi_gt
                .map { row -> tuple(row.meta.trait_name, record(
                    entry: "finemap_matrix\t${row.meta.chrom}\t${row.meta.start}\t${row.meta.end}\t${row.meta.peak}\t${row.matrix.getName()}",
                    files: [row.matrix]
                )) }
                .unique ( )
        )
        .mix (
            ch_med_results
                .map { row -> tuple(row.meta.trait_name, record(
                    entry: "mediation\t${row.meta.method}\t${row.meta.chrom}\t${row.meta.start}\t${row.meta.end}\t${row.meta.peak}\t${row.medsingle.getName()}",
                    files: [row.medsingle]
                )) }
        )
        .mix (
            ch_trait_names
                .combine (
                    ch_strain_issues
                        .map { row -> record(entry: "issues\t${row.getName()}", files: [row]) }
                        .mix (
                            ch_genotype_matrix
                                .map { row -> record(entry: "genotype_matrix\t${row.getName()}", files: [row]) }
                        )
                        .mix (
                            ch_haplotypes
                                .map { row -> record(entry: "haplotypes\t${row.getName()}", files: [row]) }
                        )
                        .mix (
                            channel.of( params.highlight_strains )
                                .map { row -> record(entry: "highlight_strains\t${row}", files: []) }
                        )
                        .mix (
                            channel.of ( record(
                                entry: "params\talpha=${params.alpha}\tsignificance_threshold=${params.significance_threshold}\tuser=${workflow.userName}\tworkflow_branch=${workflow.revision}\tworkflow_commit=${workflow.commitId}",
                                files: []) )
                        )
                        .mix (
                            ch_genes
                                .map { row -> record(entry: "genes\t${row.getName()}", files: [row]) }
                        )
                        .mix (
                            ch_git
                                .map { row -> record(entry: "params\trepo=${row.repo}\tbranch=${row.branch}\tcommit=${row.commit}", files: []) }
                        )
                )
                .map { row -> tuple(row[0], row[1]) }
        )
        .groupBy ( )

    // Reformat config collections for writing config files
    ch_config = ch_config_collection
        .map { row -> 
            def entries = []
            def files = []
            row[1].each { R ->
                entries.add(R.entry)
                files.addAll(R.files)
            }
            record(trait_name: row[0], config: entries.join("\n"), paths: files, template: file("${projectDir}/assets/report_template.html"))
        }
        .map { row -> new ConfigRecord(row) }

    ch_gwa_report = GWA_REPORT (
        ch_config
    )

    emit:
    report: Channel<ReportRecord> = ch_gwa_report
}
