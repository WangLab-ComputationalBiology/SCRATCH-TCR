#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { VDJ_QC_SW }            from './subworkflows/local/vdj_qc.nf'
include { TCELL_INTEGRATION_SW } from './subworkflows/local/tcell_integration.nf'
include { TCRI_SW }              from './subworkflows/local/tcri.nf'
include { CONGA_SW }             from './subworkflows/local/conga.nf'
include { GLIPH2_SW }            from './subworkflows/local/gliph2.nf'
include { TCRDIST3_SW }          from './subworkflows/local/tcrdist3.nf'
include { GIANA_SW }             from './subworkflows/local/giana.nf'
include { CONSENSUS_SW }         from './subworkflows/local/consensus_clustering.nf'
include { REPERTOIRE_SW }        from './subworkflows/local/repertoire.nf'
include { MASTER_SUMMARY_SW }    from './subworkflows/local/master_summary.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check mandatory parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.input_annotated_object) {
    input_annotated_object = file(params.input_annotated_object)
} else {
    exit 1, 'Please provide --input_annotated_object <PATH/TO/annotated_object.RDS> !'
}

if (params.input_vdj_contigs) {
    input_vdj_contigs = params.input_vdj_contigs
} else {
    exit 1, 'Please provide --input_vdj_contigs <PATH/GLOB/TO/VDJ/outs/*> !'
}

if (params.sample_sheet) {
    sample_sheet = file(params.sample_sheet)
} else {
    exit 1, 'Please provide --sample_sheet <PATH/TO/sample_sheet.csv> !'
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    log.info """\

        Parameters:

        Annotated object: ${input_annotated_object}
        VDJ contigs:      ${params.input_vdj_contigs}
        Sample sheet:     ${sample_sheet}

    """

    def enabled = { x -> x == null || x == true }

    ch_annotated_object = Channel.fromPath(input_annotated_object, checkIfExists: true)
    ch_sample_sheet     = Channel.fromPath(sample_sheet, checkIfExists: true)
    ch_project_name     = Channel.value(params.project_name)

    /*
     * Placeholder file for optional downstream inputs when a module is disabled.
     */
    def nofile = file("${projectDir}/assets/NO_FILE")

    /*
     * Master-summary inputs for VDJ QC curated outputs
     * Default to NO_FILE when unavailable.
     */
    def vdj_qc_per_sample_compact_for_master         = Channel.value(nofile)
    def vdj_qc_before_after_summary_for_master       = Channel.value(nofile)
    def vdj_qc_sample_sheet_resolved_for_master      = Channel.value(nofile)
    def vdj_qc_clone_rank_abundance_for_master       = Channel.value(nofile)

    def vdj_qc_before_after_retention_fig_for_master = Channel.value(nofile)
    def vdj_qc_pairing_bar_fig_for_master            = Channel.value(nofile)
    def vdj_qc_clone_rank_abundance_fig_for_master   = Channel.value(nofile)
    def vdj_qc_multiple_chains_fig_for_master        = Channel.value(nofile)

    /*
     * Step 1: VDJ QC
     */
    vdj_qc_out = VDJ_QC_SW(
        ch_sample_sheet,
        ch_project_name,
        ch_annotated_object
    )

    /*
     * Curated VDJ QC files for master summary
     */
    if (vdj_qc_out?.qc_tables) {
        vdj_qc_per_sample_compact_for_master = vdj_qc_out.qc_tables
            .flatten()
            .filter { it.getName() == 'vdj_qc_per_sample_compact.tsv' }
            .ifEmpty(nofile)

        vdj_qc_before_after_summary_for_master = vdj_qc_out.qc_tables
            .flatten()
            .filter { it.getName() == 'qc_contigs_before_after_summary.tsv' }
            .ifEmpty(nofile)

        vdj_qc_sample_sheet_resolved_for_master = vdj_qc_out.qc_tables
            .flatten()
            .filter { it.getName() == 'sample_sheet_resolved.tsv' }
            .ifEmpty(nofile)

        vdj_qc_clone_rank_abundance_for_master = vdj_qc_out.qc_tables
            .flatten()
            .filter { it.getName() == 'clone_rank_abundance.tsv' }
            .ifEmpty(nofile)
    }

    if (vdj_qc_out?.qc_figures) {
        vdj_qc_before_after_retention_fig_for_master = vdj_qc_out.qc_figures
            .flatten()
            .filter { it.getName() == 'qc_before_after_retention.png' }
            .ifEmpty(nofile)

        vdj_qc_pairing_bar_fig_for_master = vdj_qc_out.qc_figures
            .flatten()
            .filter { it.getName() == 'pairing_bar_by_sample.png' }
            .ifEmpty(nofile)

        vdj_qc_clone_rank_abundance_fig_for_master = vdj_qc_out.qc_figures
            .flatten()
            .filter { it.getName() == 'clone_rank_abundance.png' }
            .ifEmpty(nofile)

        vdj_qc_multiple_chains_fig_for_master = vdj_qc_out.qc_figures
            .flatten()
            .filter { it.getName() == 'multiple_chains_by_sample.png' }
            .ifEmpty(nofile)
    }

    /*
     * Step 2: T-cell integration
     * This is the common baseline for all downstream analyses.
     */
    tcell_out = TCELL_INTEGRATION_SW(
        vdj_qc_out.contigs_after_qc,
        ch_annotated_object,
        ch_project_name
    )

    /*
     * Step 3: Independent downstream modules
     */
    tcri_done       = Channel.empty()
    conga_done      = Channel.empty()
    gliph_done      = Channel.empty()
    tcrdist3_done   = Channel.empty()
    giana_done      = Channel.empty()
    repertoire_done = Channel.empty()

    /*
     * Consensus-specific upstream files.
     */
    gliph_export_for_consensus   = Channel.value(nofile)
    tcrdist_export_for_consensus = Channel.value(nofile)
    giana_export_for_consensus   = Channel.value(nofile)

    if (enabled(params.run_tcri)) {
        tcri_out = TCRI_SW(
            tcell_out.seurat_tcells_with_tcr,
            tcell_out.export_cells,
            ch_project_name
        )
        tcri_done = tcri_out.report_html
    }

    if (enabled(params.run_conga)) {
        conga_out = CONGA_SW(
            tcell_out.seurat_tcells_with_tcr,
            tcell_out.export_cells,
            ch_project_name
        )
        conga_done = conga_out.report_html
    }

    if (enabled(params.run_gliph2)) {
        gliph_out = GLIPH2_SW(
            tcell_out.seurat_tcells_with_tcr,
            tcell_out.export_cells,
            ch_project_name
        )
        gliph_done = gliph_out.report_html
        gliph_export_for_consensus = gliph_out.export_cells
    }

    if (enabled(params.run_tcrdist3)) {
        tcrdist3_out = TCRDIST3_SW(
            tcell_out.seurat_tcells_with_tcr,
            tcell_out.export_cells,
            ch_project_name
        )
        tcrdist3_done = tcrdist3_out.report_html
        tcrdist_export_for_consensus = tcrdist3_out.export_cells
    }

    if (enabled(params.run_giana)) {
        giana_out = GIANA_SW(
            tcell_out.seurat_tcells_with_tcr,
            tcell_out.export_cells,
            ch_project_name
        )
        giana_done = giana_out.report_html
        giana_export_for_consensus = giana_out.export_cells
    }

    if (enabled(params.run_repertoire)) {
        repertoire_out = REPERTOIRE_SW(
            tcell_out.seurat_tcells_with_tcr,
            tcell_out.export_cells,
            ch_project_name
        )
        repertoire_done = repertoire_out.report_html
    }

    /*
     * Step 4: Barrier for master summary only
     */
    downstream_done = tcri_done
        .mix(conga_done)
        .mix(gliph_done)
        .mix(tcrdist3_done)
        .mix(giana_done)
        .mix(repertoire_done)

    /*
     * Step 5: Consensus clustering
     */
    if (enabled(params.run_consensus)) {
        consensus_out = CONSENSUS_SW(
            tcell_out.seurat_tcells_with_tcr,
            tcell_out.export_cells,
            gliph_export_for_consensus,
            tcrdist_export_for_consensus,
            giana_export_for_consensus,
            ch_project_name
        )
    }

    /*
     * Step 6: Master summary
     */
    if (enabled(params.run_master_summary)) {
        master_barrier = downstream_done.collect()

        master_summary_out = MASTER_SUMMARY_SW(
            tcell_out.seurat_tcells_with_tcr,
            tcell_out.export_cells,

            vdj_qc_per_sample_compact_for_master,
            vdj_qc_before_after_summary_for_master,
            vdj_qc_sample_sheet_resolved_for_master,
            vdj_qc_clone_rank_abundance_for_master,

            vdj_qc_before_after_retention_fig_for_master,
            vdj_qc_pairing_bar_fig_for_master,
            vdj_qc_clone_rank_abundance_fig_for_master,
            vdj_qc_multiple_chains_fig_for_master,

            master_barrier,
            ch_project_name
        )
    }
}

workflow.onComplete {
    log.info(
        workflow.success
            ? """
              Done! Open results in:
              - VDJ QC:              ${launchDir}/${params.outdir}/VDJ_QC
              - TCell Integration:   ${launchDir}/${params.outdir}/TCell_Integration
              - TCRi:                ${launchDir}/${params.outdir}/TCRi
              - CoNGA:               ${launchDir}/${params.outdir}/CoNGA
              - GLIPH2:              ${launchDir}/${params.outdir}/GLIPH2
              - TCRdist3:            ${launchDir}/${params.outdir}/TCRdist3
              - GIANA:               ${launchDir}/${params.outdir}/GIANA
              - Repertoire:          ${launchDir}/${params.outdir}/Repertoire
              - Consensus:           ${launchDir}/${params.outdir}/Consensus_Clustering
              - Master Summary:      ${launchDir}/${params.outdir}/Master_Summary
              """
            : "Oops... Something went wrong"
    )
}

// #!/usr/bin/env nextflow
// nextflow.enable.dsl = 2

// include { VDJ_QC_SW }            from './subworkflows/local/vdj_qc.nf'
// include { TCELL_INTEGRATION_SW } from './subworkflows/local/tcell_integration.nf'
// include { TCRI_SW }              from './subworkflows/local/tcri.nf'
// include { CONGA_SW }             from './subworkflows/local/conga.nf'
// include { GLIPH2_SW }            from './subworkflows/local/gliph2.nf'
// include { TCRDIST3_SW }          from './subworkflows/local/tcrdist3.nf'
// include { GIANA_SW }             from './subworkflows/local/giana.nf'
// include { CONSENSUS_SW }         from './subworkflows/local/consensus_clustering.nf'
// include { REPERTOIRE_SW }        from './subworkflows/local/repertoire.nf'
// include { MASTER_SUMMARY_SW }    from './subworkflows/local/master_summary.nf'

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     Check mandatory parameters
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

// if (params.input_annotated_object) {
//     input_annotated_object = file(params.input_annotated_object)
// } else {
//     exit 1, 'Please provide --input_annotated_object <PATH/TO/annotated_object.RDS> !'
// }

// if (params.input_vdj_contigs) {
//     input_vdj_contigs = params.input_vdj_contigs
// } else {
//     exit 1, 'Please provide --input_vdj_contigs <PATH/GLOB/TO/VDJ/outs/*> !'
// }

// if (params.sample_sheet) {
//     sample_sheet = file(params.sample_sheet)
// } else {
//     exit 1, 'Please provide --sample_sheet <PATH/TO/sample_sheet.csv> !'
// }

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     RUN WORKFLOW
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

// workflow {

//     log.info """\

//         Parameters:

//         Annotated object: ${input_annotated_object}
//         VDJ contigs:      ${params.input_vdj_contigs}
//         Sample sheet:     ${sample_sheet}

//     """

//     def enabled = { x -> x == null || x == true }

//     ch_annotated_object = Channel.fromPath(input_annotated_object, checkIfExists: true)
//     ch_sample_sheet     = Channel.fromPath(sample_sheet, checkIfExists: true)
//     ch_project_name     = Channel.value(params.project_name)

//     /*
//      * Placeholder file for optional downstream inputs when a module is disabled.
//      * Make sure this file exists in the repo.
//      */
//     def nofile = file("${projectDir}/assets/NO_FILE")

//     /*
//      * Master-summary inputs for VDJ QC curated outputs
//      * Default to NO_FILE when unavailable.
//      */
//     def vdj_qc_per_sample_compact_for_master         = Channel.value(nofile)
//     def vdj_qc_before_after_summary_for_master       = Channel.value(nofile)
//     def vdj_qc_sample_sheet_resolved_for_master      = Channel.value(nofile)
//     def vdj_qc_clone_rank_abundance_for_master       = Channel.value(nofile)

//     def vdj_qc_before_after_retention_fig_for_master = Channel.value(nofile)
//     def vdj_qc_pairing_bar_fig_for_master            = Channel.value(nofile)
//     def vdj_qc_clone_rank_abundance_fig_for_master   = Channel.value(nofile)
//     def vdj_qc_multiple_chains_fig_for_master        = Channel.value(nofile)

    
//     /*
//      * Step 1: VDJ QC
//      */
//     vdj_qc_out = VDJ_QC_SW(
//         ch_sample_sheet,
//         ch_project_name,
//         ch_annotated_object
//     )

//     /*
//      * Step 2: T-cell integration
//      * This is the common baseline for all downstream analyses.
//      */
//     tcell_out = TCELL_INTEGRATION_SW(
//         vdj_qc_out.contigs_after_qc,
//         ch_annotated_object,
//         ch_project_name
//     )

//     /*
//      * Step 3: Independent downstream modules
//      * All of these run in parallel after T-cell integration finishes,
//      * because they all consume the same integrated baseline outputs.
//      */

//     tcri_done       = Channel.empty()
//     conga_done      = Channel.empty()
//     gliph_done      = Channel.empty()
//     tcrdist3_done   = Channel.empty()
//     giana_done      = Channel.empty()
//     repertoire_done = Channel.empty()
    

//     /*
//      * Consensus-specific upstream files.
//      * These default to a placeholder file when the corresponding module is disabled.
//      */
//     gliph_export_for_consensus   = Channel.value(nofile)
//     tcrdist_export_for_consensus = Channel.value(nofile)
//     giana_export_for_consensus   = Channel.value(nofile)

//     if (enabled(params.run_tcri)) {
//         tcri_out = TCRI_SW(
//             tcell_out.seurat_tcells_with_tcr,
//             tcell_out.export_cells,
//             ch_project_name
//         )
//         tcri_done = tcri_out.report_html
//     }

//     if (enabled(params.run_conga)) {
//         conga_out = CONGA_SW(
//             tcell_out.seurat_tcells_with_tcr,
//             tcell_out.export_cells,
//             ch_project_name
//         )
//         conga_done = conga_out.report_html
//     }

//     if (enabled(params.run_gliph2)) {
//         gliph_out = GLIPH2_SW(
//             tcell_out.seurat_tcells_with_tcr,
//             tcell_out.export_cells,
//             ch_project_name
//         )
//         gliph_done = gliph_out.report_html
//         gliph_export_for_consensus = gliph_out.export_cells
//     }

//     if (enabled(params.run_tcrdist3)) {
//         tcrdist3_out = TCRDIST3_SW(
//             tcell_out.seurat_tcells_with_tcr,
//             tcell_out.export_cells,
//             ch_project_name
//         )
//         tcrdist3_done = tcrdist3_out.report_html
//         tcrdist_export_for_consensus = tcrdist3_out.export_cells
//     }

//     if (enabled(params.run_giana)) {
//         giana_out = GIANA_SW(
//             tcell_out.seurat_tcells_with_tcr,
//             tcell_out.export_cells,
//             ch_project_name
//         )
//         giana_done = giana_out.report_html
//         giana_export_for_consensus = giana_out.export_cells
//     }

//     if (enabled(params.run_repertoire)) {
//         repertoire_out = REPERTOIRE_SW(
//             tcell_out.seurat_tcells_with_tcr,
//             tcell_out.export_cells,
//             ch_project_name
//         )
//         repertoire_done = repertoire_out.report_html
//     }

//     /*
//      * Step 4: Barrier for master summary only
//      * Master summary is an aggregation module and should wait for all enabled
//      * downstream modules to finish.
//      */

//     downstream_done = tcri_done
//         .mix(conga_done)
//         .mix(gliph_done)
//         .mix(tcrdist3_done)
//         .mix(giana_done)
//         .mix(repertoire_done)

//     /*
//      * Step 5: Consensus clustering
//      * Consensus is a downstream aggregation module built from clustering methods.
//      * It should depend on the actual export files from GLIPH2, TCRdist3, and GIANA.
//      * Because these are real process outputs, Nextflow will wait automatically.
//      */

//     if (enabled(params.run_consensus)) {
//         consensus_out = CONSENSUS_SW(
//             tcell_out.seurat_tcells_with_tcr,
//             tcell_out.export_cells,
//             gliph_export_for_consensus,
//             tcrdist_export_for_consensus,
//             giana_export_for_consensus,
//             ch_project_name
//         )
//     }

//     /*
//      * Step 6: Master summary
//      * Should run after all selected downstream modules have finished.
//      * NOTE: This still assumes MASTER_SUMMARY_SW accepts one extra dummy barrier input.
//      */

//     if (enabled(params.run_master_summary)) {
//         master_barrier = downstream_done.collect()

//         master_summary_out = MASTER_SUMMARY_SW(
//             tcell_out.seurat_tcells_with_tcr,
//             tcell_out.export_cells,

//             vdj_qc_per_sample_compact_for_master,
//             vdj_qc_before_after_summary_for_master,
//             vdj_qc_sample_sheet_resolved_for_master,
//             vdj_qc_clone_rank_abundance_for_master,

//             vdj_qc_before_after_retention_fig_for_master,
//             vdj_qc_pairing_bar_fig_for_master,
//             vdj_qc_clone_rank_abundance_fig_for_master,
//             vdj_qc_multiple_chains_fig_for_master,

//             master_barrier,
//             ch_project_name
//         )
//     }
// }

// workflow.onComplete {
//     log.info(
//         workflow.success
//             ? """
//               Done! Open results in:
//               - VDJ QC:              ${launchDir}/${params.outdir}/VDJ_QC
//               - TCell Integration:   ${launchDir}/${params.outdir}/TCell_Integration
//               - TCRi:                ${launchDir}/${params.outdir}/TCRi
//               - CoNGA:               ${launchDir}/${params.outdir}/CoNGA
//               - GLIPH2:              ${launchDir}/${params.outdir}/GLIPH2
//               - TCRdist3:            ${launchDir}/${params.outdir}/TCRdist3
//               - GIANA:               ${launchDir}/${params.outdir}/GIANA
//               - Repertoire:          ${launchDir}/${params.outdir}/Repertoire
//               - Consensus:           ${launchDir}/${params.outdir}/Consensus_Clustering
//               - Master Summary:      ${launchDir}/${params.outdir}/Master_Summary
//               """
//             : "Oops... Something went wrong"
//     )
// }


// // #!/usr/bin/env nextflow
// // nextflow.enable.dsl = 2

// // include { VDJ_QC_SW }            from './subworkflows/local/vdj_qc.nf'
// // include { TCELL_INTEGRATION_SW } from './subworkflows/local/tcell_integration.nf'
// // include { TCRI_SW }              from './subworkflows/local/tcri.nf'
// // include { CONGA_SW }             from './subworkflows/local/conga.nf'
// // include { GLIPH2_SW }            from './subworkflows/local/gliph2.nf'
// // include { TCRDIST3_SW }          from './subworkflows/local/tcrdist3.nf'
// // include { GIANA_SW }             from './subworkflows/local/giana.nf'
// // include { CONSENSUS_SW }         from './subworkflows/local/consensus_clustering.nf'
// // include { REPERTOIRE_SW }        from './subworkflows/local/repertoire.nf'
// // include { MASTER_SUMMARY_SW }    from './subworkflows/local/master_summary.nf'

// // /*
// // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// //     Check mandatory parameters
// // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// // */

// // if (params.input_annotated_object) {
// //     input_annotated_object = file(params.input_annotated_object)
// // } else {
// //     exit 1, 'Please provide --input_annotated_object <PATH/TO/annotated_object.RDS> !'
// // }

// // if (params.input_vdj_contigs) {
// //     input_vdj_contigs = params.input_vdj_contigs
// // } else {
// //     exit 1, 'Please provide --input_vdj_contigs <PATH/GLOB/TO/VDJ/outs/*> !'
// // }

// // if (params.sample_sheet) {
// //     sample_sheet = file(params.sample_sheet)
// // } else {
// //     exit 1, 'Please provide --sample_sheet <PATH/TO/sample_sheet.csv> !'
// // }

// // /*
// // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// //     RUN WORKFLOW
// // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// // */

// // workflow {

// //     log.info """\

// //         Parameters:

// //         Annotated object: ${input_annotated_object}
// //         VDJ contigs:      ${params.input_vdj_contigs}
// //         Sample sheet:     ${sample_sheet}

// //     """

// //     def enabled = { x -> x == null || x == true }

// //     ch_annotated_object = Channel.fromPath(input_annotated_object, checkIfExists: true)
// //     ch_sample_sheet     = Channel.fromPath(sample_sheet, checkIfExists: true)
// //     ch_project_name     = Channel.value(params.project_name)

// //     /*
// //      * Step 1: VDJ QC
// //      */
// //     vdj_qc_out = VDJ_QC_SW(
// //         ch_sample_sheet,
// //         ch_project_name,
// //         ch_annotated_object
// //     )

// //     /*
// //      * Step 2: T-cell integration
// //      * This is the common baseline for all downstream analyses.
// //      */
// //     tcell_out = TCELL_INTEGRATION_SW(
// //         vdj_qc_out.contigs_after_qc,
// //         ch_annotated_object,
// //         ch_project_name
// //     )

// //     /*
// //      * Step 3: Independent downstream modules
// //      * All of these run in parallel after T-cell integration finishes,
// //      * because they all consume the same integrated baseline outputs.
// //      */

// //     tcri_done       = Channel.empty()
// //     conga_done      = Channel.empty()
// //     gliph_done      = Channel.empty()
// //     tcrdist3_done   = Channel.empty()
// //     giana_done      = Channel.empty()
// //     repertoire_done = Channel.empty()

// //     if (enabled(params.run_tcri)) {
// //         tcri_out = TCRI_SW(
// //             tcell_out.seurat_tcells_with_tcr,
// //             tcell_out.export_cells,
// //             ch_project_name
// //         )
// //         tcri_done = tcri_out.report_html
// //     }

// //     if (enabled(params.run_conga)) {
// //         conga_out = CONGA_SW(
// //             tcell_out.seurat_tcells_with_tcr,
// //             tcell_out.export_cells,
// //             ch_project_name
// //         )
// //         conga_done = conga_out.report_html
// //     }

// //     if (enabled(params.run_gliph2)) {
// //         gliph_out = GLIPH2_SW(
// //             tcell_out.seurat_tcells_with_tcr,
// //             tcell_out.export_cells,
// //             ch_project_name
// //         )
// //         gliph_done = gliph_out.report_html
// //     }

// //     if (enabled(params.run_tcrdist3)) {
// //         tcrdist3_out = TCRDIST3_SW(
// //             tcell_out.seurat_tcells_with_tcr,
// //             tcell_out.export_cells,
// //             ch_project_name
// //         )
// //         tcrdist3_done = tcrdist3_out.report_html
// //     }

// //     if (enabled(params.run_giana)) {
// //         giana_out = GIANA_SW(
// //             tcell_out.seurat_tcells_with_tcr,
// //             tcell_out.export_cells,
// //             ch_project_name
// //         )
// //         giana_done = giana_out.report_html
// //     }

// //     if (enabled(params.run_repertoire)) {
// //         repertoire_out = REPERTOIRE_SW(
// //             tcell_out.seurat_tcells_with_tcr,
// //             tcell_out.export_cells,
// //             ch_project_name
// //         )
// //         repertoire_done = repertoire_out.report_html
// //     }

// //     /*
// //      * Step 4: Barrier channels
// //      * These collect completion signals from selected downstream modules.
// //      */

// //     downstream_done = tcri_done
// //         .mix(conga_done)
// //         .mix(gliph_done)
// //         .mix(tcrdist3_done)
// //         .mix(giana_done)
// //         .mix(repertoire_done)

// //     consensus_input_done = conga_done
// //         .mix(gliph_done)
// //         .mix(tcrdist3_done)
// //         .mix(giana_done)

// //     /*
// //      * Step 5: Consensus clustering
// //      * Should run after selected clustering-capable methods finish.
// //      * NOTE: This requires CONSENSUS_SW to accept one extra dummy barrier input.
// //      */

// //     if (enabled(params.run_consensus)) {
// //         consensus_barrier = consensus_input_done.collect()

// //         consensus_out = CONSENSUS_SW(
// //             tcell_out.seurat_tcells_with_tcr,
// //             tcell_out.export_cells,
// //             consensus_barrier,
// //             ch_project_name
// //         )
// //     }

// //     /*
// //      * Step 6: Master summary
// //      * Should run after all selected downstream modules have finished.
// //      * NOTE: This requires MASTER_SUMMARY_SW to accept one extra dummy barrier input.
// //      */

// //     if (enabled(params.run_master_summary)) {
// //         master_barrier = downstream_done.collect()

// //         master_summary_out = MASTER_SUMMARY_SW(
// //             tcell_out.seurat_tcells_with_tcr,
// //             tcell_out.export_cells,
// //             master_barrier,
// //             ch_project_name
// //         )
// //     }
// // }

// // workflow.onComplete {
// //     log.info(
// //         workflow.success
// //             ? """
// //               Done! Open results in:
// //               - VDJ QC:              ${launchDir}/${params.outdir}/VDJ_QC
// //               - TCell Integration:   ${launchDir}/${params.outdir}/TCell_Integration
// //               - TCRi:                ${launchDir}/${params.outdir}/TCRi
// //               - CoNGA:               ${launchDir}/${params.outdir}/CoNGA
// //               - GLIPH2:              ${launchDir}/${params.outdir}/GLIPH2
// //               - TCRdist3:            ${launchDir}/${params.outdir}/TCRdist3
// //               - GIANA:               ${launchDir}/${params.outdir}/GIANA
// //               - Repertoire:          ${launchDir}/${params.outdir}/Repertoire
// //               - Consensus:           ${launchDir}/${params.outdir}/Consensus_Clustering
// //               - Master Summary:      ${launchDir}/${params.outdir}/Master_Summary
// //               """
// //             : "Oops... Something went wrong"
// //     )
// // }

