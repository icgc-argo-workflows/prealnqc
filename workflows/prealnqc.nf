/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
// WorkflowPrealnqc.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// include { INPUT_CHECK       } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT ARGO MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SONG_SCORE_UPLOAD } from '../subworkflows/icgc-argo-workflows/song_score_upload/main' 
include { STAGE_INPUT       } from '../subworkflows/icgc-argo-workflows/stage_input/main'
include { CLEANUP           } from '../modules/icgc-argo-workflows/cleanup/main'
include { PREP_METRICS      } from '../modules/icgc-argo-workflows/prep/metrics/main'
include { PAYLOAD_QCMETRICS } from '../modules/icgc-argo-workflows/payload/qcmetrics/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { CUTADAPT                    } from '../modules/nf-core/cutadapt/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREALNQC {

    ch_versions = Channel.empty()
    
    // Stage input files
    STAGE_INPUT(params.study_id, params.analysis_ids, params.input)
    ch_versions = ch_versions.mix(STAGE_INPUT.out.versions)

    // MODULE: Run FastQC
    FASTQC( STAGE_INPUT.out.meta_files )
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    // // MODULE: Perform cutadpat
    CUTADAPT( STAGE_INPUT.out.meta_files )
    ch_versions = ch_versions.mix(CUTADAPT.out.versions)

    // Gather QC files    
    ch_qc_files = Channel.empty()
    ch_qc_files = ch_qc_files.mix(FASTQC.out.zip)
    ch_qc_files = ch_qc_files.mix(FASTQC.out.html)
    ch_qc_files = ch_qc_files.mix(CUTADAPT.out.log)
    ch_qc_files.collect{ meta, qc_files -> qc_files}
    .set{ ch_multiqc_files }

    // Perform MultiQC
    MULTIQC (
        ch_multiqc_files,
        ch_multiqc_config,
        ch_multiqc_custom_config.collect().ifEmpty([]),
        ch_multiqc_logo.collect().ifEmpty([])
    )
    ch_versions = ch_versions.mix(MULTIQC.out.versions)
    
    // Group the QC files by sampleId
    ch_qc_files
    .transpose()
    .map { meta, files -> [[id: meta.sample, study_id: meta.study_id], files] }
    .groupTuple()
    .set{ ch_meta_qcfiles }

    // Parse the multiqc data
    PREP_METRICS (ch_meta_qcfiles, MULTIQC.out.data.collect())

    // Collect Software Versions
    CUSTOM_DUMPSOFTWAREVERSIONS (ch_versions.unique{ it.text }.collectFile(name: 'collated_versions.yml'))
    
    // Combine channels to determine upload status and payload creation
    // make metadata and files match  
    STAGE_INPUT.out.meta_analysis.map { meta, metadata -> [[id: meta.sample, study_id: meta.study_id], metadata]}
        .unique().set{ ch_meta_metadata }
  
    ch_meta_metadata.join(ch_meta_qcfiles).join(PREP_METRICS.out.metrics_json)
    .set { ch_metadata_files }

    STAGE_INPUT.out.upRdpc.combine(ch_metadata_files)
    .map{upRdpc, meta, metadata, files, metrics -> 
    [[id: meta.id, study_id: meta.study_id, upRdpc: upRdpc],
      metadata, files, metrics]}
    .branch{
      upload: it[0].upRdpc
    }.set{ch_metadata_files_status}

    // generate payload
    PAYLOAD_QCMETRICS(
        ch_metadata_files_status.upload, CUSTOM_DUMPSOFTWAREVERSIONS.out.yml.collect()) 

    SONG_SCORE_UPLOAD(PAYLOAD_QCMETRICS.out.payload_files)

    if (params.cleanup) {
      // cleanup
      // Gather files to remove   
      ch_files = Channel.empty()
      ch_files = ch_files.mix(STAGE_INPUT.out.meta_files)
      ch_files = ch_files.mix(STAGE_INPUT.out.meta_analysis)
      ch_files = ch_files.mix(FASTQC.out.zip)
      ch_files = ch_files.mix(FASTQC.out.html)
      ch_files = ch_files.mix(CUTADAPT.out.log)
      ch_files = ch_files.mix(CUTADAPT.out.reads)
      ch_files.map{ meta, files -> files}
      .unique()
      .set { ch_files_to_remove1 }

      PAYLOAD_QCMETRICS.out.payload_files
      .map {meta, payload, files -> files}
      .unique()
      .set { ch_files_to_remove2 }

      ch_files_to_remove = Channel.empty()
      ch_files_to_remove = ch_files_to_remove.mix(MULTIQC.out.report)
      ch_files_to_remove = ch_files_to_remove.mix(MULTIQC.out.data)
      ch_files_to_remove = ch_files_to_remove.mix(ch_files_to_remove1)
      ch_files_to_remove = ch_files_to_remove.mix(ch_files_to_remove2)
      CLEANUP(ch_files_to_remove.unique().collect(), SONG_SCORE_UPLOAD.out.analysis_id)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// workflow.onComplete {
//     if (params.email || params.email_on_fail) {
//         NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
//     }
//     NfcoreTemplate.summary(workflow, params, log)
//     if (params.hook_url) {
//         NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
//     }
// }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
