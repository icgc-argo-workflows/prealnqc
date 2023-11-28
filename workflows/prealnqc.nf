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
include { INPUT_CHECK       } from '../subworkflows/local/input_check'

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
    
    // Read in samplesheet, validate and stage input files
    if (params.local_mode) {
      if (params.input) {
        ch_input = Channel.fromPath(params.input)
        ch_input_sample = INPUT_CHECK (ch_input).reads
      } 
      else { exit 1, 'Input samplesheet must be specified for local mode!' }
    } else if (params.study_id && params.analysis_ids) {
      ch_study = Channel.of(params.study_id)
      ch_analysis_ids = Channel.fromList(params.analysis_ids.split(',') as List)
      ch_input = ch_study.combine(ch_analysis_ids)

      STAGE_INPUT(ch_input)
      ch_input_sample = STAGE_INPUT.out.sample_files
      ch_metadata = STAGE_INPUT.out.meta_analysis
      ch_versions = ch_versions.mix(STAGE_INPUT.out.versions)
    
    } else { exit 1, 'study_id & analysis_ids must be specified for rdpc mode!' }


    // MODULE: Run FastQC
    FASTQC( ch_input_sample )
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    // MODULE: Perform cutadpat
    CUTADAPT( ch_input_sample )
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
    .map { meta, files -> [[id: meta.sample], files] }
    .groupTuple()
    .set{ ch_meta_qcfiles }

    // Parse the multiqc data
    PREP_METRICS (ch_meta_qcfiles, MULTIQC.out.data.collect())

    // Collect Software Versions
    CUSTOM_DUMPSOFTWAREVERSIONS (ch_versions.unique{ it.text }.collectFile(name: 'collated_versions.yml'))

    // upload QC files and metadata to song/score
    if (!params.local_mode) {
      // make metadata and files match  
      ch_metadata.map { meta, metadata -> [[id: meta.sample], metadata]}
          .unique().set{ ch_meta_metadata }

      ch_meta_metadata.join(ch_meta_qcfiles).join(PREP_METRICS.out.metrics_json)
      .set { ch_metadata_upload }

      // // generate payload
      PAYLOAD_QCMETRICS(
        ch_metadata_upload, CUSTOM_DUMPSOFTWAREVERSIONS.out.yml.collect()) 

      // SONG_SCORE_UPLOAD(PAYLOAD_QCMETRICS.out.payload_files)

      // // cleanup
      // // Gather files to remove   
      // ch_files = Channel.empty()
      // ch_files = ch_files.mix(STAGE_INPUT.out.sample_files)
      // ch_files = ch_files.mix(STAGE_INPUT.out.analysis_meta)
      // ch_files = ch_files.mix(FASTQC.out.zip)
      // ch_files = ch_files.mix(FASTQC.out.html)
      // ch_files = ch_files.mix(CUTADAPT.out.log)
      // ch_files = ch_files.mix(CUTADAPT.out.reads)
      // ch_files.map{ meta, files -> files}
      // .unique()
      // .set { ch_files_to_remove1 }

      // PAYLOAD_QCMETRICS.out.payload_files
      // .map {meta, payload, files -> files}
      // .unique()
      // .set { ch_files_to_remove2 }

      // ch_files_to_remove = Channel.empty()
      // ch_files_to_remove = ch_files_to_remove.mix(STAGE_INPUT.out.input_files)
      // ch_files_to_remove = ch_files_to_remove.mix(MULTIQC.out.report)
      // ch_files_to_remove = ch_files_to_remove.mix(MULTIQC.out.data)
      // ch_files_to_remove = ch_files_to_remove.mix(ch_files_to_remove1)
      // ch_files_to_remove = ch_files_to_remove.mix(ch_files_to_remove2)
      // CLEANUP(ch_files_to_remove.unique().collect(), SONG_SCORE_UPLOAD.out.analysis_id)
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
