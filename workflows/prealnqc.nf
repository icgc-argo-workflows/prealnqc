/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowPrealnqc.initialise(params, log)

// Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
// def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
// for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if ( params.local_mode ) {
  if (params.input) {
    ch_input = file(params.input, checkIfExists: true)
  } 
  else { exit 1, 'Input samplesheet must be specified for local mode!' }
} else if (params.study_id && params.analysis_id) {
  ch_input = [params.study_id, params.analysis_id]
} else { exit 1, 'study_id & analysis_id must be specified for rdpc mode!' }

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
include { PAYLOAD_QCMETRICS } from '../modules/icgc-argo-workflows/payload/qcmetrics/main'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { SONG_SCORE_UPLOAD } from '../subworkflows/icgc-argo-workflows/song_score_upload/main' 
include { STAGE_INPUT       } from '../subworkflows/icgc-argo-workflows/stage_input/main'
include { INPUT_CHECK       } from '../subworkflows/local/input_check'
include { CLEANUP           } from '../modules/icgc-argo-workflows/cleanup/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
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
      ch_input_sample = INPUT_CHECK (ch_input).reads

    } else {
      STAGE_INPUT(ch_input)
      ch_input_sample = STAGE_INPUT.out.sample_files
      ch_metadata = STAGE_INPUT.out.analysis_meta
      ch_versions = ch_versions.mix(STAGE_INPUT.out.versions)
    
    } 


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
    
    // Match the QC files with the metadata info
    ch_qc_files
    .map { meta, files -> [[id: meta.sample, study_id: meta.study_id], files] }
    .groupTuple()
    .map { meta, files -> meta }
    .set{ ch_meta }

    ch_meta.combine(ch_multiqc_files.concat(MULTIQC.out.report, MULTIQC.out.data).collect().toList())
    .set {ch_qc_files_sample}

    // Collect Software Versions
    CUSTOM_DUMPSOFTWAREVERSIONS (ch_versions.unique{ it.text }.collectFile(name: 'collated_versions.yml'))

    // upload QC files and metadata to song/score
    if (!params.local_mode) {
      // make metadata and files match  
      ch_metadata.map { meta, metadata -> [[id: meta.sample, study_id: meta.study_id], metadata]}
          .unique().set{ ch_metadata_sample }
          
      ch_qc_files_sample.join(ch_metadata_sample)
      .set { ch_metadata_upload }

      // generate payload
      PAYLOAD_QCMETRICS(
        ch_metadata_upload, '', '', CUSTOM_DUMPSOFTWAREVERSIONS.out.yml.collect()) 

      SONG_SCORE_UPLOAD(PAYLOAD_QCMETRICS.out.payload_files)

      // cleanup
      // Gather files to remove   
      ch_files = Channel.empty()
      ch_files = ch_files.mix(STAGE_INPUT.out.sample_files)
      ch_files = ch_files.mix(STAGE_INPUT.out.analysis_meta)
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

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
