## Introduction

**icgc-argo-workflows/prealnqc** is a reproducible bioinformatics best-practice analysis pipeline of ICGC ARGO Pre Alignment QC Workflow for DNA/RNA Sequencing Reads.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, many processes have been installed from [nf-core/modules](https://github.com/nf-core/modules). Specifically, ICGC ARGO specific modules have been submitted and installed form [icgc-argo-workflows/argo-modules](https://github.com/icgc-argo-workflows/argo-modules), in order to make them available to all ICGC ARGO pipelines!


## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`)

2. Install [`Docker`](https://docs.docker.com/engine/installation/).

3. Test the workflow running in `Local` mode on a minimal dataset with a single command:

   ```bash
   nextflow run icgc-argo-workflows/prealnqc -profile test,standard
   ```

4. Test the workflow running in `RDPC` mode with a single command if you have access to `RDPC-QA` env and have your valid api_token available:
   ```bash
   nextflow run icgc-argo-workflows/prealnqc -profile rdpc_qa,test_rdpc_qa,standard --api_token <YOUR_API_TOKEN>
   ```

5. Start running your own analysis!
   
   If you are getting the input data from & sending output data to ICGC-ARGO data center, and you have valid api_token, you can run the workflow with:
   ```bash
   nextflow run icgc-argo-workflows/prealnqc -profile <rdpc,rdpc_qa,rdpc_dev>,standard --api_token <YOUR_API_TOKEN> --study_id <STUDY_ID> --analysis_ids <ANALYSIS_IDS>
   ```
   Otherwise, you can provide the path to the input data in `samplesheet.csv` and run the workflow with:
   ```bash
   nextflow run icgc-argo-workflows/prealnqc -profile standard --input samplesheet.csv --outdir <OUTDIR>
   ```

## Pipeline summary
Depending on where the input data are coming from and output data are sending to, the workflow can be running in two modes: `Local` and `RDPC` . The major tasks performed in the workflow are:
- (`RDPC` mode only) Download input sequencing metadata/data from data center using SONG/SCORE client tools
- (`RDPC` mode only) Preprocess input sequencing reads (in FASTQ or BAM) into FASTQ file(s) per read group
- Perform FastQC analysis for FASTQ file(s) per read group
- Perform Cutadapt analysis for FASTQ file(s) per read group
- Perform MultiQC analysis to generate aggregated results
- (`RDPC` mode only) Generate SONG metadata for all collected QC metrics files and upload them to SONG/SCORE

### Inputs
#### Local mode
First, prepare a sample sheet with your input data that looks as following:

`sample_sheet.csv`:

```csv
sample,lane,fastq_1,fastq_2,single_end(optional)
TEST,C0HVY.2,C0HVY.2_r1.fq.gz,C0HVY.2_r2.fq.gz,false
TEST,D0RE2.1,D0RE2.1_r1.fq.gz,D0RE2.1_r2.fq.gz
TEST,D0RH0.2,D0RH0.2_r1.fq.gz,D0RH0.2_r2.fq.gz
```

Each row represents a read_group of sequencing reads from a sample.

Now, you can run the workflow using:

```bash
nextflow run icgc-argo-workflows/prealnqc \
   -profile <standard/singularity> \
   --local_mode true \
   --input sample_sheet.csv \
   --outdir <OUTDIR>
```

#### RDPC mode
You can run the workflow in RDPC mode by using:
```bash
nextflow run icgc-argo-workflows/prealnqc \
  -profile <rdpc,rdpc_qa,rdpc_dev>,<standard/singularity> \
  --local_mode false \
  --study_id <STUDY_ID> \
  --analysis_ids <ANALYSIS_IDS> \
  --api_token <YOUR_API_TOKEN> \ 
  --outdir <OUTDIR>
```

> **NOTE**
> Please provide workflow parameters via the CLI or Nextflow `-params-file` option. 

### Outputs
Upon completion, you can find the aggregated QC metrics under directory:
```
/path/to/outdir/prep_metrics/<sample_id>.argo_metrics.json
```

## Credits

icgc-argo-workflows/prealnqc was mostly written by Linda Xiang (@lindaxiang), with contributions from 
Andrej Benjak, Charlotte Ng, Desiree Schnidrig, Edmund Su, Miguel Vazquez, Morgan Taschuk, Raquel Manzano Garcia, Romina Royo and ICGC-ARGO Quality Control Working Group.  

Authors (alphabetical)
- Andrej Benjak
- Charlotte Ng
- Desiree Schnidrig
- Edmund Su
- Linda Xiang
- Miguel Vazquez
- Morgan Taschuk
- Raquel Manzano Garcia
- Romina Royo

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  icgc-argo-workflows/prealnqc for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
