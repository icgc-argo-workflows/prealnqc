## Introduction

**icgc-argo-workflows/prealnqc** is a reproducible bioinformatics best-practice analysis pipeline of ICGC ARGO Pre Alignment QC Workflow for DNA/RNA Sequencing Reads.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, many processes have been installed from [nf-core/modules](https://github.com/nf-core/modules). Specifically, ICGC ARGO specific modules have been submitted and installed form [icgc-argo-workflows/argo-modules](https://github.com/icgc-argo-workflows/argo-modules), in order to make them available to all ICGC ARGO pipelines!


## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`)

2. Install [`Docker`](https://docs.docker.com/engine/installation/).

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```bash
   nextflow run icgc-argo-workflows/prealnqc -profile test,standard
   ```

4. Start running your own analysis!
   ```bash
   nextflow run icgc-argo-workflows/prealnqc --input samplesheet.csv --outdir <OUTDIR> -profile standard
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
*Local mode*
- input

*RDPC mode*
- study_id
- analysis_id

### Outputs
- Sequencing quality control metrics (FastQC)
- Sequencing adapter removal, adapter contamination metrics (Cutadapt)
- Overall pipeline statistics summaries (MultiQC)

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
