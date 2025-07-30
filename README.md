# ![folkehelseinstituttet/hcv_illumina](docs/images/logo-engelsk-hele-navnet.jpg#gh-light-mode-only) ![folkehelseinstituttet/hcv_illumina](docs/images/logo-engelsk-hele-navnet-hvit.png#gh-dark-mode-only)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

## About HCV Illumina

**folkehelseinstituttet/hcv_illumina** is a bioinformatics pipeline used at the [Norwegian Institute of Public Health](https://www.fhi.no/en/) that is designed for highly variable viruses, and viruses that are likely to appear as co-infections between multiple strains, such as Hepatitis C Virus. The pipeline will identify the most likely major and minor strain in a sample sequenced with the Illumina platform. It will map the reads to these references using [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and create consensus sequences. For Hepatitis C Viruses the pipeline can also run a [GLUE-analysis](http://hcv-glue.cvr.gla.ac.uk/#/home) to identify drug resistance mutations.
maps Illumina reads to a reference genome and creates a consensus sequence.

## Requirements
The pipeline only requires [Nextflow](https://nextflow.io/) and [Docker](https://www.docker.com/) in order to run. Note that you must be able to run Docker as a non-root user as described [here](https://docs.docker.com/engine/install/linux-postinstall/#manage-docker-as-a-non-root-user).

## Test the pipeline
To run a minimal test:
```
nextflow run folkehelseinstituttet/hcv_illumina -profile docker,test
```
This is only to see if you can get the pipeline up and running and will not run the entire pipeline such as HCV-GLUE. The results will be in a directory called `minimal_test`.

To run a full test on a real dataset type:
```
nextflow run folkehelseinstituttet/hcv_illumina -profile docker,test_full
```
This will download a HCV Illumina dataset from SRA and run the entire pipeline. The results will be in a directory called `full_test`.

## Required parameters
### Samplesheet input
You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown below. The sample names can contain numbers and underscores (_), but not spaces, dots (.) or other symbols. The fastq_1 and fastq_2 columns must contain the full path to the gzipped paired fastq files corresponding to the same sample.

```
sample,fastq_1,fastq_2
Sample_1,/path/to/sample1_fastq_R1.fastq.gz,/path/to/sample1_fastq_R2.fastq.gz
Sample_2,/path/to/sample2_fastq_R1.fastq.gz,/path/to/sample2_fastq_R2.fastq.gz
```

The samplesheet is input to the pipeline using the `--input` parameter, e.g.:
`--input assets/samplesheet_illumina.csv`

An [example samplesheet](assets/samplesheet_illumina.csv) has been provided with the pipeline in the assets directory.

### Output directory
The output directory is specified using the `--outdir` parameter, e.g.:
`--outdir results`

### Profiles
The pipeline can be run using different profiles, which will determine how the pipeline is executed. The default profile is `docker`, which uses Docker containers to run the pipeline. You can also use `singularity` or `conda` profiles if you prefer those environments. To set the profile use the `-profile` parameter, e.g.: `-profile docker/singularity/conda`

> [!IMPORTANT]
> HCV-GLUE is only available with the Docker profile. We recommend that you always run the pipeline with Docker.

## Optional parameters


## Usage

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  folkehelseinstituttet/hcv_illumina for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

