# ![folkehelseinstituttet/hcv_illumina](docs/images/logo-engelsk-hele-navnet.jpg#gh-light-mode-only) ![folkehelseinstituttet/hcv_illumina](docs/images/logo-engelsk-hele-navnet-hvit.png#gh-dark-mode-only)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

## Table of Contents
- [About HCV Illumina](#about-hcv-illumina)
- [Requirements](#requirements)
- [Run the pipeline](#run-the-pipeline)
- [Test the pipeline](#test-the-pipeline)
- [Required parameters](#required-parameters)
  - [Samplesheet input](#samplesheet-input)
  - [Output directory](#output-directory)
  - [Profiles](#profiles)
  - [Provide parameters in a file](#provide-parameters-in-a-file)
- [Optional parameters](#optional-parameters)
  - [Kraken2 databases](#kraken2-databases)
  - [HCV reference sequences](#hcv-reference-sequences)
  - [Co-infections (major and minor strains)](#co-infections-major-and-minor-strains)
- [Starting and stopping the pipeline](#starting-and-stopping-the-pipeline)
- [Customizing the pipeline](#customizing-the-pipeline)
- [Output files](#output-files)
- [Citations](#citations)

## About HCV Illumina

**folkehelseinstituttet/hcv_illumina** is a bioinformatics pipeline used at the [Norwegian Institute of Public Health](https://www.fhi.no/en/) that is designed for highly variable viruses, and viruses that are likely to appear as co-infections between multiple strains, such as Hepatitis C Virus. The pipeline will identify the most likely major and minor strain in a sample sequenced with the Illumina platform. It will map the reads to these references using [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and create consensus sequences. For Hepatitis C Viruses the pipeline can also run a [GLUE-analysis](http://hcv-glue.cvr.gla.ac.uk/#/home) to identify drug resistance mutations.
maps Illumina reads to a reference genome and creates a consensus sequence.

## Requirements
The pipeline only requires [Nextflow](https://nextflow.io/) and [Docker](https://www.docker.com/) in order to run. Note that you must be able to run Docker as a non-root user as described [here](https://docs.docker.com/engine/install/linux-postinstall/#manage-docker-as-a-non-root-user).

> [!IMPORTANT]
> HCV-GLUE is currently only available with the Docker profile. We recommend that you always run the pipeline with Docker.

## Run the pipeline
The pipeline does not require any installation, only an internet connection. The pipeline is typically run with the following command:
```
nextflow run folkehelseinstituttet/hcv_illumina -r v1.0.6 \
    --input samplesheet.csv \
    --outdir <OUTDIR> \
    -profile docker
```
Nextflow will pull the pipeline from the GitHub repo automatically when it is launched. Here, v1.0.6 release is downloaded and run. You can omit `-r` and the code from the master branch will be used. But we always recommend that you specify either branch or release using `-r`.

If you want to download a local copy of the pipeline you can run:
```
nextflow pull folkehelseinstituttet/hcv_illumina -r v1.0.6
```
Again, `-r` is optional.


## Test the pipeline
To run a minimal test:
```
nextflow run folkehelseinstituttet/hcv_illumina -profile docker,test
```
This is only to see if you can get the pipeline up and running and will not run the entire pipeline such as HCV-GLUE. The results will be in a directory called `minimal_test`.

To run a full test on a real dataset type:
```
# First download the test dataset using nf-core/fetchngs
nextflow run nf-core/fetchngs -profile docker --input 'https://raw.githubusercontent.com/folkehelseinstituttet/hcv_illumina/refs/heads/dev/assets/test_ids.csv' --outdir full_test

# Then run the pipeline on the downloaded dataset
nextflow run folkehelseinstituttet/hcv_illumina -profile docker,test_full
```
This will download a HCV Illumina dataset from SRA and run the entire pipeline. The results will be in a directory called `full_test`.
Note that the pipeline will by default download and use the Kraken 2 [PlusPFP-8](https://benlangmead.github.io/aws-indexes/k2) database. This reqires at least 5 GB of free disk space and will take a few minutes to download and unpack. In addition, the default memory and cpu requirements of 12 cpus and 72 GB have been overridden to `50.GB` and `8`.

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
The pipeline can be run using different profiles, which will determine how the pipeline is executed. The default profile is `docker`, which uses Docker containers to run the pipeline. You can also use `singularity` or `conda` profiles if you prefer those environments. To set the profile use the `-profile` parameter, e.g.: `-profile docker/singularity/conda`.

### Provide parameters in a file
The different parameters can be provided in a file using the argument `-params-file path/to/params-file.yml`. The file can be either YAML-formatted:

```yml
input: 'samplesheet.csv'
outdir: 'results'
```

or JSON-formatted:
```json
{
  "input": "samplesheet.csv",
  "outdir": "results"
}
```

## Optional parameters
### Kraken2 databases
The pipeline uses [Kraken2](https://github.com/DerrickWood/kraken2) for two purposes. One is to classify the reads against a general database to get a broad overview of the taxonomic diversity within the sample (e.g., are there a lot of human reads?). The second is to classify the reads against a specific HCV-database and then use only the classified reads for the rest of the pipeline. This is done to reduce the computational load and time needed to run mapping and _de novo_ assembly.

By default, the pipeline will download and use the [PlusPFP-8 database](https://benlangmead.github.io/aws-indexes/k2) compiled by Ben Langmead for the broad classification. This requires the download and upacking of a fairly large file (>5 GB) and we recommend that you download and unpack this yourself and specify the path to the database using the `--kraken_all_db` parameter.

For the HCV-specific classification, the pipeline will use a very small and provided database which consists of around 200 different HCV strains. You can specify a custom HCV-datavase using the `--kraken_focused` paramter.

### HCV reference sequences
The database comes with a provided set of about 200 HCV reference sequences downloaded from NCBI. See the file [data/blast_db/HCVgenosubtypes_8.5.19_clean.fa](data/blast_db/HCVgenosubtypes_8.5.19_clean.fa). The fasta headers have been modified to begin with the genotype and subtype information (e.g., `1a`, `3b`, etc.) followed by an underscore and the NCBI accession number (e.g, `1a_AF009606`).  You can for example add or remove HCV strains by modifying this file. Remember to format the fasta headers accordingly. This file will then be used in the mapping and analysis of the de novo assembled contigs to identify genotype and subtype. You need to provide the path to this file like this: `--references /path/to/HCV-sequences.fasta`.

### Co-infections (major and minor strains)
The pipeline will first map all HCV-classified reads against all HCV reference sequences. Then it will identify the reference sequence with the most mapped reads and use the genotype and subtype information from this reference sequence to call major genotype and subtype. To identify a potential co-infection (minor strain), the pipeline will identify the reference that belongs to a different genotype than the major strain (expect for genotypes 1a and 1b which are considered different enough so that we can distinguish them in a co-infection) and has the highest coverage (i.e., percent of the genome covered by 5 or more reads). By default we have set a threshold of minimum 500 reads and 30% genome coverage in order to consider a strain as a minor strain at all. This can be overridden using the parameters `--minRead` and `--minCov`.

Note that there is a recombinant strain between subtypes 2k and 1b present in the database. If this is detected, the pipeline will not allow for a co-infection with either genotypes 1 or 2.

## Starting and stopping the pipeline
If the pipeline crashes, or stopped deliberately, it can be restarted from the last completed step by running the same command but with the `-resume` option. Read more about resuming a Nextflow pipeline [here](https://www.nextflow.io/docs/latest/cache-and-resume.html).

## Customizing the pipeline
Changing the arguments given to the various sub-tools can be done in several ways, perhaps the easiest is to create a custom config file. Described in more detail [here](https://nf-co.re/docs/usage/configuration#custom-configuration-files).

## Output files

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

