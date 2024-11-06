# ![folkehelseinstituttet/viralseq](docs/images/logo-engelsk-hele-navnet.jpg#gh-light-mode-only) ![niph/viralseq](docs/images/logo-engelsk-hele-navnet-hvit.png#gh-dark-mode-only)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

## Introduction

**folkehelseinstituttet/viralseq** is a bioinformatics pipeline used at the [Norwegian Institute of Public Health](https://www.fhi.no/en/) that is designed for highly variable viruses, and viruses that are likely to appear as co-infections between multiple strains, such as Hepatitis C Virus. The pipeline will identify the most likely major and minor strain in a sample sequenced with the Illumina platform. It will map the reads to these references using [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and create consensus sequences. For Hepatitis C Viruses the pipeline can also run a [GLUE-analysis](http://hcv-glue.cvr.gla.ac.uk/#/home) to identify drug resistance mutations.
maps Illumina reads to a reference genome and creates a consensus sequence.

## Requirements
The pipeline only requires [Nextflow](https://nextflow.io/) and [Docker](https://www.docker.com/) in order to run. Note that you must be able to run Docker as a non-root user as described [here](https://docs.docker.com/engine/install/linux-postinstall/#manage-docker-as-a-non-root-user).

## Usage
To run a minimal test:
```
nextflow run folkehelseinstituttet/viralseq/main.nf -r master -profile docker,test_illumina
```


## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  niph/viralseq for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

