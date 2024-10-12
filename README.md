# ![niph/viralseq](docs/images/logo-engelsk-hele-navnet.jpg#gh-light-mode-only) ![niph/viralseq](docs/images/logo-engelsk-hele-navnet-hvit.png#gh-dark-mode-only)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

## TODO
- [ ] Rotavirus: Could include all contigs for each segment in the first alignment and phylogeny step (maybe add a length cutoff to remove very short contigs). And then use that tree to identify co-infection with multiple genotypes/strains.

## Introduction

**niph/viralseq** is a bioinformatics pipeline that is designed for highly variable viruses, and viruses that are likely to appear as co-infections between multiple strains, such as Hepatitis C Virus. The pipeline will identify the most likely major and minor strain in a sample sequenced with the Illumina platform. It will map the reads to these references using the [Tanoti mapper](https://github.com/vbsreenu/Tanoti), which is designed for highly variable datasets, and create consensus sequences. For Hepatitis C Viruses the pipeline can also run a [GLUE-analysis](http://hcv-glue.cvr.gla.ac.uk/#/home) to identify drug resistance mutations.
maps Illumina reads to a reference genome and creates a consensus sequence.

## Requirements
The pipeline only requires [Nextflow](https://nextflow.io/) and [Docker](https://www.docker.com/) in order to run. Note that you must be able to run Docker as a non-root user as described [here](https://docs.docker.com/engine/install/linux-postinstall/#manage-docker-as-a-non-root-user).

## Usage
First you should clone this repo:
```
git clone https://github.com/jonbra/viralseq.git
```

To run a minimal test:
```
nextflow run viralseq/main.nf -profile docker,test
```

Now, you can run the pipeline using:


```bash
nextflow run viralseq main.f \
   -profile docker \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

:::warning
Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).
:::




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

## Notes

Rotavirus IQTREE evolutionary models based on analyzing 50 samples:

   gene  best_fit_model count
   <chr> <chr>          <int>
 1 NSP1  GTR+F+I+G4        41
 2 NSP2  GTR+F+I+G4        38
 3 NSP2  GTR+F+G4           1
 4 NSP2  TIM3+F+G4          1
 5 NSP2  TIM3+F+I+G4        1
 6 NSP3  GTR+F+G4          31
 7 NSP3  GTR+F+I+G4         9
 8 NSP4  GTR+F+G4          34
 9 NSP5  TPM2u+F+G4        40
10 VP1   TIM3+F+I+G4       39
11 VP2   GTR+F+I+G4        37
12 VP3   TIM3+F+I+G4       40
13 VP4   TIM3+F+I+G4       32
14 VP4   GTR+F+I+G4         7
15 VP6   GTR+F+I+G4        42
16 VP7   GTR+F+I+G4        41
