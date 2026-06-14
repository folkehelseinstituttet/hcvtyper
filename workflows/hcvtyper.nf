/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { softwareVersionsToYAML                         } from '../subworkflows/nf-core/utils_nfcore_pipeline/main.nf'
include { GET_MAPPING_STATS as GET_MAPPING_STATS_WITHDUP } from '../subworkflows/local/get_mapping_stats'
include { GET_MAPPING_STATS as GET_MAPPING_STATS_MARKDUP } from '../subworkflows/local/get_mapping_stats'
include { TARGETED_MAPPING                               } from '../subworkflows/local/targeted_mapping'
include { CONTAMINATION_CHECK                            } from '../subworkflows/local/contamination_check/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
//include { fromSamplesheet                    } from 'plugin/nf-validation'
include { INPUT_CHECK                        } from '../subworkflows/local/input_check'
include { BOWTIE2_BUILD                      } from '../modules/nf-core/bowtie2/build/main'
include { BLAST_MAKEBLASTDB                  } from '../modules/nf-core/blast/makeblastdb/main'
include { FASTP                              } from '../modules/nf-core/fastp/main'
include { FASTQC as FASTQC_RAW               } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIM              } from '../modules/nf-core/fastqc/main'
include { CUTADAPT                           } from '../modules/nf-core/cutadapt/main'
include { PRINSEQPLUSPLUS                    } from '../modules/nf-core/prinseqplusplus/main'
include { MULTIQC                            } from '../modules/nf-core/multiqc/main'
include { KRAKEN2_KRAKEN2                    } from '../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_FOCUSED } from '../modules/nf-core/kraken2/kraken2/main'
include { SPADES                             } from '../modules/nf-core/spades/main'
include { BLAST_BLASTN                       } from '../modules/nf-core/blast/blastn/main'
include { BOWTIE2_ALIGN                      } from '../modules/nf-core/bowtie2/align/main'
include { SAMTOOLS_SORMADUP                  } from '../modules/nf-core/samtools/sormadup/main'
include { PICARD_COLLECTINSERTSIZEMETRICS    } from '../modules/nf-core/picard/collectinsertsizemetrics/main'
include { UNTAR as UNTAR_KRAKEN_ALL          } from '../modules/nf-core/untar/main'
include { UNTAR as UNTAR_KRAKEN_FOCUSED      } from '../modules/nf-core/untar/main'

//
// Local modules
//
include { INSTRUMENTID                       } from '../modules/local/instrumentid/main'
include { BLASTPARSE                         } from '../modules/local/blastparse/main'
include { PARSEFIRSTMAPPING                  } from '../modules/local/parsefirstmapping/main'
include { GLUEPARSE as HCV_GLUE_PARSER       } from '../modules/local/glueparse/main'
include { HCVGLUE                            } from '../modules/local/hcvglue/main'
include { SUMMARIZE                          } from '../modules/local/summarize/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HCVTYPER {

    main:

    // GUARD (CR-01): The de novo-informed candidate model ranks N candidates, but the
    // on-disk compatibility shim only writes two filename slots (`.major.`/`.minor.`,
    // backed by `_major.fa`/`_minor.fa`). With n_candidates > 2, a passing 3rd+ candidate
    // would be silently mapped against the 2nd candidate's reference — wrong-reference data,
    // no crash. The N-slot filename migration (`.cand1.`/`.cand2.`/...) is deferred to
    // Phase 9 (COMPAT-02). Until then, fail loudly rather than emit wrong results.
    if (params.n_candidates > 2) {
        error "params.n_candidates = ${params.n_candidates} is not yet supported: the candidate-to-filename shim only has two slots (major/minor). Set --n_candidates to 1 or 2. N>2 support arrives with the Phase 9 filename-slot migration (COMPAT-02)."
    }

    ch_versions = Channel.empty()

    //
    // CONFIG FILES
    //
    ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )
    reads = INPUT_CHECK.out.reads
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // Prepare Kraken2 database for all domains of life
    //
    ch_kraken_all_db = Channel.empty()
    if (params.kraken_all_db) {
        if (params.kraken_all_db.endsWith('.tar.gz')) {
            UNTAR_KRAKEN_ALL (
                [ [:], params.kraken_all_db ] // Add empty meta map
            )
            ch_kraken_all_db = UNTAR_KRAKEN_ALL.out.untar.map { it[1] } // Do not extract the meta map which is emitted by default
            ch_versions = ch_versions.mix(UNTAR_KRAKEN_ALL.out.versions.first())
        } else {
            ch_kraken_all_db = Channel.value(file(params.kraken_all_db))
        }
    }

    //
    // Prepare Kraken2 database for HCV only
    //
    ch_kraken_focused_db = Channel.empty()
        if (params.kraken_focused_db.endsWith('.tar.gz')) {
            UNTAR_KRAKEN_FOCUSED (
                [ [:], params.kraken_focused_db ] // Add empty meta map
            )
            ch_kraken_focused_db = UNTAR_KRAKEN_FOCUSED.out.untar.map { it[1] } // Do not extract the meta map which is emitted by default
            ch_versions = ch_versions.mix(UNTAR_KRAKEN_FOCUSED.out.versions.first())
        } else {
            ch_kraken_focused_db = Channel.value(file(params.kraken_focused_db))
        }

    //
    // MODULE: Identify the instrument ID
    //
    INSTRUMENTID (
        reads
    )
    ch_versions = ch_versions.mix(INSTRUMENTID.out.versions.first())

    //
    // MODULE: Run Bowtie2_build to create a reference index
    //

    BOWTIE2_BUILD(
        [ [], file(params.references) ] // Add empty meta map before the reference file path
    )
    ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)

    //
    // MODULE: Create Blast database from reference sequences
    //
    BLAST_MAKEBLASTDB (
        [ [id:"blastdb/"], file(params.references) ] // Add meta map before the reference file path
    )
    ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC_RAW (
        reads
    )
    ch_versions = ch_versions.mix(FASTQC_RAW.out.versions)

    //
    // MODULE: Trim reads with either Cutadapt or Fastp
    //
    if (params.trimmer == "fastp") {
        FASTP(
            reads
                .map { meta, fastq ->
                    def n = fastq[0].countFastq() // Count fastq reads in the R1 fastq file
                    return [meta, fastq, n] // Add the count as the last element in the tuple
                }
                .filter { _meta, _fastq, n -> n > 0 } // Filter out empty fastq files because Fastp does not produce correct gzipped file if input fastq is empty
                .map { meta, fastq, _n -> [ meta, fastq, [] ] }, // Remove the count and add empty element to get the structure correct for Fastp
            false,
            false,
            false
        )
        ch_trimmed_reads = FASTP.out.reads
        ch_versions = ch_versions.mix(FASTP.out.versions)
    } else if (params.trimmer == "cutadapt") {
        CUTADAPT(
            reads
        )
        ch_trimmed_reads = CUTADAPT.out.reads
        ch_versions = ch_versions.mix(CUTADAPT.out.versions)
    }

    //
    // MODULE: Run FastQC on trimmed reads
    //
    FASTQC_TRIM (
        ch_trimmed_reads
    )
    ch_versions = ch_versions.mix(FASTQC_TRIM.out.versions.first())

    //
    // MODULE: Remove low complexity reads with Prinseq++ (optional)
    //

    // Initialise a default empty channel
    ch_prinseq_log = Channel.empty()
    if (params.filter_low_complexity) {
        // NOTE:
        // In some cases there are empty fastq files after trimming. Remove these before PRINSEQPLUSPLUS
        ch_prinseq = ch_trimmed_reads
            .map { meta, fastq ->
                def n = fastq[0].countFastq() // Count fastq reads in the R1 fastq file
                return [meta, fastq, n] // Add the count as the last element in the tuple
            }
            .filter { _meta, _fastq, n -> n > 0 } // Filter out empty fastq files
            .map { meta, fastq, _n -> [meta, fastq] } // Remove the count to get the channel structure correct for PRINSEQPLUSPLUS

        PRINSEQPLUSPLUS (
            ch_prinseq
        )
        ch_versions = ch_versions.mix(PRINSEQPLUSPLUS.out.versions.first())

        // Kraken2 input is from PRINSEQPLUSPLUS
        ch_kraken_input = PRINSEQPLUSPLUS.out.good_reads

        // For MultiQC later
        ch_prinseq_log = PRINSEQPLUSPLUS.out.log.collect { it[1] }.ifEmpty([])
    } else {
        // Skip PRINSEQPLUSPLUS and use trimmed reads as input to Kraken2
        ch_kraken_input = ch_trimmed_reads
    }

    //
    // MODULE: Run Kraken2 to classify reads
    //
    KRAKEN2_KRAKEN2 (
        ch_kraken_input,
        ch_kraken_all_db,
        false,
        false
    )
    ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions.first())

    //
    // MODULE: Run Kraken2 to identify target viral reads
    //
    KRAKEN2_FOCUSED (
        ch_kraken_input,
        ch_kraken_focused_db,
        params.save_output_fastqs,
        params.save_reads_assignment
    )

    //
    // MODULE: Run Spades to assemble classified reads
    //

    // Create input read channel for SPADES.
    // A tuple with meta, paired Illumina reads, and empty elements for pacbio and nanopore reads
    // Filter in case of no classified reads emitted by KRAKEN2_FOCUSED
    ch_reads = KRAKEN2_FOCUSED.out.classified_reads_fastq
        .map { meta, fastq ->
            def n = fastq[0].countFastq() // Count fastq reads in the R1 fastq file
            return [meta, fastq, n] // Add the count as the last element in the tuple
        }
        .filter { _meta, _fastq, n -> n > 0 } // Filter out empty fastq files
        .map { meta, fastq, _n -> [ meta, fastq, [], [] ] } // Recreate the channel structure correct for SPADES

    if (!params.skip_assembly) {
            SPADES (
                ch_reads,
                [], // Empty input channel. Can be used to specify hmm profile
                []  // Empty input channel. Placeholder for separate specification of reads.
            )
            ch_versions = ch_versions.mix(SPADES.out.versions.first())

            //
            // MODULE: Blast assembled contigs against viral references.
            //
            // NOTE:
            // In some cases there is an empty contig file produced by Spades. Filter out these
            ch_blastn = SPADES.out.contigs
                .map { meta, contigs ->
                def n = contigs.countFasta() // Count fasta records
                return [meta, contigs, n] // Add the count as the last element in the tuple
            }
            .filter { _meta, _contigs, n -> n > 0 } // Filter out empty fasta files
            .map { meta, contigs, _n -> [meta, contigs] } // Return the count to get the channel structure correct for BLASTN_BLASTN
            BLAST_BLASTN (
                ch_blastn,
                BLAST_MAKEBLASTDB.out.db,
                [], // taxidlist - empty, no taxonomic filtering
                "", // taxids - empty string, no taxonomic filtering
                false // negative_tax - false, not using negative filtering
            )
            ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())

            //
            // MODULE: Parse blast output
            //
            ch_blastparse = BLAST_BLASTN.out.txt.join(SPADES.out.contigs) // Create input channel that holds val(meta), path(blast_out), path(contigs)
            BLASTPARSE (
                ch_blastparse,
                file(params.references),
                params.agens
            )
            ch_versions = ch_versions.mix(BLASTPARSE.out.versions.first())

            //
            // SUBWORKFLOW: Detect cross-sample contamination via all-vs-all BLAST
            //
            if (!params.skip_contamination_check) {
                CONTAMINATION_CHECK(
                    SPADES.out.contigs,
                    Channel.empty(),  // fastp JSONs — not wired in main pipeline
                    Channel.empty()   // GLUE JSONs  — not wired in main pipeline
                )
                ch_versions = ch_versions.mix(CONTAMINATION_CHECK.out.versions)
            }
    }

    //
    // MODULE: Map classified reads against all references
    //
    BOWTIE2_ALIGN (
        KRAKEN2_FOCUSED.out.classified_reads_fastq,
        BOWTIE2_BUILD.out.index,
        [ [], file(params.references) ], // Add empty meta map and reference fasta for CRAM support
        false, // Do not save unmapped reads
        true // Sort bam file
    )
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())
    ch_aligned = BOWTIE2_ALIGN.out.bam

    //
    // SUBWORKFLOW: Get mapping statistics with duplicates included
    //
    GET_MAPPING_STATS_WITHDUP (
        ch_aligned, // Channel: [ val(meta), [ bam ] ]
        Channel.value(file(params.references)).map { [ [:], it ] } // Add empty meta map before the reference file path
    )
    versions = GET_MAPPING_STATS_WITHDUP.out.versions

    // Remove duplicate reads
    SAMTOOLS_SORMADUP (
        ch_aligned,
        [ [], file(params.references) ]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORMADUP.out.versions.first())

    //
    // MODULE: Collect insert size metrics with duplicates removed
    //
    PICARD_COLLECTINSERTSIZEMETRICS (
        SAMTOOLS_SORMADUP.out.bam
    )
    ch_versions = ch_versions.mix(PICARD_COLLECTINSERTSIZEMETRICS.out.versions.first())

    //
    // SUBWORKFLOW: Get mapping statistics with duplicates removed
    //
    GET_MAPPING_STATS_MARKDUP (
        SAMTOOLS_SORMADUP.out.bam,
        Channel.value(file(params.references)).map { [ [:], it ] } // Add empty meta map before the reference file path
    )
    versions = GET_MAPPING_STATS_MARKDUP.out.versions

    //
    // MODULE: Identify the two references with most mapped reads, duplicates included
    //
    ch_parsefirstmapping = GET_MAPPING_STATS_WITHDUP.out.idxstats.join(GET_MAPPING_STATS_WITHDUP.out.tsv) // val(meta), path(idxstats), path(tsv)
        .filter { _meta, _idxstats, tsv ->
            tsv.size() > 0 // Filter out empty tsv files
        }

    PARSEFIRSTMAPPING (
        // Join idxstats and depth on the meta map
        ch_parsefirstmapping,
        file(params.references),
        file("${projectDir}/bin/genotype_utils.R", checkIfExists: true)
    )

    //
    // SUBWORKFLOW: Map reads against EACH neutrally-ranked candidate reference (D-04 / REFSEL-03)
    //
    // The two asymmetric major/minor alias routes are collapsed into ONE
    // uniform per-candidate fan-out over the long-format candidates CSV. We splitCsv ALL rows
    // (NOT elements[0] — the legacy code only read row[0] because the wide CSV was single-row;
    // the long-format candidates CSV is one row PER candidate, Assumption A2) and emit one
    // channel element per candidate. Each candidate carries its own per-rank meta + FASTA.
    //
    // The per-rank FASTA paths (_major.fa / _minor.fa) come from the legacy major_mapping /
    // minor_mapping emits (D-06 shim) so the TARGETED_MAPPING `meta.reference` enrichment
    // (fasta basename split, e.g. <ref>_major) stays byte-identical to the legacy filenames.
    //
    // Join all per-sample inputs by meta.id: the candidates CSV, the per-rank FASTAs, and the
    // classified reads. The legacy major_mapping/minor_mapping emits are `optional: true`
    // (a single-candidate sample emits no `_minor.fa`, a no-candidate sample emits neither),
    // so both legacy joins use `remainder: true` — otherwise a missing optional emit would
    // silently DROP the whole sample (including its passing major candidate). A null FASTA is
    // tolerated and the candidate guarded out below.
    ch_candidate_mapping = PARSEFIRSTMAPPING.out.candidates
        .join(PARSEFIRSTMAPPING.out.major_mapping, remainder: true)        // meta, candidates_csv, wide_csv?, major_fasta?
        .join(PARSEFIRSTMAPPING.out.minor_mapping, remainder: true)        // ..., wide_csv?, minor_fasta?
        .join(KRAKEN2_FOCUSED.out.classified_reads_fastq)                  // ..., classified_reads
        .flatMap { meta, candidates_csv, _wide1, major_fasta, _wide2, minor_fasta, classified_reads ->
            // Iterate ALL candidate rows (one element per candidate), not just row[0].
            def rows = candidates_csv.splitCsv( header: true, sep:',' )
            rows.collect { row ->
                // Lift the candidate row into the meta map. Carry candidate_rank,
                // candidate_ref and confirmation_status (R-emitted STRINGS — never coerced
                // to Integer here, so a single/no-candidate NA field can never crash, Pitfall 3).
                def new_meta = meta + row

                // Fail loudly if the candidate row's sample disagrees with the pipeline meta.id.
                assert new_meta.id == new_meta.sample : "Metadata mismatch: id=${new_meta.id}, sample=${new_meta.sample}"

                // Per-candidate meta uniqueness: same-sample candidates share meta.id, so an
                // id-only join inside TARGETED_MAPPING would cross-pair index/fasta/reads (Pitfall 2).
                // We KEEP meta.id == sample (so output filenames stay <sample>.<ref>... for summarize.R),
                // and rely on TARGETED_MAPPING joining by the FULL meta map: candidate_rank + candidate_ref
                // (and the `reference` enrichment inside the subworkflow) differ per candidate, so the
                // whole-meta join key is distinct for every candidate of a sample.
                def rank = new_meta.candidate_rank.toString()

                // Pick the per-rank FASTA from the legacy shim emits (rank 1 -> _major.fa, else _minor.fa).
                def fasta = (rank == '1') ? major_fasta : minor_fasta

                tuple(new_meta, fasta, classified_reads)
            }
        }
        // Route on the R-emitted per-candidate STRING (confirmation_status), never a Groovy
        // numeric coercion of possibly-NA fields (Pitfall 3 — avoids NA.toInteger() crash).
        // confirmation_status == 'pass' generalizes the legacy gate: a candidate is mapped only
        // when its own reads>minRead && cov>minCov comparison passed (rank 1 == the major gate,
        // rank 2 == a passing second candidate). At default N=2 on a single-strain fixture the
        // second candidate is 'below_threshold', so the two-slot topology is preserved (D-06).
        // A passing candidate always has its per-rank FASTA written by the selection script, so
        // the null-FASTA guard only drops below-threshold/absent candidates (defensive).
        .filter { entry -> entry[0]['confirmation_status'] == 'pass' && entry[1] != null }

    TARGETED_MAPPING(
        ch_candidate_mapping, // val(meta), path(fasta), path(reads)
    )
    ch_versions = ch_versions.mix(TARGETED_MAPPING.out.versions)

    //
    // MODULE: Run GLUE genotyping and resistance annotation for HCV
    //
    if (!params.skip_hcvglue) {
        HCVGLUE (
            TARGETED_MAPPING.out.aligned.collect({it[1]}).collect(), // Collect all candidate BAMs (T-2 lockstep). Can only have one GLUE process running
            params.hcvglue_threshold
        )
        ch_versions = ch_versions.mix(HCVGLUE.out.versions)

        // Collect all glue reports and parse them
        HCV_GLUE_PARSER (
            HCVGLUE.out.GLUE_json.collect()
        )
        ch_versions = ch_versions.mix(HCV_GLUE_PARSER.out.versions)
    }

    //
    // MODULE: Summarize
    //
    // Create channel with this structure: path(stats), path(depth), path(blast), path(json)
    // Collect all the files in separate channels for clarity. Don't need the meta
    ch_sequence_id      = INSTRUMENTID.out.id.collect({it[1]})
    if (params.trimmer == "fastp") {
        ch_trimmed_reads         = FASTP.out.log.collect({it[1]})
    } else if (params.trimmer == "cutadapt") {
        ch_trimmed_reads         = CUTADAPT.out.log.collect({it[1]})
    }
    ch_classified_reads = KRAKEN2_FOCUSED.out.report.collect({it[1]})
    // Phase 7 (ASUP-02): stage the Phase-6 long-format *.candidates.csv alongside
    // the legacy *.parsefirstmapping.csv into parsefirst_mapping/ so summarize.R
    // can read it for the genotype-level assembly-support join.
    ch_summarize_first_mapping = PARSEFIRSTMAPPING.out.csv.collect({it[1]}).mix(PARSEFIRSTMAPPING.out.candidates.collect({it[1]})).collect()
    // T-2 lockstep: the single per-candidate fan-out already contains ALL candidate
    // stats/depth/consensus, so each former .mix(MAJOR..., MINOR...) pair collapses to the
    // single TARGETED_MAPPING.out.* . Missing any one would silently halve the stats.
    ch_stats_withdup    = TARGETED_MAPPING.out.stats_withdup.collect({it[1]})
    ch_stats_markdup    = TARGETED_MAPPING.out.stats_markdup.collect({it[1]})
    ch_depth            = TARGETED_MAPPING.out.depth.collect({it[1]})
    // De novo / BLAST evidence (PLUMB-01/PLUMB-02): collect the parsed BLASTPARSE
    // CSVs (*.blastparse.csv) and the per-contig table (*_blast_out.csv) into one
    // staged channel. BLASTPARSE is invoked only inside if (!params.skip_assembly),
    // so its .out attribute is undefined on a skip-assembly run -- referencing it
    // unconditionally is a hard Nextflow error (process not invoked), which .ifEmpty
    // cannot rescue. Guard the channel construction with the same condition (mirroring
    // the ch_glue if/else below): skip-assembly yields [] -> empty denovo/ staging dir
    // -> NA de novo columns + no dropped rows (the PLUMB-02 path).
    if (!params.skip_assembly) {
        // Phase 7 (ASUP-02): also stage the per-subtype *.assembly_support.csv into
        // denovo/ so summarize.R can join it to candidates at genotype level.
        ch_denovo = BLASTPARSE.out.csv.collect({it[1]}).mix(BLASTPARSE.out.blast_res.collect({it[1]})).mix(BLASTPARSE.out.support.collect({it[1]})).collect().ifEmpty([])
    } else {
        ch_denovo = []
    }
    if (params.agens == "HCV" && !params.skip_hcvglue) {
        ch_glue = HCV_GLUE_PARSER.out.GLUE_summary
    } else {
        ch_glue = []
    }
    ch_variation = TARGETED_MAPPING.out.variation.collect()
    ch_consensus_distance = TARGETED_MAPPING.out.consensus_distance.collect({it[1]})

    SUMMARIZE (
        workflow.manifest.version,
        workflow.manifest.name,
        file(params.input),
        ch_trimmed_reads.collect(),
        ch_classified_reads.collect(),
        ch_summarize_first_mapping,
        ch_stats_withdup.collect(),
        ch_stats_markdup.collect(),
        ch_depth.collect(),
        ch_denovo,
        ch_glue,
        ch_sequence_id.collect(),
        ch_variation.collect(),
        ch_consensus_distance.collect(),
        file("${projectDir}/bin/genotype_utils.R"),
        file("${projectDir}/bin/denovo_confirm.R"),
        file("${projectDir}/bin/denovo_layer.R"),
        file("${projectDir}/bin/assembly_support_join.R"),
        file("${projectDir}/bin/classify_roles.R"),
    )
    ch_versions = ch_versions.mix(SUMMARIZE.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: "software_mqc_versions.yml",
            sort: true,
            newLine: true
        )
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    def summary_params = paramsSummaryMap(workflow)
    workflow_summary    = WorkflowHCVTyper.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowHCVTyper.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{it[1]}.ifEmpty([]))
    if (params.trimmer == "cutadapt") {
        ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT.out.log.collect{it[1]}.ifEmpty([]))
    } else if (params.trimmer == "fastp") {
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]}.ifEmpty([]))
    }
    ch_multiqc_files = ch_multiqc_files.mix(ch_prinseq_log)
    ch_multiqc_files = ch_multiqc_files.mix(ch_trimmed_reads.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIM.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_KRAKEN2.out.report.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_FOCUSED.out.report.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(SUMMARIZE.out.mqc.collect())
    if (!params.skip_assembly && !params.skip_contamination_check) {
        ch_multiqc_files = ch_multiqc_files.mix(CONTAMINATION_CHECK.out.mqc)
    }

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
