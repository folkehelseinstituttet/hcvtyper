process PREPARE_BOWTIE2_BUILD {
    tag "$meta.id"
    label 'process_single'

    conda "YOUR-TOOL-HERE"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'docker.io/pegi3s/biopython:1.78' }"

    input:
    tuple val(meta),  path(fasta)

    output:
    tuple val(meta), path("*_cov_*.fasta"), emit: contig, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gene_name = "${meta.gene}"
    """
    #!/usr/bin/env python3
    from Bio import SeqIO

    def extract_node_sequences(input_fasta):
        with open(input_fasta, "r") as infile:
            for record in SeqIO.parse(infile, "fasta"):
                if record.id.startswith("NODE"):
                    # Sanitize the sequence name to keep only the first two elements
                    # sanitized_contig_name = "_".join(record.id.split("_")[:2])
                    output_fasta = f"${prefix}.${gene_name}.{record.id}.fasta"
                    with open(output_fasta, "w") as outfile:
                        SeqIO.write(record, outfile, "fasta")

    if __name__ == "__main__":
        input_fasta = "$fasta"
        extract_node_sequences(input_fasta)
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """

}
