process HCVGLUE {

    tag "$meta.id"
    label 'process_low'

    // Single hermetic all-in-one image (MySQL 5.7 + GLUE engine). One task per BAM.
    // The HCV project data is NOT baked in — it is staged at runtime via params.hcvglue_db
    // and loaded by run-glue.sh (baked into the image at /usr/local/bin).
    // SECURITY: pin this image by tag+digest once published to GHCR (RESEARCH Open Q1 / T-05-02).
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/folkehelseinstituttet/hcvglue-allinone:1.1.114':
        'ghcr.io/folkehelseinstituttet/hcvglue-allinone:1.1.114' }"

    stageInMode 'copy' // Can't mount symlinked files into docker containers; dump is ~59 MB so copy is cheap

    input:
    tuple val(meta), path(bam), path(hcvglue_db)
    val hcvglue_threshold

    output:
    path("*.json")     , optional: true, emit: GLUE_json
    path("*.html")     , optional: true, emit: GLUE_html
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/bin/bash

    # Hermetic per-BAM HCV-GLUE: init task-local MySQL, load staged dump, run GLUE, exit.
    run-glue.sh ${hcvglue_db} ${hcvglue_threshold} ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GLUE project version: \$(ls *.json 2>/dev/null | head -n 1 | xargs -I {} grep -oP '"projectVersion"\\s*:\\s*"\\K[^"]+' {} 2>/dev/null || echo "unknown")
        GLUE engine version: \$(ls *.json 2>/dev/null | head -n 1 | xargs -I {} grep -oP '"engineVersion"\\s*:\\s*"\\K[^"]+' {} 2>/dev/null || echo "unknown")
        GLUE extension version: \$(ls *.json 2>/dev/null | head -n 1 | xargs -I {} grep -oP '"extensionVersion"\\s*:\\s*"\\K[^"]+' {} 2>/dev/null || echo "unknown")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/bin/bash

    # Create dummy outputs for stub run (named to match the live BAM-derived outputs)
    echo "Creating dummy outputs for stub run"

    cat > ${prefix}.major.nodup.json << 'EOF'
    {
        "projectVersion": "0.1.63",
        "engineVersion": "1.1.113",
        "extensionVersion": "0.1.33",
        "genotype": "3a",
        "subtype": "3a",
        "coverage": 95.5
    }
    EOF

    cat > ${prefix}.major.nodup.html << 'EOF'
    <!DOCTYPE html>
    <html>
    <head><title>HCV-GLUE Report</title></head>
    <body>
        <h1>HCV-GLUE Analysis Report</h1>
        <p>Genotype: 3a</p>
        <p>Coverage: 95.5%</p>
    </body>
    </html>
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GLUE project version: \$(ls *.json | head -n 1 | xargs -I {} grep -oP '"projectVersion"\\s*:\\s*"\\K[^"]+' {})
        GLUE engine version: \$(ls *.json | head -n 1 | xargs -I {} grep -oP '"engineVersion"\\s*:\\s*"\\K[^"]+' {})
        GLUE extension version: \$(ls *.json | head -n 1 | xargs -I {} grep -oP '"extensionVersion"\\s*:\\s*"\\K[^"]+' {})
    END_VERSIONS
    """
}
