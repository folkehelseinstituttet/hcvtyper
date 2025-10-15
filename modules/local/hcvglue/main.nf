process HCVGLUE {

    label 'process_low'

    // Environment with Docker created using the podman package from conda-forge. Created using seqera containers.
    // Singularity image: https://wave.seqera.io/view/builds/bd-e170c468aba99710_1?_gl=1*j1kzwl*_gcl_au*MTM5MTA4NDk2NS4xNzUzNjg2MzUxLjY0MTQxNDc2Ni4xNzU2MzA2NTExLjE3NTYzMDY1MjU.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/podman:5.6.2--e170c468aba99710':
        'docker.io/ubuntu:22.04' }"

    stageInMode 'copy' // Can't mount symlinked files into docker containers

    input:
    path '*'
    val hcvglue_threshold

    output:
    path("*.json")     , optional: true, emit: GLUE_json
    path("*.html")     , optional: true, emit: GLUE_html
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/bin/bash

    # Run HCV-GLUE analysis using external script
    run_hcvglue.sh ${hcvglue_threshold}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GLUE project version: \$(ls *.json 2>/dev/null | head -n 1 | xargs -I {} grep -oP '"projectVersion"\\s*:\\s*"\\K[^"]+' {} 2>/dev/null || echo "unknown")
        GLUE engine version: \$(ls *.json 2>/dev/null | head -n 1 | xargs -I {} grep -oP '"engineVersion"\\s*:\\s*"\\K[^"]+' {} 2>/dev/null || echo "unknown")
        GLUE extension version: \$(ls *.json 2>/dev/null | head -n 1 | xargs -I {} grep -oP '"extensionVersion"\\s*:\\s*"\\K[^"]+' {} 2>/dev/null || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    #!/bin/bash

    # Create dummy outputs for stub run
    echo "Creating dummy outputs for stub run"

    cat > sample.json << 'EOF'
    {
        "projectVersion": "0.1.63",
        "engineVersion": "1.1.113",
        "extensionVersion": "0.1.33",
        "genotype": "3a",
        "subtype": "3a",
        "coverage": 95.5
    }
    EOF

    cat > sample.html << 'EOF'
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
