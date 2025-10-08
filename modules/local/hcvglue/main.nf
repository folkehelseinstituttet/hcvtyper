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

    # Detect execution environment and setup container runtime
    if [ "\${NXF_EXECUTOR}" = "local" ] && [ -z "\${SINGULARITY_CONTAINER}" ] && [ -z "\${DOCKER_CONTAINER_ID}" ]; then
        echo "Running in conda/local environment"
        # Check if Docker is available on host
        if command -v docker &> /dev/null && docker info &> /dev/null 2>&1; then
            CONTAINER_CMD="docker"
            echo "Using host Docker daemon"
        elif command -v podman &> /dev/null; then
            CONTAINER_CMD="podman"
            echo "Using host Podman"
        else
            echo "ERROR: No container runtime available. Please install Docker or Podman."
            echo "For conda: conda install -c conda-forge docker-compose"
            exit 1
        fi
    elif [ -n "\${SINGULARITY_CONTAINER}" ]; then
        echo "Running in Singularity container"
        # Install Docker CLI in Singularity (if not already available)
        if command -v apt-get &> /dev/null; then
            apt-get update && apt-get install -y docker.io
        elif command -v apk &> /dev/null; then
            apk add --no-cache docker-cli
        elif command -v yum &> /dev/null; then
            yum install -y docker
        fi
        CONTAINER_CMD="docker"
        echo "Using Docker CLI in Singularity"
    else
        echo "Running in Docker container"
        # Install Docker CLI in container
        if command -v apt-get &> /dev/null; then
            apt-get update && apt-get install -y docker.io
        elif command -v apk &> /dev/null; then
            apk add --no-cache docker-cli
        fi
        CONTAINER_CMD="docker"
        echo "Using Docker CLI in container"
    fi

    # Verify container runtime works
    if ! \$CONTAINER_CMD info &> /dev/null; then
        echo "ERROR: Container runtime '\$CONTAINER_CMD' not accessible"
        echo "Make sure Docker daemon is running and accessible"
        exit 1
    fi

    echo "Container runtime: \$CONTAINER_CMD"

    # Remove the container in case it is already running
    if \$CONTAINER_CMD ps -a --filter "name=gluetools-mysql" --format '{{.Names}}' | grep -q "^gluetools-mysql\$"; then
        echo "Container 'gluetools-mysql' is running or exists."
        # Stop the container
        \$CONTAINER_CMD stop gluetools-mysql || echo "Failed to stop container or it is already stopped."

        # Wait for removal if already in progress
        while \$CONTAINER_CMD ps -a --filter "name=gluetools-mysql" --format '{{.State}}' | grep -q "removing"; do
            echo "Container 'gluetools-mysql' is being removed. Waiting..."
            sleep 1
        done

        # Remove the container
        \$CONTAINER_CMD rm gluetools-mysql || echo "Failed to remove container or it has already been removed."
        echo "Container 'gluetools-mysql' has been stopped and removed."
    else
        echo "Container 'gluetools-mysql' is not running or does not exist."
    fi

    # Pull the latest images
    \$CONTAINER_CMD pull cvrbioinformatics/gluetools-mysql:latest
    \$CONTAINER_CMD pull cvrbioinformatics/gluetools:latest

    # Start the gluetools-mysql container
    \$CONTAINER_CMD run --detach --name gluetools-mysql cvrbioinformatics/gluetools-mysql:latest

    # Install the pre-built GLUE HCV project
    TIMEOUT=300
    START_TIME=\$(date +%s)

    until \$CONTAINER_CMD exec gluetools-mysql mysql --user=root --password=root123 -e "status" &> /dev/null
    do
        echo "Waiting for database connection..."
        # Wait for two seconds before checking again
        sleep 2

        # Check if the timeout has been reached
        CURRENT_TIME=\$(date +%s)
        ELAPSED_TIME=\$((CURRENT_TIME - START_TIME))
        if [ \$ELAPSED_TIME -ge \$TIMEOUT ]; then
            echo "Timeout reached. Exiting script."
            exit 1
        fi
    done
    \$CONTAINER_CMD exec gluetools-mysql installGlueProject.sh ncbi_hcv_glue

    # Make a for loop over all bam files and run HCV-GLUE
    ## Adding || true to the end of the command to prevent the pipeline from failing if the bam file is not valid
    for bam in \$(ls *.bam 2>/dev/null || echo ""); do
        if [ -z "\$bam" ]; then
            echo "No BAM files found"
            break
        fi

        echo "Processing \$bam for JSON report..."
        # First create json report
        \$CONTAINER_CMD run --rm \\
           --name gluetools \\
            -v \$PWD:/opt/bams \\
            -w /opt/bams \\
            --link gluetools-mysql \\
            cvrbioinformatics/gluetools:latest gluetools.sh \\
             -p cmd-result-format:json \\
            -EC \\
            -i project hcv module phdrReportingController invoke-function reportBam \${bam} ${hcvglue_threshold} > \${bam%".bam"}.json || true
    done

    # Then create html report
    for bam in \$(ls *.bam 2>/dev/null || echo ""); do
        if [ -z "\$bam" ]; then
            echo "No BAM files found for HTML reports"
            break
        fi

        echo "Processing \$bam for HTML report..."
        \$CONTAINER_CMD run --rm \\
            --name gluetools \\
            -v \$PWD:/opt/bams \\
            -w /opt/bams \\
            --link gluetools-mysql \\
            cvrbioinformatics/gluetools:latest gluetools.sh \\
            --console-option log-level:FINEST \\
            --inline-cmd project hcv module phdrReportingController invoke-function reportBamAsHtml \${bam} ${hcvglue_threshold} \${bam%".bam"}.html || true
    done

    # Cleanup: stop the MySQL container
    echo "Cleaning up containers..."
    \$CONTAINER_CMD stop gluetools-mysql || true
    \$CONTAINER_CMD rm gluetools-mysql || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        container_runtime: \$(\$CONTAINER_CMD --version 2>/dev/null || echo "unknown")
        execution_env: "\${NXF_EXECUTOR:-local}"
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
