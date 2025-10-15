#!/bin/bash

# run_hcvglue.sh - Script to run HCV-GLUE analysis using Docker containers
# Arguments: $1 = hcvglue_threshold

set -euo pipefail

HCVGLUE_THRESHOLD=$1

# Detect execution environment and setup container runtime
echo "DEBUG: NXF_EXECUTOR=${NXF_EXECUTOR:-not_set}"
echo "DEBUG: SINGULARITY_CONTAINER=${SINGULARITY_CONTAINER:-not_set}"
echo "DEBUG: DOCKER_CONTAINER_ID=${DOCKER_CONTAINER_ID:-not_set}"

if [ "${NXF_EXECUTOR:-}" = "local" ] && [ -z "${SINGULARITY_CONTAINER:-}" ] && [ -z "${DOCKER_CONTAINER_ID:-}" ]; then
    echo "Running in conda/local environment"
    # Check if Docker is available on host
    if command -v docker &> /dev/null && docker info &> /dev/null 2>&1; then
        CONTAINER_CMD="docker"
        echo "Using host Docker daemon"
    elif command -v podman &> /dev/null && podman info &> /dev/null 2>&1; then
        CONTAINER_CMD="podman"
        echo "Using host Podman daemon"
    else
        echo "ERROR: Neither Docker nor Podman found on host"
        exit 1
    fi
elif [ -n "${SINGULARITY_CONTAINER:-}" ]; then
    echo "Running in Singularity container"
    # Install Docker CLI inside Singularity container
    if command -v apt-get &> /dev/null; then
        apt-get update && apt-get install -y docker.io
    elif command -v apk &> /dev/null; then
        apk add --no-cache docker-cli
    fi
    CONTAINER_CMD="docker"
    echo "Using Docker CLI in Singularity"
elif [ -n "${DOCKER_CONTAINER_ID:-}" ]; then
    echo "Running inside Docker container"
    # Install Docker CLI inside Docker container
    if command -v apt-get &> /dev/null; then
        apt-get update && apt-get install -y docker.io
    elif command -v apk &> /dev/null; then
        apk add --no-cache docker-cli
    fi
    CONTAINER_CMD="docker"
    echo "Using Docker CLI in container"
else
    echo "Execution environment not clearly detected, checking for available container runtimes..."
    # Check if we're in any kind of container and install Docker CLI
    if [ -f /.dockerenv ] || [ -n "${container:-}" ] || [ -n "${SINGULARITY_NAME:-}" ]; then
        echo "Detected container environment, installing Docker CLI"
        if command -v apt-get &> /dev/null; then
            apt-get update && apt-get install -y docker.io
        elif command -v apk &> /dev/null; then
            apk add --no-cache docker-cli
        elif command -v yum &> /dev/null; then
            yum install -y docker
        fi
        CONTAINER_CMD="docker"
        echo "Using Docker CLI in container"
    else
        # Fallback: try to detect container runtime on host
        if command -v docker &> /dev/null && docker info &> /dev/null; then
            CONTAINER_CMD="docker"
            echo "Using Docker on host"
        elif command -v podman &> /dev/null && podman info &> /dev/null; then
            CONTAINER_CMD="podman"
            echo "Using Podman on host"
        else
            echo "ERROR: Neither Docker nor Podman found or accessible"
            echo "Make sure a container runtime is installed and running"
            exit 1
        fi
    fi
fi

# Verify container runtime is accessible
if ! $CONTAINER_CMD info &> /dev/null; then
    echo "ERROR: Container runtime '$CONTAINER_CMD' not accessible"
    echo "Make sure Docker daemon is running and accessible"
    exit 1
fi

echo "Container runtime: $CONTAINER_CMD"

# Remove the container in case it is already running
if $CONTAINER_CMD ps -a --filter "name=gluetools-mysql" --format '{{.Names}}' | grep -q "^gluetools-mysql$"; then
    echo "Container 'gluetools-mysql' is running or exists."
    # Stop the container
    $CONTAINER_CMD stop gluetools-mysql || echo "Failed to stop container or it is already stopped."

    # Wait for removal if already in progress
    while $CONTAINER_CMD ps -a --filter "name=gluetools-mysql" --format '{{.State}}' | grep -q "removing"; do
        echo "Container 'gluetools-mysql' is being removed. Waiting..."
        sleep 1
    done

    # Remove the container
    $CONTAINER_CMD rm gluetools-mysql || echo "Failed to remove container or it has already been removed."
    echo "Container 'gluetools-mysql' has been stopped and removed."
else
    echo "Container 'gluetools-mysql' is not running or does not exist."
fi

# Pull the latest images
$CONTAINER_CMD pull cvrbioinformatics/gluetools-mysql:latest
$CONTAINER_CMD pull cvrbioinformatics/gluetools:latest

# Start the gluetools-mysql container
$CONTAINER_CMD run --detach --name gluetools-mysql cvrbioinformatics/gluetools-mysql:latest

# Install the pre-built GLUE HCV project
TIMEOUT=300
START_TIME=$(date +%s)

until $CONTAINER_CMD exec gluetools-mysql mysql --user=root --password=root123 -e "status" &> /dev/null
do
    echo "Waiting for database connection..."
    # Wait for two seconds before checking again
    sleep 2

    # Check if the timeout has been reached
    CURRENT_TIME=$(date +%s)
    ELAPSED_TIME=$((CURRENT_TIME - START_TIME))
    if [ $ELAPSED_TIME -ge $TIMEOUT ]; then
        echo "Timeout reached. Exiting script."
        exit 1
    fi
done
$CONTAINER_CMD exec gluetools-mysql installGlueProject.sh ncbi_hcv_glue

# Make a for loop over all bam files and run HCV-GLUE
## Adding || true to the end of the command to prevent the pipeline from failing if the bam file is not valid
for bam in $(ls *.bam 2>/dev/null || echo ""); do
    if [ -z "$bam" ]; then
        echo "No BAM files found"
        break
    fi

    echo "Processing $bam for JSON report..."
    # First create json report
    $CONTAINER_CMD run --rm \
       --name gluetools \
        -v $PWD:/opt/bams \
        -w /opt/bams \
        --link gluetools-mysql \
        cvrbioinformatics/gluetools:latest gluetools.sh \
         -p cmd-result-format:json \
        -EC \
        -i project hcv module phdrReportingController invoke-function reportBam ${bam} $HCVGLUE_THRESHOLD > ${bam%".bam"}.json || true
done

# Then create html report
for bam in $(ls *.bam 2>/dev/null || echo ""); do
    if [ -z "$bam" ]; then
        echo "No BAM files found for HTML reports"
        break
    fi

    echo "Processing $bam for HTML report..."
    $CONTAINER_CMD run --rm \
        --name gluetools \
        -v $PWD:/opt/bams \
        -w /opt/bams \
        --link gluetools-mysql \
        cvrbioinformatics/gluetools:latest gluetools.sh \
        --console-option log-level:FINEST \
        --inline-cmd project hcv module phdrReportingController invoke-function reportBamAsHtml ${bam} $HCVGLUE_THRESHOLD ${bam%".bam"}.html || true
done

# Cleanup: stop the MySQL container
echo "Cleaning up containers..."
$CONTAINER_CMD stop gluetools-mysql || true
$CONTAINER_CMD rm gluetools-mysql || true

echo "HCV-GLUE analysis completed"