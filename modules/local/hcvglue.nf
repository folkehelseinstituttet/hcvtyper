process HCVGLUE {

    label 'process_low'

    // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'docker.io/docker:24.0.7-cli-alpine3.18' }"

    stageInMode = 'copy' // Can't mount symlinked files into docker containers

    input:
    path '*'

    output:
    path("*.json")     , optional: true, emit: GLUE_json
    path("*.html")     , optional: true, emit: GLUE_html
    path "versions.yml", emit: versions

    script:
    """
    # Remove the container in case it is already running
    if docker ps -a --filter "name=gluetools-mysql" --format '{{.Names}}' | grep -q "^gluetools-mysql\$"; then
        echo "Container 'gluetools-mysql' is running or exists."
        # Stop the container
        docker stop gluetools-mysql
        # Remove the container
        docker rm gluetools-mysql
        echo "Container 'gluetools-mysql' has been stopped and removed."
    else
        echo "Container 'gluetools-mysql' is not running or does not exist."
    fi

    # Pull the latest images
    docker pull cvrbioinformatics/gluetools-mysql:latest
    docker pull cvrbioinformatics/gluetools:latest

    # Start the gluetools-mysql containter
    docker run --detach --name gluetools-mysql cvrbioinformatics/gluetools-mysql:latest

    # Install the pre-built GLUE HCV project
    TIMEOUT=300
    START_TIME=\$(date +%s)

    until docker exec gluetools-mysql mysql --user=root --password=root123 -e "status" &> /dev/null
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
    docker exec gluetools-mysql installGlueProject.sh ncbi_hcv_glue

    # Make a for loop over all bam files and run HCV-GLUE
    ## Adding || true to the end of the command to prevent the pipeline from failing if the bam file is not valid

    for bam in \$(ls *.bam)
    do
    docker run --rm \
       --name gluetools \
        -v \$PWD:/opt/bams \
        -w /opt/bams \
        --link gluetools-mysql \
        cvrbioinformatics/gluetools:latest gluetools.sh \
         -p cmd-result-format:json \
        -EC \
        -i project hcv module phdrReportingController invoke-function reportBam \${bam} 15.0 > \${bam%".bam"}.json || true
    done

    # Then create html files
    for bam in \$(ls *.bam)
    do
    docker run --rm \
        --name gluetools \
        -v \$PWD:/opt/bams \
        -w /opt/bams \
        --link gluetools-mysql \
        cvrbioinformatics/gluetools:latest gluetools.sh \
    	--console-option log-level:FINEST \
        --inline-cmd project hcv module phdrReportingController invoke-function reportBamAsHtml \${bam} 15.0 \${bam%".bam"}.html || true
    done


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GLUE project version: \$(ls *.json | head -n 1 | xargs -I {} grep -oP '"projectVersion"\\s*:\\s*"\\K[^"]+' {})
        GLUE engine version: \$(ls *.json | head -n 1 | xargs -I {} grep -oP '"engineVersion"\\s*:\\s*"\\K[^"]+' {})
        GLUE extension version: \$(ls *.json | head -n 1 | xargs -I {} grep -oP '"extensionVersion"\\s*:\\s*"\\K[^"]+' {})
    END_VERSIONS
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
