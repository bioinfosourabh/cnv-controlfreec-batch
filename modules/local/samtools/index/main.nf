process SAMTOOLS_INDEX {
    tag "${bam.name}"
    label 'process_low'
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    tuple val(key), path(bam)

    output:
    tuple val(key), path(bam), path("*.bai"), emit: bam
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools index -@ ${task.cpus} ${bam}

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n 1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${bam}.bai

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.21
    END_VERSIONS
    """
}
