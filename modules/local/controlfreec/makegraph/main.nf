// Genome-wide ratio (and BAF) plots with FREEC's makeGraph2.0.R.
process CONTROLFREEC_MAKEGRAPH {
    tag "${meta.id}"
    label 'process_low'
    container 'quay.io/biocontainers/control-freec:11.6b--hde5307d_3'

    publishDir "${params.outdir}/samples/${meta.id}/plots", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.yml' ? null : fn }

    input:
    tuple val(meta), path(ratio), path(baf)

    output:
    tuple val(meta), path("*.png"), emit: png
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def baf_arg = baf ? "${baf}" : ''
    """
    cat \$(which makeGraph2.0.R) | R --slave --args ${ratio} ${baf_arg}

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        controlfreec: 11.6b
    END_VERSIONS
    """

    stub:
    """
    touch ${ratio}.png ${ratio}.log2.png

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        controlfreec: 11.6b
    END_VERSIONS
    """
}
