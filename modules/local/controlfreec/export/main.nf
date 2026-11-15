// Copy-number tracks for genome browsers (BED) and for Circos.
process CONTROLFREEC_EXPORT {
    tag "${meta.id}"
    label 'process_single'
    container 'quay.io/biocontainers/control-freec:11.6b--hde5307d_3'

    publishDir "${params.outdir}/samples/${meta.id}", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.yml' ? null : fn }

    input:
    tuple val(meta), path(ratio)

    output:
    tuple val(meta), path("*.CNV.bed")   , emit: bed
    tuple val(meta), path("*.circos.txt"), emit: circos
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ploidy = meta.ploidy ?: 2
    """
    freec2bed.pl -f ${ratio} -p ${ploidy} > ${meta.id}.CNV.bed
    freec2circos.pl -f ${ratio} -p ${ploidy} > ${meta.id}.circos.txt

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        controlfreec: 11.6b
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.CNV.bed ${meta.id}.circos.txt

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        controlfreec: 11.6b
    END_VERSIONS
    """
}
