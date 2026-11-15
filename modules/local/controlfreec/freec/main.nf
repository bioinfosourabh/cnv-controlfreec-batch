process CONTROLFREEC_FREEC {
    tag "${meta.id}"
    label 'process_medium'
    // Override with --freec_container, e.g. for an image that also provides bedtools.
    container "${params.freec_container}"

    publishDir "${params.outdir}/samples/${meta.id}", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.yml' ? null : fn }

    input:
    tuple val(meta), path(config), path(tumour), path(tumour_bai), path(normal), path(normal_bai)
    path chr_len
    path chr_files
    path capture_bed
    path snp_file
    path snp_tbi
    path fasta
    path fai
    path mappability

    output:
    tuple val(meta), path("${meta.id}_CNVs")        , emit: cnvs
    tuple val(meta), path("${meta.id}_ratio.txt")   , emit: ratio
    tuple val(meta), path("${meta.id}_BAF.txt")     , emit: baf, optional: true
    tuple val(meta), path("${meta.id}_info.txt")    , emit: info
    tuple val(meta), path("${meta.id}_*.BedGraph")  , emit: bedgraph, optional: true
    tuple val(meta), path("${meta.id}*.cpn")        , emit: cpn, optional: true
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def id = meta.id
    def n  = normal ? normal.name : ''
    """
    freec -conf ${config}

    # FREEC names outputs after the tumour BAM; rename them after the sample ID.
    t=${tumour.name}
    mv "\${t}_CNVs"       ${id}_CNVs
    mv "\${t}_ratio.txt"  ${id}_ratio.txt
    mv "\${t}_info.txt"   ${id}_info.txt
    [ -f "\${t}_BAF.txt" ]           && mv "\${t}_BAF.txt"           ${id}_BAF.txt
    [ -f "\${t}_ratio.BedGraph" ]    && mv "\${t}_ratio.BedGraph"    ${id}_ratio.BedGraph
    [ -f "\${t}_sample.cpn" ]        && mv "\${t}_sample.cpn"        ${id}_sample.cpn
    if [ -n "${n}" ]; then
        [ -f "${n}_control.cpn" ] && mv "${n}_control.cpn" ${id}_control.cpn
        [ -f "${n}_BAF.txt" ]     && mv "${n}_BAF.txt"     ${id}_normal_BAF.txt
    fi

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        controlfreec: \$(freec 2>&1 | sed -n 's/^Control-FREEC v\\([0-9.a-z]*\\).*/\\1/p' | head -n 1)
    END_VERSIONS
    """

    stub:
    def id = meta.id
    """
    printf '1\\t80000\\t140000\\t3\\tgain\\n' > ${id}_CNVs
    printf 'Chromosome\\tStart\\tRatio\\tMedianRatio\\tCopyNumber\\n' > ${id}_ratio.txt
    touch ${id}_info.txt ${id}_ratio.BedGraph ${id}_sample.cpn

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        controlfreec: 11.6b
    END_VERSIONS
    """
}
