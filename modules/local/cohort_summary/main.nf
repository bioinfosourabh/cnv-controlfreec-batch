process COHORT_SUMMARY {
    label 'process_single'
    container 'quay.io/biocontainers/python:3.9--1'

    publishDir "${params.outdir}/cohort", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.yml' ? null : fn }

    input:
    path cnv_files, stageAs: "calls/*"
    path chr_len

    output:
    path "cohort_cnv_segments.tsv"   , emit: segments
    path "cohort_cnv_summary.tsv"    , emit: summary
    path "freec_cnv_summary_mqc.tsv" , emit: mqc
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    freec_cohort_summary.py --chr-len ${chr_len} --alpha ${params.alpha} --outdir . calls/*

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    freec_cohort_summary.py --chr-len ${chr_len} --outdir . calls/*

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.9.1
    END_VERSIONS
    """
}
