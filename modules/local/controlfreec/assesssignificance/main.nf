// Wilcoxon rank-sum and Kolmogorov-Smirnov p-values for each predicted CNV.
process CONTROLFREEC_ASSESSSIGNIFICANCE {
    tag "${meta.id}"
    label 'process_low'
    // Same image as nf-core/modules uses for this script.
    container 'quay.io/biocontainers/control-freec:11.6--h1b792b2_1'

    publishDir "${params.outdir}/samples/${meta.id}", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.yml' ? null : fn }

    input:
    tuple val(meta), path(cnvs), path(ratio)

    output:
    tuple val(meta), path("*.p.value.txt"), emit: p_value
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat \$(which assess_significance.R) | R --slave --args ${cnvs} ${ratio}

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        controlfreec: 11.6
        r-base: \$(R --version | head -n 1 | sed 's/R version \\([0-9.]*\\).*/\\1/')
    END_VERSIONS
    """

    stub:
    """
    printf 'chr\\tstart\\tend\\tcopy number\\tstatus\\tWilcoxonRankSumTestPvalue\\tKolmogorovSmirnovPvalue\\n' > ${cnvs}.p.value.txt

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        controlfreec: 11.6
        r-base: 4.0
    END_VERSIONS
    """
}
