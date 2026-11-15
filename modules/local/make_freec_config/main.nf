// Writes the Control-FREEC config for one sample with bin/make_freec_configs.py.
// File references in the config are the names under which the FREEC task
// stages its inputs, so the config is valid inside that task's directory.
process MAKE_FREEC_CONFIG {
    tag "${meta.id}"
    label 'process_single'
    container 'quay.io/biocontainers/python:3.9--1'

    publishDir "${params.outdir}/configs", mode: params.publish_dir_mode,
        saveAs: { fn -> fn.endsWith('.config.txt') ? fn : null }

    input:
    tuple val(meta), val(row)   // row: sample-sheet columns with BAM paths replaced by file names
    path template
    val  ref_names              // map of reference-file names, see workflows/controlfreec_batch.nf

    output:
    tuple val(meta), path("*.config.txt"), emit: config
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def header = row.keySet().join(',')
    def values = row.values().join(',')
    def optional = [
        '--capture-bed': ref_names.capture_bed,
        '--snp-file'   : ref_names.snp_file,
        '--fasta'      : ref_names.fasta,
        '--mappability': ref_names.mappability
    ].findAll { _flag, v -> v }.collect { flag, v -> "${flag} ${v}" }.join(' ')
    def sets = (params.freec_set ? params.freec_set.tokenize(',') : []).collect { kv -> "--set '${kv.trim()}'" }.join(' ')
    """
    printf '%s\\n%s\\n' '${header}' '${values}' > sheet.csv

    make_freec_configs.py \\
        --samplesheet sheet.csv \\
        --template ${template} \\
        --outdir . \\
        --chr-len ${ref_names.chr_len} \\
        --chr-files ${ref_names.chr_files} \\
        ${optional} \\
        --threads ${params.freec_threads} \\
        --output-dir . \\
        ${sets}

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    // The generator is plain Python and fast, so stub runs execute it for real.
    def header_s = row.keySet().join(',')
    def values_s = row.values().join(',')
    def optional_s = [
        '--capture-bed': ref_names.capture_bed,
        '--snp-file'   : ref_names.snp_file,
        '--fasta'      : ref_names.fasta,
        '--mappability': ref_names.mappability
    ].findAll { _flag, v -> v }.collect { flag, v -> "${flag} ${v}" }.join(' ')
    def sets_s = (params.freec_set ? params.freec_set.tokenize(',') : []).collect { kv -> "--set '${kv.trim()}'" }.join(' ')
    """
    printf '%s\\n%s\\n' '${header_s}' '${values_s}' > sheet.csv

    make_freec_configs.py \\
        --samplesheet sheet.csv \\
        --template ${template} \\
        --outdir . \\
        --chr-len ${ref_names.chr_len} \\
        --chr-files ${ref_names.chr_files} \\
        ${optional_s} \\
        --threads ${params.freec_threads} \\
        --output-dir . \\
        ${sets_s}

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

}
