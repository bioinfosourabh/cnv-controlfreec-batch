//
// Sample-sheet parsing, parameter validation and run summary.
//

// Sample sheet: sample_id,tumor_bam,normal_bam,sex plus optional per-sample
// columns that override template placeholders (e.g. ploidy, contamination).
def readSamplesheet(String path) {
    def csv = file(path, checkIfExists: true)
    def rows = csv.splitCsv(header: true, strip: true)
    if (!rows) {
        error("Sample sheet ${path} has no rows")
    }
    def required = ['sample_id', 'tumor_bam', 'normal_bam', 'sex']
    def missing = required.findAll { c -> !rows[0].containsKey(c) }
    if (missing) {
        error("Sample sheet ${path} is missing column(s): ${missing.join(', ')}")
    }
    def seen = [] as Set
    return rows.withIndex().collect { row, i ->
        def line = i + 2
        def id = row.sample_id
        if (!(id ==~ /^[A-Za-z0-9][A-Za-z0-9._-]*$/)) {
            error("Sample sheet line ${line}: invalid sample_id '${id}'")
        }
        if (!seen.add(id)) {
            error("Sample sheet line ${line}: duplicate sample_id '${id}'")
        }
        def sex = (row.sex ?: '').toUpperCase()
        if (!(sex in ['XX', 'XY'])) {
            error("Sample sheet line ${line}: sex must be XX or XY, got '${row.sex}'")
        }
        row.values().each { v ->
            if (v.contains(',') || v.contains("'")) {
                error("Sample sheet line ${line}: values must not contain commas or quotes")
            }
        }
        def tumour = resolveBam(row.tumor_bam, csv, line, 'tumor_bam')
        def normal = row.normal_bam ? resolveBam(row.normal_bam, csv, line, 'normal_bam') : null
        def extra = row.findAll { k, _v -> !(k in required) }
        def meta = [id: id, sex: sex, paired: normal != null, ploidy: row.ploidy ?: null]
        def sheet_row = [sample_id: id, tumor_bam: tumour.name, normal_bam: normal ? normal.name : '', sex: sex] + extra
        [meta, sheet_row, tumour, normal]
    }
}

def resolveBam(String value, csv, int line, String column) {
    if (!value) {
        error("Sample sheet line ${line}: ${column} is empty")
    }
    if (!value.endsWith('.bam')) {
        error("Sample sheet line ${line}: ${column} must be a .bam file: ${value}")
    }
    def f = (value ==~ /^[a-zA-Z][a-zA-Z0-9+.-]*:\/\/.*/ || value.startsWith('/')) ? file(value) : csv.parent.resolve(value)
    if (!f.exists()) {
        error("Sample sheet line ${line}: ${column} not found: ${f}")
    }
    return f
}

// BAM index next to the BAM (<x>.bam.bai or <x>.bai), or null.
def findBai(bam) {
    def candidates = [file("${bam}.bai"), file(bam.toString().replaceAll(/\.bam$/, '.bai'))]
    return candidates.find { c -> c.exists() }
}

def validateParams(String schema_path) {
    def errors = []
    def schema = new groovy.json.JsonSlurper().parse(new File(schema_path))
    def props = [:]
    schema['$defs'].each { _group, definition -> props.putAll(definition.properties) }
    params.each { name, _value ->
        if (!props.containsKey(name) && !name.toString().contains('-')) {
            log.warn("Parameter --${name} is not defined in nextflow_schema.json and will be ignored")
        }
    }
    props.each { name, spec ->
        def value = params[name]
        if (value == null) {
            return
        }
        def allowed = spec['enum']
        if (allowed && !(value in allowed)) {
            errors.add("--${name} must be one of ${allowed.join(', ')} (got '${value}')")
        }
        if (spec.type == 'integer' && !(value instanceof Integer)) {
            errors.add("--${name} must be an integer (got '${value}')")
        }
        if (spec.type == 'boolean' && !(value instanceof Boolean)) {
            errors.add("--${name} must be true or false (got '${value}')")
        }
    }
    ['input', 'outdir', 'fasta'].each { p ->
        if (!params[p]) {
            errors.add("--${p} is required")
        }
    }
    if (params.freec_set) {
        params.freec_set.tokenize(',').each { kv ->
            if (!(kv.trim() ==~ /^[A-Za-z0-9_]+=.*$/)) {
                errors.add("--freec_set expects NAME=VALUE[,NAME=VALUE...], got '${kv}'")
            }
        }
    }
    if (params.snp_vcf && params.capture_bed && params.freec_container == 'quay.io/biocontainers/control-freec:11.6b--hde5307d_3') {
        errors.add("BAF on targeted data (--snp_vcf with --capture_bed) needs bedtools inside the FREEC container, " +
            "which the default BioContainers image lacks. Build containers/freec-bedtools/Dockerfile and pass --freec_container, " +
            "or leave out --snp_vcf.")
    }
    if (errors) {
        error("Invalid parameters:\n  - " + errors.join("\n  - "))
    }
}

def pipelineSummary() {
    def keys = ['input', 'outdir', 'fasta', 'chr_len', 'capture_bed', 'snp_vcf', 'freec_set', 'freec_container']
    def width = keys.collect { k -> k.size() }.max()
    def lines = keys.collect { k -> "  ${k.padRight(width)} : ${params[k] ?: '-'}" }
    def rule = '-' * 60
    def mode = params.capture_bed ? 'targeted (capture BED)' : 'whole genome'
    return (["${workflow.manifest.name} v${workflow.manifest.version}", rule] + lines + ["  mode : ${mode}", rule]).join('\n')
}
