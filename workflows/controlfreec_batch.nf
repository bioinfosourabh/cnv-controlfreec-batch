//
// Control-FREEC for a batch of samples: one config per sample generated from
// the template, FREEC, significance testing, plots, BED/Circos export and a
// merged cohort table.
//
include { readSamplesheet; findBai           } from '../subworkflows/local/utils_pipeline'
include { PREPARE_FREEC_REFERENCE            } from '../modules/local/prepare_freec_reference'
include { SAMTOOLS_INDEX                     } from '../modules/local/samtools/index'
include { MAKE_FREEC_CONFIG                  } from '../modules/local/make_freec_config'
include { CONTROLFREEC_FREEC                 } from '../modules/local/controlfreec/freec'
include { CONTROLFREEC_ASSESSSIGNIFICANCE    } from '../modules/local/controlfreec/assesssignificance'
include { CONTROLFREEC_MAKEGRAPH             } from '../modules/local/controlfreec/makegraph'
include { CONTROLFREEC_EXPORT                } from '../modules/local/controlfreec/export'
include { COHORT_SUMMARY                     } from '../modules/local/cohort_summary'
include { MULTIQC                            } from '../modules/local/multiqc'

def optionalFile(path) {
    return path ? file(path, checkIfExists: true) : []
}

workflow CONTROLFREEC_BATCH {
    main:
    ch_versions = channel.empty()

    def samples = readSamplesheet(params.input)

    //
    // Reference files
    //
    ch_fasta = channel.value(file(params.fasta, checkIfExists: true))
    def fai_path = params.fasta_fai ?: "${params.fasta}.fai"
    def fai = file(fai_path).exists() ? file(fai_path) : []
    ch_fai = channel.value(fai)

    ch_chr_len   = channel.empty()
    ch_chr_files = channel.empty()
    if (params.chr_len && params.chr_files) {
        ch_chr_len   = channel.value(file(params.chr_len, checkIfExists: true))
        ch_chr_files = channel.value(file(params.chr_files, checkIfExists: true))
    } else {
        PREPARE_FREEC_REFERENCE(ch_fasta, ch_fai)
        ch_chr_len   = params.chr_len ? channel.value(file(params.chr_len, checkIfExists: true)) : PREPARE_FREEC_REFERENCE.out.chr_len
        ch_chr_files = PREPARE_FREEC_REFERENCE.out.chr_files
        ch_versions  = ch_versions.mix(PREPARE_FREEC_REFERENCE.out.versions)
    }
    def capture     = optionalFile(params.capture_bed)
    def snp_vcf     = optionalFile(params.snp_vcf)
    def snp_tbi     = params.snp_vcf && file("${params.snp_vcf}.tbi").exists() ? file("${params.snp_vcf}.tbi") : []
    def mappability = optionalFile(params.mappability)

    // Names under which the FREEC task stages the reference files.
    def ref_names = [
        chr_len    : params.chr_len ? file(params.chr_len).name : 'genome.len',
        chr_files  : params.chr_files ? file(params.chr_files).name : 'chromosomes',
        capture_bed: capture ? capture.name : '',
        snp_file   : snp_vcf ? snp_vcf.name : '',
        fasta      : snp_vcf ? file(params.fasta).name : '',
        mappability: mappability ? mappability.name : ''
    ]

    //
    // BAM indexes: use existing ones, build missing ones
    //
    def all_bams = samples.collect { s -> [s[2]] + (s[3] ? [s[3]] : []) }.flatten().unique { b -> b.toString() }
    ch_bams = channel.fromList(all_bams).map { bam -> [bam.toString(), bam, findBai(bam)] }
        .branch { _key, _bam, bai ->
            indexed: bai
            missing: true
        }
    SAMTOOLS_INDEX(ch_bams.missing.map { key, bam, _bai -> [key, bam] })
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
    ch_indexed = ch_bams.indexed.mix(SAMTOOLS_INDEX.out.bam)   // [key, bam, bai]

    //
    // Per-sample config from the template
    //
    ch_samples = channel.fromList(samples)
    MAKE_FREEC_CONFIG(
        ch_samples.map { meta, row, _t, _n -> [meta, row] },
        file(params.template, checkIfExists: true),
        ref_names
    )
    ch_versions = ch_versions.mix(MAKE_FREEC_CONFIG.out.versions.first())

    // [meta, config, tumour, tumour_bai, normal, normal_bai]
    ch_freec_in = ch_samples
        .map { meta, _row, t, n -> [t.toString(), meta, n ? n.toString() : 'NO_NORMAL'] }
        .combine(ch_indexed, by: 0)
        .map { _tkey, meta, nkey, tbam, tbai -> [nkey, meta, tbam, tbai] }
        .branch { nkey, _meta, _tbam, _tbai ->
            paired: nkey != 'NO_NORMAL'
            tumour_only: true
        }
    ch_paired = ch_freec_in.paired
        .combine(ch_indexed, by: 0)
        .map { _nkey, meta, tbam, tbai, nbam, nbai -> [meta, tbam, tbai, nbam, nbai] }
    ch_tonly = ch_freec_in.tumour_only
        .map { _nkey, meta, tbam, tbai -> [meta, tbam, tbai, [], []] }
    ch_run = MAKE_FREEC_CONFIG.out.config
        .join(ch_paired.mix(ch_tonly))

    CONTROLFREEC_FREEC(
        ch_run, ch_chr_len, ch_chr_files,
        capture, snp_vcf, snp_tbi, ch_fasta, ch_fai, mappability
    )
    ch_versions = ch_versions.mix(CONTROLFREEC_FREEC.out.versions.first())

    //
    // Post-processing
    //
    ch_calls = CONTROLFREEC_FREEC.out.cnvs.map { _meta, f -> f }
    if (params.run_significance) {
        CONTROLFREEC_ASSESSSIGNIFICANCE(CONTROLFREEC_FREEC.out.cnvs.join(CONTROLFREEC_FREEC.out.ratio))
        ch_calls    = ch_calls.mix(CONTROLFREEC_ASSESSSIGNIFICANCE.out.p_value.map { _meta, f -> f })
        ch_versions = ch_versions.mix(CONTROLFREEC_ASSESSSIGNIFICANCE.out.versions.first())
    }
    if (params.run_plots) {
        ch_graph_in = CONTROLFREEC_FREEC.out.ratio
            .join(CONTROLFREEC_FREEC.out.baf, remainder: true)
            .map { meta, ratio, baf -> [meta, ratio, baf ?: []] }
        CONTROLFREEC_MAKEGRAPH(ch_graph_in)
        ch_versions = ch_versions.mix(CONTROLFREEC_MAKEGRAPH.out.versions.first())
    }
    if (params.run_export) {
        CONTROLFREEC_EXPORT(CONTROLFREEC_FREEC.out.ratio)
        ch_versions = ch_versions.mix(CONTROLFREEC_EXPORT.out.versions.first())
    }

    COHORT_SUMMARY(ch_calls.collect(sort: { a, b -> a.name <=> b.name }), ch_chr_len)
    ch_versions = ch_versions.mix(COHORT_SUMMARY.out.versions)

    ch_versions_yml = ch_versions
        .unique { f -> f.text }
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info", sort: true, newLine: true)

    if (!params.skip_multiqc) {
        MULTIQC(COHORT_SUMMARY.out.mqc.collect(), file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true))
    }

    emit:
    summary  = COHORT_SUMMARY.out.summary
    versions = ch_versions_yml
}
