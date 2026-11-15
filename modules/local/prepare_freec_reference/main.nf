// Chromosome length file and per-chromosome FASTA directory for Control-FREEC,
// restricted to contigs matching params.contig_regex (default: 1-22, X, Y).
process PREPARE_FREEC_REFERENCE {
    tag "${fasta.name}"
    label 'process_single'
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    publishDir "${params.outdir}/reference", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.yml' || !params.save_reference ? null : fn }

    input:
    path fasta
    path fai

    output:
    path "genome.len"   , emit: chr_len
    path "chromosomes"  , emit: chr_files
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def index_cmd = fai ? '' : "samtools faidx ${fasta}"
    """
    ${index_cmd}
    awk -v re='${params.contig_regex}' 'BEGIN{OFS="\\t"} \$1 ~ re {print \$1, \$2}' ${fasta}.fai > genome.len
    if [ ! -s genome.len ]; then
        echo "No contig in ${fasta}.fai matches --contig_regex '${params.contig_regex}'" >&2
        exit 1
    fi
    mkdir chromosomes
    cut -f1 genome.len | while read -r contig; do
        samtools faidx ${fasta} "\$contig" > "chromosomes/\${contig}.fa"
    done

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n 1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    printf 'chr1\\t248956422\\nchr2\\t242193529\\n' > genome.len
    mkdir chromosomes
    touch chromosomes/chr1.fa chromosomes/chr2.fa

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.21
    END_VERSIONS
    """
}
