#!/usr/bin/env bash
# Run Control-FREEC for every sample in a sample sheet, without Nextflow.
#
# For each sample: generate its config from the template, run freec, then
# assess_significance.R, makeGraph2.0.R, freec2bed.pl and freec2circos.pl.
# Finished steps leave a marker in <outdir>/.steps/<sample>/ and are skipped
# on the next run, so an interrupted batch can simply be started again.
# Finally all samples are merged into cohort tables.
#
# Usage: run_controlfreec_batch.sh -s samples.csv -c freec_batch.conf -o results [-j 2] [-n]
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(dirname "$script_dir")

usage() {
    cat <<USAGE
Usage: $(basename "$0") -s SAMPLES.csv -c CONFIG -o OUTDIR [options]

  -s FILE   sample sheet: sample_id,tumor_bam,normal_bam,sex[,...]   (required)
  -c FILE   batch settings, see conf/freec_batch.conf.example         (required)
  -o DIR    output directory                                           (required)
  -j N      samples run in parallel                                    [1]
  -t N      threads per FREEC run (overrides THREADS in the config)
  -n        dry run: write the configs only
  -f        force: ignore step markers and re-run everything
  -h        show this help
USAGE
}

samplesheet=""; config=""; outdir=""; jobs=1; threads=""; dry_run=false; force=false
while getopts ":s:c:o:j:t:nfh" opt; do
    case "$opt" in
        s) samplesheet="$OPTARG" ;;
        c) config="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        j) jobs="$OPTARG" ;;
        t) threads="$OPTARG" ;;
        n) dry_run=true ;;
        f) force=true ;;
        h) usage; exit 0 ;;
        :) echo "ERROR: -$OPTARG needs a value" >&2; exit 2 ;;
        *) echo "ERROR: unknown option -$OPTARG" >&2; usage >&2; exit 2 ;;
    esac
done
[[ -n "$samplesheet" && -n "$config" && -n "$outdir" ]] || { usage >&2; exit 2; }
[[ "$jobs" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: -j must be a positive integer" >&2; exit 2; }

# Settings (paths to tools and reference files) come from the config file.
FREEC_BIN="freec"
FREEC_SCRIPTS_DIR=""
TEMPLATE="${repo_dir}/assets/config_template.txt"
CHR_LEN_FILE=""
CHR_FILES=""
CAPTURE_BED=""
SNP_FILE=""
FASTA=""
MAPPABILITY=""
THREADS=4
RUN_POSTPROCESS=true
EXTRA_SET=()
# shellcheck source=/dev/null
source "$config"
[[ -n "$threads" ]] && THREADS="$threads"

mkdir -p "$outdir/configs" "$outdir/logs" "$outdir/.steps"
outdir=$(cd "$outdir" && pwd)
main_log="$outdir/logs/run_controlfreec_batch.$(date +%Y%m%d_%H%M%S).log"

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$main_log" >&2; }
die() { log "ERROR: $*"; exit 1; }

# ---------------------------------------------------------------- checks
[[ -f "$samplesheet" ]] || die "sample sheet not found: $samplesheet"
[[ -f "$TEMPLATE" ]] || die "template not found: $TEMPLATE"
[[ -f "$CHR_LEN_FILE" ]] || die "CHR_LEN_FILE not found: '$CHR_LEN_FILE'"
[[ -z "$CHR_FILES" || -d "$CHR_FILES" ]] || die "CHR_FILES is not a directory: $CHR_FILES"
for f in "$CAPTURE_BED" "$SNP_FILE" "$FASTA" "$MAPPABILITY"; do
    [[ -z "$f" || -f "$f" ]] || die "file not found: $f"
done
if [[ -n "$SNP_FILE" && -z "$FASTA" ]]; then
    die "SNP_FILE is set, so FASTA is needed for FREEC to build pileups"
fi
command -v python3 >/dev/null 2>&1 || die "python3 not found"
if [[ "$dry_run" == false ]]; then
    command -v "$FREEC_BIN" >/dev/null 2>&1 || die "FREEC binary not found: $FREEC_BIN"
    command -v samtools >/dev/null 2>&1 || die "samtools not found (FREEC reads BAM files with it)"
    if [[ -n "$SNP_FILE" && -n "$CAPTURE_BED" ]] && ! command -v bedtools >/dev/null 2>&1; then
        die "bedtools not found: FREEC needs it for BAF on targeted data"
    fi
fi

find_script() {
    local name="$1"
    if [[ -n "$FREEC_SCRIPTS_DIR" && -f "$FREEC_SCRIPTS_DIR/$name" ]]; then
        echo "$FREEC_SCRIPTS_DIR/$name"
    elif command -v "$name" >/dev/null 2>&1; then
        command -v "$name"
    fi
}
assess_r=$(find_script assess_significance.R)
graph_r=$(find_script makeGraph2.0.R)
bed_pl=$(find_script freec2bed.pl)
circos_pl=$(find_script freec2circos.pl)
if [[ "$RUN_POSTPROCESS" == true && "$dry_run" == false ]]; then
    if ! command -v R >/dev/null 2>&1; then
        log "WARNING: R not found; significance tests and plots will be skipped"
        assess_r=""; graph_r=""
    elif [[ -z "$assess_r" ]]; then
        log "WARNING: assess_significance.R not found; set FREEC_SCRIPTS_DIR"
    fi
fi

# ---------------------------------------------------------------- configs
set_args=()
for kv in "${EXTRA_SET[@]}"; do set_args+=(--set "$kv"); done
log "writing configs to $outdir/configs"
python3 "$script_dir/make_freec_configs.py" \
    --samplesheet "$samplesheet" --template "$TEMPLATE" --outdir "$outdir/configs" \
    --chr-len "$CHR_LEN_FILE" ${CHR_FILES:+--chr-files "$CHR_FILES"} \
    ${CAPTURE_BED:+--capture-bed "$CAPTURE_BED"} ${SNP_FILE:+--snp-file "$SNP_FILE"} \
    ${FASTA:+--fasta "$FASTA"} ${MAPPABILITY:+--mappability "$MAPPABILITY"} \
    --threads "$THREADS" --output-dir "$outdir/samples/{sample_id}" --check-files \
    "${set_args[@]}" | tee "$outdir/configs/manifest.tsv" >&2

if [[ "$dry_run" == true ]]; then
    log "dry run: configs written, nothing else run"
    exit 0
fi

# ---------------------------------------------------------------- per sample
# step NAME SAMPLE CMD...  runs CMD unless the marker exists; logs to logs/SAMPLE/NAME.log
step() {
    local name="$1" sample="$2"; shift 2
    local marker="$outdir/.steps/$sample/$name.done" slog="$outdir/logs/$sample/$name.log"
    mkdir -p "$outdir/.steps/$sample" "$outdir/logs/$sample"
    if [[ -f "$marker" && "$force" == false ]]; then
        log "$sample: skip $name (done)"
        return 0
    fi
    log "$sample: $name"
    if "$@" > "$slog" 2>&1; then
        touch "$marker"
    else
        log "$sample: $name FAILED, last lines of $slog:"
        tail -n 5 "$slog" | sed 's/^/    /' >&2
        return 1
    fi
}

# Post-processing helpers; FREEC's R scripts read the script from stdin.
assess_significance() { R --slave --args "$2" "$3" < "$1"; }
make_graph() { (cd "$1" && R --slave --args "$3" ${4:+"$4"} < "$2"); }
freec_perl() { perl "$1" -f "$2" -p "$3" > "$4"; }

run_sample() {
    local sample="$1"
    local cfg="$outdir/configs/${sample}.config.txt"
    local sdir="$outdir/samples/$sample"
    mkdir -p "$sdir"
    # A changed config invalidates the previous FREEC run.
    if [[ -f "$outdir/.steps/$sample/freec.done" ]] && ! cmp -s "$cfg" "$outdir/.steps/$sample/config.used"; then
        log "$sample: config changed, re-running"
        rm -f "$outdir/.steps/$sample/"*.done
    fi
    step freec "$sample" "$FREEC_BIN" -conf "$cfg" || return 1
    cp "$cfg" "$outdir/.steps/$sample/config.used"

    local cnvs ratio baf
    cnvs=$(find "$sdir" -maxdepth 1 -name '*_CNVs' | head -n 1)
    ratio=$(find "$sdir" -maxdepth 1 -name '*_ratio.txt' | head -n 1)
    baf=$(find "$sdir" -maxdepth 1 -name '*_BAF.txt' ! -name '*_normal_BAF.txt' | head -n 1)
    [[ -n "$cnvs" && -n "$ratio" ]] || { log "$sample: FREEC produced no _CNVs/_ratio.txt"; return 1; }
    [[ "$RUN_POSTPROCESS" == true ]] || return 0

    local ploidy
    ploidy=$(awk -F' *= *' '$1 == "ploidy" {print $2}' "$cfg")
    if [[ -n "$assess_r" ]]; then
        step significance "$sample" assess_significance "$assess_r" "$cnvs" "$ratio" || true
    fi
    if [[ -n "$graph_r" ]]; then
        step graph "$sample" make_graph "$sdir" "$graph_r" "$ratio" "$baf" || true
    fi
    if [[ -n "$bed_pl" ]]; then
        step bed "$sample" freec_perl "$bed_pl" "$ratio" "${ploidy:-2}" "$sdir/${sample}.CNV.bed" || true
    fi
    if [[ -n "$circos_pl" ]]; then
        step circos "$sample" freec_perl "$circos_pl" "$ratio" "${ploidy:-2}" "$sdir/${sample}.circos.txt" || true
    fi
}

mapfile -t samples < <(cut -f1 "$outdir/configs/manifest.tsv")
log "${#samples[@]} sample(s), $jobs in parallel"

failed=()
declare -A pid_sample=()
for s in "${samples[@]}"; do
    while (( $(jobs -rp | wc -l) >= jobs )); do sleep 1; done
    run_sample "$s" &
    pid_sample[$!]="$s"
done
for pid in "${!pid_sample[@]}"; do
    if ! wait "$pid"; then failed+=("${pid_sample[$pid]}"); fi
done

# ---------------------------------------------------------------- cohort
mapfile -t cnv_files < <(find "$outdir/samples" -maxdepth 2 \( -name '*_CNVs' -o -name '*_CNVs.p.value.txt' \) | sort)
if (( ${#cnv_files[@]} > 0 )); then
    python3 "$script_dir/freec_cohort_summary.py" --name-from-dir --chr-len "$CHR_LEN_FILE" --outdir "$outdir/cohort" "${cnv_files[@]}"
    log "cohort tables written to $outdir/cohort"
fi

if (( ${#failed[@]} > 0 )); then
    log "finished with ${#failed[@]} failed sample(s): ${failed[*]}"
    exit 1
fi
log "finished: ${#samples[@]} sample(s) completed"
