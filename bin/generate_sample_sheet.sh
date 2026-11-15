#!/usr/bin/env bash
# Build a Control-FREEC sample sheet from a directory of BAM files by pairing
# tumour and normal BAMs that share a sample prefix, e.g.
#   S12_SO.bam (tumour) + S12_GL.bam (normal)  ->  S12,<tumour>,<normal>,<sex>
#
# Usage: generate_sample_sheet.sh -d BAM_DIR [-o samples.csv] [-t _SO.bam]
#                                 [-n _GL.bam] [-s auto|XX|XY] [-y 0.0005] [-a]
set -euo pipefail

usage() {
    cat <<USAGE
Usage: $(basename "$0") -d BAM_DIR [options]

  -d DIR     directory searched (recursively) for BAM files        (required)
  -o FILE    output CSV                                            [samples.csv]
  -t SUFFIX  tumour BAM suffix                                     [_SO.bam]
  -n SUFFIX  normal BAM suffix                                     [_GL.bam]
  -s SEX     auto, XX or XY; auto uses the chrY read fraction      [auto]
  -y FRAC    chrY fraction of mapped reads at or above which a
             sample is called XY in auto mode                      [0.0005]
  -a         also list tumours without a normal (tumour-only rows)
  -h         show this help

Sex is inferred from the normal BAM when present, otherwise from the tumour.
The chrY heuristic depends on the assay; check the result, and set -s or edit
the sheet when in doubt.
USAGE
}

bam_dir=""
out="samples.csv"
tumour_suffix="_SO.bam"
normal_suffix="_GL.bam"
sex_mode="auto"
y_threshold="0.0005"
allow_tumour_only=false

while getopts ":d:o:t:n:s:y:ah" opt; do
    case "$opt" in
        d) bam_dir="$OPTARG" ;;
        o) out="$OPTARG" ;;
        t) tumour_suffix="$OPTARG" ;;
        n) normal_suffix="$OPTARG" ;;
        s) sex_mode="$OPTARG" ;;
        y) y_threshold="$OPTARG" ;;
        a) allow_tumour_only=true ;;
        h) usage; exit 0 ;;
        :) echo "ERROR: -$OPTARG needs a value" >&2; usage >&2; exit 2 ;;
        *) echo "ERROR: unknown option -$OPTARG" >&2; usage >&2; exit 2 ;;
    esac
done

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }

[[ -n "$bam_dir" ]] || { usage >&2; exit 2; }
[[ -d "$bam_dir" ]] || { log "ERROR: not a directory: $bam_dir"; exit 1; }
case "$sex_mode" in auto|XX|XY) ;; *) log "ERROR: -s must be auto, XX or XY"; exit 2 ;; esac
if [[ "$sex_mode" == auto ]] && ! command -v samtools >/dev/null 2>&1; then
    log "ERROR: samtools is needed for -s auto (or set -s XX / -s XY)"
    exit 1
fi

infer_sex() {
    local bam="$1" stats
    if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
        log "indexing $bam"
        samtools index "$bam"
    fi
    stats=$(samtools idxstats "$bam")
    awk -v thr="$y_threshold" '
        { total += $3 }
        $1 == "chrY" || $1 == "Y" { y += $3 }
        END {
            frac = (total > 0) ? y / total : 0
            print (frac >= thr) ? "XY" : "XX"
        }' <<<"$stats"
}

n_pairs=0
n_tumour_only=0
n_skipped=0
tmp=$(mktemp)
trap 'rm -f "$tmp"' EXIT

while IFS= read -r -d '' tumour; do
    base=$(basename "$tumour" "$tumour_suffix")
    normal="$(dirname "$tumour")/${base}${normal_suffix}"
    tumour_abs=$(realpath "$tumour")
    if [[ -f "$normal" ]]; then
        normal_abs=$(realpath "$normal")
        sex_source="$normal_abs"
        n_pairs=$((n_pairs + 1))
        log "pair: $base"
    elif [[ "$allow_tumour_only" == true ]]; then
        normal_abs=""
        sex_source="$tumour_abs"
        n_tumour_only=$((n_tumour_only + 1))
        log "tumour only: $base"
    else
        log "skipped (no ${base}${normal_suffix}; use -a for tumour-only): $base"
        n_skipped=$((n_skipped + 1))
        continue
    fi
    if [[ "$sex_mode" == auto ]]; then
        sex=$(infer_sex "$sex_source")
    else
        sex="$sex_mode"
    fi
    printf '%s,%s,%s,%s\n' "$base" "$tumour_abs" "$normal_abs" "$sex" >> "$tmp"
done < <(find "$bam_dir" -type f -name "*${tumour_suffix}" -print0 | sort -z)

if [[ ! -s "$tmp" ]]; then
    log "ERROR: no tumour BAMs matching *${tumour_suffix} were usable in $bam_dir"
    exit 1
fi

{ echo "sample_id,tumor_bam,normal_bam,sex"; sort -t, -k1,1 "$tmp"; } > "$out"
log "wrote $out: $n_pairs tumour-normal pair(s), $n_tumour_only tumour-only, $n_skipped skipped"
