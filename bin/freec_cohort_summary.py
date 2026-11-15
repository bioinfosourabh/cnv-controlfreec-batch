#!/usr/bin/env python3
"""Merge Control-FREEC CNV calls from many samples into cohort tables.

Inputs are FREEC `*_CNVs` files or the `*_CNVs.p.value.txt` files written by
assess_significance.R (which add Wilcoxon and Kolmogorov-Smirnov p-values).
When both exist for a sample, the p-value file is used.

Outputs (in --outdir):
  cohort_cnv_segments.tsv   one row per CNV segment, all samples
  cohort_cnv_summary.tsv    one row per sample: counts, altered bases, FGA
  freec_cnv_summary_mqc.tsv the per-sample table as MultiQC custom content

    freec_cohort_summary.py --chr-len hg38.len --outdir summary results/*/*_CNVs*
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys

SEGMENT_COLUMNS = ["sample", "chrom", "start", "end", "length", "copy_number", "status",
                   "genotype", "uncertainty", "wilcoxon_p", "ks_p", "significant"]
SUMMARY_COLUMNS = ["sample", "n_segments", "n_gain", "n_loss", "n_cn_loh", "bp_gain", "bp_loss",
                   "bp_altered", "fraction_genome_altered", "n_significant"]


def sample_name(path: str) -> str:
    name = os.path.basename(path)
    name = re.sub(r"\.p\.value\.txt$", "", name)
    name = re.sub(r"_CNVs$", "", name)
    name = re.sub(r"\.(bam|cram|pileup)$", "", name)
    return name


def read_genome_size(chr_len: str | None) -> int | None:
    if not chr_len:
        return None
    total = 0
    with open(chr_len) as fh:
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            # FREEC accepts "chr<TAB>len" and "index<TAB>chr<TAB>len".
            total += int(parts[-1])
    return total


def _float(v: str) -> float | None:
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def parse_cnv_file(path: str, sample: str, alpha: float) -> list[dict]:
    rows = []
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]
    if not lines:
        return rows
    header = None
    first = lines[0].split("\t")
    if first[0].lower() in {"chr", "chromosome"}:
        header = [h.strip().lower().replace(" ", "_") for h in first]
        lines = lines[1:]
    for line in lines:
        f = line.split("\t")
        rec = dict(zip(header, f)) if header else {}
        chrom = rec.get("chr", f[0])
        start, end = int(rec.get("start", f[1])), int(rec.get("end", f[2]))
        if end <= start:
            # FREEC can report zero-length records at contig ends; they carry no CNV.
            continue
        cn = rec.get("copy_number", f[3])
        status = rec.get("status", f[4])
        genotype = rec.get("genotype", f[5] if len(f) > 6 and not header else "")
        uncertainty = rec.get("uncertainty", f[6] if len(f) > 6 and not header else "")
        w = rec.get("wilcoxonranksumtestpvalue", "")
        ks = rec.get("kolmogorovsmirnovpvalue", "")
        wp, kp = _float(w), _float(ks)
        significant = "" if wp is None and kp is None else str(
            (wp is not None and wp < alpha) and (kp is not None and kp < alpha)).lower()
        rows.append({
            "sample": sample, "chrom": chrom, "start": start, "end": end,
            "length": max(0, end - start), "copy_number": cn, "status": status,
            "genotype": genotype if genotype not in ("-", "") else "",
            "uncertainty": uncertainty if uncertainty not in ("-1", "-") else "",
            "wilcoxon_p": w, "ks_p": ks, "significant": significant,
        })
    return rows


def is_cn_loh(seg: dict) -> bool:
    """Copy-neutral LOH: normal copy number but a homozygous genotype (e.g. AA)."""
    g = seg["genotype"]
    return status_of(seg) == "normal" and len(g) > 1 and len(set(g)) == 1


def status_of(seg: dict) -> str:
    s = (seg["status"] or "").lower()
    if s in {"gain", "loss"}:
        return s
    return "normal" if s in {"normal", "neutral", ""} else s


def summarise(segments: list[dict], samples: list[str], genome_size: int | None) -> list[dict]:
    out = []
    for s in samples:
        segs = [x for x in segments if x["sample"] == s]
        gains = [x for x in segs if status_of(x) == "gain"]
        losses = [x for x in segs if status_of(x) == "loss"]
        bp_gain = sum(x["length"] for x in gains)
        bp_loss = sum(x["length"] for x in losses)
        fga = "" if not genome_size else f"{(bp_gain + bp_loss) / genome_size:.4f}"
        out.append({
            "sample": s, "n_segments": len(segs), "n_gain": len(gains), "n_loss": len(losses),
            "n_cn_loh": sum(1 for x in segs if is_cn_loh(x)),
            "bp_gain": bp_gain, "bp_loss": bp_loss, "bp_altered": bp_gain + bp_loss,
            "fraction_genome_altered": fga,
            "n_significant": sum(1 for x in segs if x["significant"] == "true"),
        })
    return out


def write_tsv(path: str, rows: list[dict], columns: list[str], comment: str = "") -> None:
    with open(path, "w", newline="") as fh:
        if comment:
            fh.write(comment)
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


MQC_HEADER = """# id: 'freec_cnv_summary'
# section_name: 'Control-FREEC CNV summary'
# description: 'Per-sample copy-number segments called by Control-FREEC. FGA = fraction of the genome (chrLenFile) in gain or loss segments.'
# plot_type: 'table'
"""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("cnv_files", nargs="+", help="*_CNVs or *_CNVs.p.value.txt files")
    ap.add_argument("--chr-len", help="chromosome length file, for the fraction of genome altered")
    ap.add_argument("--alpha", type=float, default=0.05, help="p-value threshold for 'significant' [0.05]")
    ap.add_argument("--outdir", default=".", help="output directory [.]")
    ap.add_argument("--name-from-dir", action="store_true",
                    help="use each file's parent directory as the sample name (layout <sample>/<file>)")
    args = ap.parse_args(argv)

    chosen: dict[str, str] = {}
    for p in args.cnv_files:
        s = os.path.basename(os.path.dirname(os.path.abspath(p))) if args.name_from_dir else sample_name(p)
        if s not in chosen or p.endswith(".p.value.txt"):
            chosen[s] = p
    samples = sorted(chosen)
    segments = []
    for s in samples:
        segments.extend(parse_cnv_file(chosen[s], s, args.alpha))
    segments.sort(key=lambda x: (x["sample"], x["chrom"], x["start"]))

    os.makedirs(args.outdir, exist_ok=True)
    summary = summarise(segments, samples, read_genome_size(args.chr_len))
    write_tsv(os.path.join(args.outdir, "cohort_cnv_segments.tsv"), segments, SEGMENT_COLUMNS)
    write_tsv(os.path.join(args.outdir, "cohort_cnv_summary.tsv"), summary, SUMMARY_COLUMNS)
    write_tsv(os.path.join(args.outdir, "freec_cnv_summary_mqc.tsv"), summary, SUMMARY_COLUMNS, MQC_HEADER)
    print(f"{len(samples)} samples, {len(segments)} segments", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
