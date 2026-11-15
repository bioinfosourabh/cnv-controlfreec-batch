#!/usr/bin/env python3
"""Generate one Control-FREEC configuration file per sample from a sample sheet.

The sample sheet is a CSV file with at least the columns

    sample_id,tumor_bam,normal_bam,sex

`normal_bam` may be empty (tumour-only analysis). Any further column, for
example `ploidy`, `contamination` or `capture_bed`, overrides the template
placeholder of the same name (upper-cased) for that sample only.

Placeholders in the template are written as {{NAME}}. Values are resolved in
this order: sample sheet column, --set NAME=VALUE, command-line options for
reference files, built-in mode-dependent defaults. A line whose placeholders
all resolve to empty is dropped; a section whose `#@requires NAME` marker
resolves to empty is dropped with all its lines.

Examples
--------
Exome, tumour-normal, configs written to configs/<sample>.config.txt and
FREEC output directories results/<sample>:

    make_freec_configs.py -s samples.csv -t assets/config_template.txt \\
        -o configs --chr-len hg38.len --chr-files hg38_chroms/ \\
        --capture-bed exome_targets.bed --output-dir 'results/{sample_id}'

Whole genome (no capture BED) with BAF from dbSNP positions:

    make_freec_configs.py -s samples.csv -t assets/config_template.txt \\
        -o configs --chr-len hg38.len --chr-files hg38_chroms/ \\
        --snp-file dbsnp_common.vcf.gz --fasta hg38.fa --set WINDOW=50000
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from dataclasses import dataclass, field

REQUIRED_COLUMNS = ("sample_id", "tumor_bam", "normal_bam", "sex")
PLACEHOLDER = re.compile(r"\{\{([A-Z0-9_]+)\}\}")
REQUIRES = re.compile(r"^\s*#@requires\s+([A-Z0-9_]+)\s*$")
SAMPLE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
VALID_SEX = {"XX", "XY"}

# Placeholders that may legitimately be empty and have no default.
OPTIONAL_NAMES = {
    "CHR_FILES", "MAPPABILITY", "SAMTOOLS", "BEDTOOLS", "SAMBAMBA",
    "SAMBAMBA_THREADS", "STEP", "SNP_FILE",
    "SNP_POSITIONS", "FASTA", "CAPTURE_BED", "NORMAL_BAM",
}


class SheetError(ValueError):
    """Problem with the sample sheet or template."""


@dataclass
class Sample:
    sample_id: str
    values: dict = field(default_factory=dict)

    @property
    def has_normal(self) -> bool:
        return bool(self.values.get("NORMAL_BAM"))


def read_sample_sheet(path: str, check_files: bool = False) -> list[Sample]:
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        header = [h.strip() for h in (reader.fieldnames or [])]
        rows = [r for r in reader if any((v or "").strip() for v in r.values())]
    missing = [c for c in REQUIRED_COLUMNS if c not in header]
    if missing:
        raise SheetError(f"{path}: missing column(s): {', '.join(missing)}")
    if not rows:
        raise SheetError(f"{path}: no samples listed")

    samples, seen = [], set()
    for n, row in enumerate(rows, start=2):
        row = {k.strip(): (v or "").strip() for k, v in row.items() if k}
        sid = row["sample_id"]
        if not SAMPLE_ID.match(sid):
            raise SheetError(f"{path} line {n}: invalid sample_id '{sid}' (letters, digits, '.', '_', '-')")
        if sid in seen:
            raise SheetError(f"{path} line {n}: duplicate sample_id '{sid}'")
        seen.add(sid)
        sex = row["sex"].upper()
        if sex not in VALID_SEX:
            raise SheetError(f"{path} line {n}: sex must be XX or XY, got '{row['sex']}'")
        row["sex"] = sex
        if not row["tumor_bam"]:
            raise SheetError(f"{path} line {n}: tumor_bam is empty")
        for col in ("tumor_bam", "normal_bam"):
            p = row[col]
            if p and not p.endswith((".bam", ".cram")):
                raise SheetError(f"{path} line {n}: {col} must be a .bam or .cram file: {p}")
            if p and check_files and not os.path.exists(p):
                raise SheetError(f"{path} line {n}: {col} not found: {p}")
        values = {k.upper(): v for k, v in row.items()}
        values["SAMPLE_ID"] = sid
        samples.append(Sample(sid, values))
    return samples


def mode_defaults(sample: Sample, globals_: dict) -> dict:
    """Defaults that depend on whether the data are targeted and paired."""
    targeted = bool(sample.values.get("CAPTURE_BED") or globals_.get("CAPTURE_BED"))
    paired = sample.has_normal
    return {
        "PLOIDY": "2",
        "THREADS": "4",
        "CONTAMINATION": "0.25",
        # Exome: one window per capture region; genome: 50 kb windows.
        "WINDOW": "0" if targeted else "50000",
        "NOISY_DATA": "TRUE" if targeted else "FALSE",
        "READ_COUNT_THRESHOLD": "50" if targeted else "10",
        # Control-based normalisation uses a linear fit through the origin;
        # GC-content normalisation (tumour-only) a cubic with intercept.
        "DEGREE": "1" if paired else "3",
        "INTERCEPT": "0" if paired else "1",
        "OUTPUT_DIR": "{sample_id}",
    }


def resolve(name: str, sample: Sample, globals_: dict, defaults: dict) -> str | None:
    for source in (sample.values, globals_, defaults):
        if name in source and source[name] is not None and source[name] != "":
            return str(source[name])
    if name in sample.values or name in globals_ or name in OPTIONAL_NAMES:
        return ""
    return None


def render(template: str, sample: Sample, globals_: dict) -> str:
    defaults = mode_defaults(sample, globals_)
    out_sections: list[list[str]] = [[]]
    unknown: set[str] = set()

    def substitute(line: str) -> str | None:
        names = PLACEHOLDER.findall(line)
        if not names:
            return line
        values = {}
        for nm in names:
            v = resolve(nm, sample, globals_, defaults)
            if v is None:
                unknown.add(nm)
                v = ""
            values[nm] = v.replace("{sample_id}", sample.sample_id)
        if all(values[nm] == "" for nm in names):
            return None
        return PLACEHOLDER.sub(lambda m: values[m.group(1)], line)

    current_required_ok = True
    for raw in template.splitlines():
        line = raw.rstrip()
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            out_sections.append([stripped])
            current_required_ok = True
            continue
        m = REQUIRES.match(line)
        if m:
            v = resolve(m.group(1), sample, globals_, defaults)
            if v is None:
                unknown.add(m.group(1))
            current_required_ok = bool(v)
            if not current_required_ok:
                out_sections[-1] = []
            continue
        if not current_required_ok:
            continue
        if not stripped or stripped.startswith("#"):
            continue
        new = substitute(line)
        if new is not None:
            out_sections[-1].append(new)

    if unknown:
        raise SheetError(
            f"template placeholder(s) with no value for sample {sample.sample_id}: "
            f"{', '.join(sorted(unknown))} (add a sample-sheet column or --set NAME=VALUE)")

    blocks = [s for s in out_sections if len(s) > 1 or (s and not s[0].startswith("["))]
    return "\n\n".join("\n".join(b) for b in blocks) + "\n"


def parse_set(items: list[str]) -> dict:
    out = {}
    for item in items or []:
        if "=" not in item:
            raise SheetError(f"--set expects NAME=VALUE, got '{item}'")
        k, v = item.split("=", 1)
        k = k.strip().upper()
        if not re.fullmatch(r"[A-Z0-9_]+", k):
            raise SheetError(f"--set: invalid placeholder name '{k}'")
        out[k] = v.strip()
    return out


def build_globals(args) -> dict:
    g = {
        "CHR_LEN_FILE": args.chr_len,
        "CHR_FILES": args.chr_files,
        "CAPTURE_BED": args.capture_bed,
        "SNP_FILE": args.snp_file,
        "SNP_POSITIONS": args.snp_positions or args.snp_file,
        "FASTA": args.fasta,
        "MAPPABILITY": args.mappability,
        "THREADS": str(args.threads) if args.threads else None,
        "OUTPUT_DIR": args.output_dir,
    }
    g = {k: v for k, v in g.items() if v}
    g.update(parse_set(args.set))
    return g


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog=__doc__[__doc__.index("Examples"):])
    p.add_argument("-s", "--samplesheet", required=True, help="CSV: sample_id,tumor_bam,normal_bam,sex[,...]")
    p.add_argument("-t", "--template", required=True, help="Control-FREEC config template")
    p.add_argument("-o", "--outdir", default="configs", help="directory for the generated configs [configs]")
    p.add_argument("--sample", action="append", help="only write configs for these sample IDs (repeatable)")
    p.add_argument("--chr-len", help="chromosome length file (chrLenFile)")
    p.add_argument("--chr-files", help="directory of per-chromosome FASTA files (chrFiles)")
    p.add_argument("--capture-bed", help="capture regions BED for exome/panel data; omit for WGS")
    p.add_argument("--snp-file", help="SNP VCF for BAF (SNPfile); enables the [BAF] section")
    p.add_argument("--snp-positions", help="SNP positions for FREEC's own pileup (makePileup); defaults to --snp-file")
    p.add_argument("--fasta", help="reference FASTA (fastaFile, needed by makePileup)")
    p.add_argument("--mappability", help="GEM mappability file (gemMappabilityFile)")
    p.add_argument("--threads", type=int, help="maxThreads [4]")
    p.add_argument("--output-dir", help="FREEC outputDir; '{sample_id}' is replaced [{sample_id}]")
    p.add_argument("--set", action="append", metavar="NAME=VALUE", help="set any placeholder (repeatable)")
    p.add_argument("--check-files", action="store_true", help="fail if BAM files in the sheet do not exist")
    args = p.parse_args(argv)

    try:
        samples = read_sample_sheet(args.samplesheet, check_files=args.check_files)
        if args.sample:
            wanted = set(args.sample)
            missing = wanted - {s.sample_id for s in samples}
            if missing:
                raise SheetError(f"sample(s) not in sheet: {', '.join(sorted(missing))}")
            samples = [s for s in samples if s.sample_id in wanted]
        with open(args.template) as fh:
            template = fh.read()
        globals_ = build_globals(args)
        os.makedirs(args.outdir, exist_ok=True)
        for s in samples:
            text = render(template, s, globals_)
            path = os.path.join(args.outdir, f"{s.sample_id}.config.txt")
            with open(path, "w") as fh:
                fh.write(text)
            mode = "tumour-normal" if s.has_normal else "tumour-only"
            print(f"{s.sample_id}\t{mode}\t{path}")
    except (SheetError, OSError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
