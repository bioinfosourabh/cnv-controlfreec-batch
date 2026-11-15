# Usage

## Inputs

| Parameter | Description |
|---|---|
| `--input` | sample sheet (`sample_id,tumor_bam,normal_bam,sex`, extra columns override placeholders). Relative BAM paths are resolved against the sheet's directory. BAM indexes (`.bam.bai` or `.bai`) are used if present and built otherwise. |
| `--fasta` | the reference the BAMs were aligned to |
| `--chr_len`, `--chr_files` | optional pre-built `chrLenFile` and `chrFiles` directory; otherwise built from `--fasta`, keeping contigs that match `--contig_regex` |
| `--capture_bed` | exome or panel targets; switches to targeted defaults (`window = 0`) |
| `--snp_vcf` | SNP positions (bgzipped VCF with `.tbi`) for BAF and genotype calls, e.g. common dbSNP or gnomAD SNPs restricted to the targets |
| `--mappability` | optional GEM mappability file |

## Changing FREEC settings

- For all samples: `--freec_set WINDOW=50000,STEP=10000,CONTAMINATION=0.3`
- For one sample: add a column to the sample sheet, e.g. `contamination`
- Permanently: edit `assets/config_template.txt`, or copy it and pass `--template my_template.txt`

The placeholders and their defaults are listed in the README. Any FREEC option can be added to the template as a fixed value or as a new placeholder.

## Resources

FREEC runs with `--freec_threads` threads (default 4) and 16 GB of memory; other steps are small. Limits are set in `conf/base.config` (`process.resourceLimits`).

## Running the Bash version

`bin/run_controlfreec_batch.sh -s samples.csv -c freec_batch.conf -o results [-j N] [-t threads] [-n] [-f]`

The config file (see `conf/freec_batch.conf.example`) holds the FREEC binary, the helper-script directory and the reference paths. The script checks that the tools and files exist before starting, writes one log per sample and step under `results/logs/`, and exits non-zero if any sample failed.
