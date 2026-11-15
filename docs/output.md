# Output

| Path | Content |
|---|---|
| `configs/<sample>.config.txt` | the FREEC configuration used for the sample |
| `samples/<sample>/<sample>_CNVs` | predicted CNVs: chromosome, start, end, copy number, gain/loss; with BAF also genotype (e.g. AAB) and uncertainty |
| `samples/<sample>/<sample>_CNVs.p.value.txt` | the same with Wilcoxon rank-sum and Kolmogorov-Smirnov p-values (`assess_significance.R`) |
| `samples/<sample>/<sample>_ratio.txt` | per-window ratio, median ratio and copy number |
| `samples/<sample>/<sample>_ratio.BedGraph` | ratio track for genome browsers |
| `samples/<sample>/<sample>_BAF.txt`, `_normal_BAF.txt` | B-allele frequencies of tumour (and normal) |
| `samples/<sample>/<sample>_info.txt` | FREEC run information, including estimated ploidy and contamination |
| `samples/<sample>/<sample>_sample.cpn`, `_control.cpn` | raw read counts per window |
| `samples/<sample>/<sample>.CNV.bed` | copy-number segments as BED (`freec2bed.pl`) |
| `samples/<sample>/<sample>.circos.txt` | Circos track (`freec2circos.pl`) |
| `samples/<sample>/plots/` | `makeGraph2.0.R` plots of ratio (linear and log2) and BAF |
| `cohort/cohort_cnv_segments.tsv` | all segments from all samples; `significant` is true when both p-values are below `--alpha` |
| `cohort/cohort_cnv_summary.tsv` | per sample: number of segments, gains, losses, copy-neutral LOH segments, bases gained and lost, fraction of genome altered (relative to `chrLenFile`), significant segments |
| `multiqc/multiqc_report.html` | the cohort summary as a table |
| `pipeline_info/` | execution report, timeline, trace, `software_versions.yml` |

FREEC writes chromosome names without the `chr` prefix in `_CNVs` files; the cohort tables keep them as written by FREEC.
