#!/usr/bin/env python3
"""Simulate a small tumour/normal data set with known copy-number changes.

Reference: two contigs (chr1 240 kb, chr2 160 kb) whose GC content varies in
20 kb blocks, so that GC-content normalisation has something to fit.
Truth (tumour only, 100% purity in the simulation):
  chr1:80,001-140,000   gain  (3 copies)
  chr2:40,001-100,000   loss  (1 copy, with loss of heterozygosity)
Normal and tumours share 300 heterozygous SNPs (listed in snps.vcf.gz).

Writes genome.fasta(.fai), snps.vcf.gz(.tbi), targets.bed (pseudo-exons
every 2 kb) and coordinate-sorted, indexed BAMs:
  S1_SO.bam / S1_GL.bam   tumour / normal of patient S1
  S2_SO.bam               tumour without a normal

Requirements: wgsim, bwa, samtools, bgzip, tabix on PATH.
    python3 tests/make_test_data.py tests/data
"""
import os
import random
import subprocess
import sys

SEED = 7
CONTIGS = {"chr1": 240_000, "chr2": 160_000}
GAIN = ("chr1", 80_000, 140_000)
LOSS = ("chr2", 40_000, 100_000)
COVERAGE = 12
READ_LEN = 100


def run(cmd, **kw):
    subprocess.run(cmd, check=True, **kw)


def make_reference(rng):
    seqs = {}
    for name, length in CONTIGS.items():
        blocks = []
        for start in range(0, length, 20_000):
            gc = rng.uniform(0.34, 0.56)
            at = (1 - gc) / 2
            blocks.append("".join(rng.choices("ACGT", weights=[at, gc / 2, gc / 2, at], k=min(20_000, length - start))))
        seqs[name] = "".join(blocks)
    return seqs


def write_fasta(path, seqs):
    with open(path, "w") as fh:
        for n, s in seqs.items():
            fh.write(f">{n}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")


def mutate(seqs, snps):
    out = {k: list(v) for k, v in seqs.items()}
    for c, p, r, a in snps:
        out[c][p - 1] = a
    return {k: "".join(v) for k, v in out.items()}


def simulate(fa, pairs, seed, tag, r1, r2):
    t1, t2 = r1 + ".tmp", r2 + ".tmp"
    run(["wgsim", "-N", str(pairs), "-1", str(READ_LEN), "-2", str(READ_LEN), "-e", "0.001",
         "-r", "0", "-R", "0", "-d", "300", "-s", "30", "-S", str(seed), fa, t1, t2],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for t, f in ((t1, r1), (t2, r2)):
        with open(t) as src, open(f, "a") as dst:
            for i, line in enumerate(src):
                dst.write("@" + tag + "_" + line[1:] if i % 4 == 0 else line)
        os.remove(t)


def region_fasta(seqs, region, path):
    c, s, e = region
    write_fasta(path, {f"{c}_{s}": seqs[c][s:e]})


def main(outdir):
    rng = random.Random(SEED)
    os.makedirs(outdir, exist_ok=True)
    ref = make_reference(rng)
    fa = os.path.join(outdir, "genome.fasta")
    write_fasta(fa, ref)
    run(["samtools", "faidx", fa])
    run(["bwa", "index", fa], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    snps = []
    for _ in range(300):
        c = rng.choice(list(CONTIGS))
        p = rng.randint(500, CONTIGS[c] - 500)
        r = ref[c][p - 1]
        snps.append((c, p, r, rng.choice([b for b in "ACGT" if b != r])))
    snps = sorted(set(snps), key=lambda v: (list(CONTIGS).index(v[0]), v[1]))
    vcf = os.path.join(outdir, "snps.vcf")
    with open(vcf, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        for c, length in CONTIGS.items():
            fh.write(f"##contig=<ID={c},length={length}>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for i, (c, p, r, a) in enumerate(snps):
            fh.write(f"{c}\t{p}\tsnp{i}\t{r}\t{a}\t.\tPASS\t.\n")
    run(["bgzip", "-f", vcf])
    run(["tabix", "-f", "-p", "vcf", vcf + ".gz"])

    with open(os.path.join(outdir, "targets.bed"), "w") as fh:
        for c, length in CONTIGS.items():
            for s in range(1_000, length - 1_000, 2_000):
                fh.write(f"{c}\t{s}\t{s + 600}\n")

    hap_a = os.path.join(outdir, "hapA.tmp.fa")
    hap_b = os.path.join(outdir, "hapB.tmp.fa")
    write_fasta(hap_a, ref)
    hb = mutate(ref, snps)
    write_fasta(hap_b, hb)
    # Tumour haplotype B lacks the LOSS region; the GAIN region is added once more from hap A.
    c, s, e = LOSS
    hb_del = dict(hb)
    hb_del[c] = hb[c][:s] + hb[c][e:]
    hap_b_del = os.path.join(outdir, "hapBdel.tmp.fa")
    write_fasta(hap_b_del, hb_del)
    gain_fa = os.path.join(outdir, "gain.tmp.fa")
    region_fasta(ref, GAIN, gain_fa)

    genome = sum(CONTIGS.values())
    per_hap = genome * COVERAGE // (2 * 2 * READ_LEN)
    gain_pairs = (GAIN[2] - GAIN[1]) * COVERAGE // (2 * 2 * READ_LEN)
    del_len = LOSS[2] - LOSS[1]
    plans = {
        "S1_GL": [(hap_a, per_hap), (hap_b, per_hap)],
        "S1_SO": [(hap_a, per_hap), (hap_b_del, per_hap * (genome - del_len) // genome), (gain_fa, gain_pairs)],
        "S2_SO": [(hap_a, per_hap), (hap_b_del, per_hap * (genome - del_len) // genome), (gain_fa, gain_pairs)],
    }
    for k, (sample, parts) in enumerate(plans.items()):
        r1 = os.path.join(outdir, f"{sample}_R1.fq")
        r2 = os.path.join(outdir, f"{sample}_R2.fq")
        for f in (r1, r2):
            if os.path.exists(f):
                os.remove(f)
        for j, (src, n) in enumerate(parts):
            simulate(src, n, SEED * 100 + k * 10 + j, f"{sample}h{j}", r1, r2)
        bam = os.path.join(outdir, f"{sample}.bam")
        rg = f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA"
        bwa = subprocess.Popen(["bwa", "mem", "-t", "2", "-R", rg, fa, r1, r2],
                               stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        run(["samtools", "sort", "-o", bam, "-"], stdin=bwa.stdout, stderr=subprocess.DEVNULL)
        bwa.wait()
        run(["samtools", "index", bam])
        os.remove(r1)
        os.remove(r2)

    for f in os.listdir(outdir):
        if ".tmp." in f or f.startswith("genome.fasta.") and f.split(".")[-1] in {"amb", "ann", "bwt", "pac", "sa"}:
            os.remove(os.path.join(outdir, f))
    with open(os.path.join(outdir, "truth_cnv.tsv"), "w") as fh:
        fh.write("chrom\tstart\tend\tcopy_number\tsamples\n")
        fh.write(f"{GAIN[0]}\t{GAIN[1]}\t{GAIN[2]}\t3\tS1_SO,S2_SO\n")
        fh.write(f"{LOSS[0]}\t{LOSS[1]}\t{LOSS[2]}\t1\tS1_SO,S2_SO\n")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "tests/data")
