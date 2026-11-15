import csv

CNVS_WITH_BAF = "1\t81000\t135600\t3\tgain\tAAB\t22.0\n2\t35000\t99600\t1\tloss\tA\t-1\n3\t1000\t50000\t2\tnormal\tAA\t1.0\n"
CNVS_PLAIN = "1\t80000\t140000\t3\tgain\n"
PVALUES = ("chr\tstart\tend\tcopy number\tstatus\tWilcoxonRankSumTestPvalue\tKolmogorovSmirnovPvalue\n"
           "1\t80000\t140000\t3\tgain\t1e-10\t2e-8\n2\t10\t20\t1\tloss\t0.3\t0.01\n")


def read(path):
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def test_summary_counts_and_fga(summary, tmp_path):
    (tmp_path / "A.bam_CNVs").write_text(CNVS_WITH_BAF)
    (tmp_path / "B_CNVs").write_text(CNVS_PLAIN)
    chr_len = tmp_path / "genome.len"
    chr_len.write_text("1\t200000\n2\t100000\n3\t100000\n")
    rc = summary.main([str(tmp_path / "A.bam_CNVs"), str(tmp_path / "B_CNVs"),
                       "--chr-len", str(chr_len), "--outdir", str(tmp_path / "out")])
    assert rc == 0
    rows = {r["sample"]: r for r in read(tmp_path / "out" / "cohort_cnv_summary.tsv")}
    assert set(rows) == {"A", "B"}
    assert rows["A"]["n_gain"] == "1" and rows["A"]["n_loss"] == "1"
    assert rows["A"]["n_cn_loh"] == "1"
    assert rows["A"]["bp_altered"] == str((135600 - 81000) + (99600 - 35000))
    assert rows["A"]["fraction_genome_altered"] == f"{((135600 - 81000) + (99600 - 35000)) / 400000:.4f}"
    segs = read(tmp_path / "out" / "cohort_cnv_segments.tsv")
    assert len(segs) == 4
    assert segs[0]["genotype"] == "AAB"


def test_pvalue_file_preferred_and_significance(summary, tmp_path):
    d = tmp_path / "S1"
    d.mkdir()
    (d / "S1_SO.bam_CNVs").write_text(CNVS_PLAIN)
    (d / "S1_SO.bam_CNVs.p.value.txt").write_text(PVALUES)
    rc = summary.main([str(d / "S1_SO.bam_CNVs"), str(d / "S1_SO.bam_CNVs.p.value.txt"),
                       "--outdir", str(tmp_path / "out"), "--name-from-dir"])
    assert rc == 0
    segs = read(tmp_path / "out" / "cohort_cnv_segments.tsv")
    assert [s["sample"] for s in segs] == ["S1", "S1"]
    assert [s["significant"] for s in segs] == ["true", "false"]
    rows = read(tmp_path / "out" / "cohort_cnv_summary.tsv")
    assert rows[0]["n_significant"] == "1"
    assert rows[0]["fraction_genome_altered"] == ""
    mqc = (tmp_path / "out" / "freec_cnv_summary_mqc.tsv").read_text()
    assert mqc.startswith("# id: 'freec_cnv_summary'")


def test_sample_name_parsing(summary):
    assert summary.sample_name("/x/S1_SO.bam_CNVs") == "S1_SO"
    assert summary.sample_name("S1_CNVs.p.value.txt") == "S1"
    assert summary.sample_name("P7.bam_CNVs.p.value.txt") == "P7"


def test_zero_length_records_skipped(summary, tmp_path):
    (tmp_path / "Z_CNVs").write_text("1\t240000\t240000\t0\tloss\n1\t10\t20\t3\tgain\n")
    summary.main([str(tmp_path / "Z_CNVs"), "--outdir", str(tmp_path / "out")])
    segs = read(tmp_path / "out" / "cohort_cnv_segments.tsv")
    assert len(segs) == 1 and segs[0]["status"] == "gain"
