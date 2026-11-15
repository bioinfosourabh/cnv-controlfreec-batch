import configparser

import pytest

HEADER = "sample_id,tumor_bam,normal_bam,sex\n"


def parse(text):
    cp = configparser.ConfigParser(delimiters=("=",))
    cp.optionxform = str
    cp.read_string(text)
    return cp


def run(mfc, tmp_path, sheet_path, template_path, *extra):
    out = tmp_path / "configs"
    rc = mfc.main(["-s", str(sheet_path), "-t", str(template_path), "-o", str(out),
                   "--chr-len", "hg38.len", "--chr-files", "chroms", *extra])
    return rc, out


def test_tumour_normal_exome(mfc, tmp_path, sheet, template_path):
    p = sheet(HEADER + "S1,S1_SO.bam,S1_GL.bam,XY\n")
    rc, out = run(mfc, tmp_path, p, template_path, "--capture-bed", "targets.bed")
    assert rc == 0
    cfg = parse((out / "S1.config.txt").read_text())
    assert cfg["general"]["window"] == "0"
    assert cfg["general"]["sex"] == "XY"
    assert cfg["general"]["degree"] == "1"
    assert cfg["general"]["outputDir"] == "S1"
    assert cfg["sample"]["mateFile"] == "S1_SO.bam"
    assert cfg["control"]["mateFile"] == "S1_GL.bam"
    assert cfg["target"]["captureRegions"] == "targets.bed"
    assert "BAF" not in cfg


def test_tumour_only_drops_control_and_uses_gc_defaults(mfc, tmp_path, sheet, template_path):
    p = sheet(HEADER + "S2,S2_SO.bam,,XX\n")
    rc, out = run(mfc, tmp_path, p, template_path)
    assert rc == 0
    text = (out / "S2.config.txt").read_text()
    cfg = parse(text)
    assert "control" not in cfg
    assert "target" not in cfg
    assert cfg["general"]["window"] == "50000"
    assert cfg["general"]["degree"] == "3"
    assert cfg["general"]["intercept"] == "1"
    assert "{{" not in text


def test_empty_placeholders_drop_lines(mfc, tmp_path, sheet, template_path):
    p = sheet(HEADER + "S1,S1_SO.bam,S1_GL.bam,XY\n")
    rc, out = run(mfc, tmp_path, p, template_path)
    cfg = parse((out / "S1.config.txt").read_text())
    for key in ("sambamba", "bedtools", "gemMappabilityFile", "step"):
        assert key not in cfg["general"]


def test_baf_section_with_snp_file(mfc, tmp_path, sheet, template_path):
    p = sheet(HEADER + "S1,S1_SO.bam,S1_GL.bam,XY\n")
    rc, out = run(mfc, tmp_path, p, template_path, "--snp-file", "snps.vcf.gz", "--fasta", "ref.fa")
    cfg = parse((out / "S1.config.txt").read_text())
    assert cfg["BAF"]["SNPfile"] == "snps.vcf.gz"
    # FREEC expects a file of SNP positions here, not TRUE.
    assert cfg["BAF"]["makePileup"] == "snps.vcf.gz"
    assert cfg["BAF"]["fastaFile"] == "ref.fa"


def test_sheet_columns_override_per_sample(mfc, tmp_path, sheet, template_path):
    p = sheet("sample_id,tumor_bam,normal_bam,sex,contamination,ploidy\n"
              "S1,S1_SO.bam,S1_GL.bam,XY,0.4,3\nS2,S2_SO.bam,S2_GL.bam,XX,,\n")
    rc, out = run(mfc, tmp_path, p, template_path, "--set", "CONTAMINATION=0.1")
    c1 = parse((out / "S1.config.txt").read_text())
    c2 = parse((out / "S2.config.txt").read_text())
    assert c1["general"]["contamination"] == "0.4"
    assert c1["general"]["ploidy"] == "3"
    assert c2["general"]["contamination"] == "0.1"   # --set used when the cell is empty
    assert c2["general"]["ploidy"] == "2"            # built-in default


def test_output_dir_pattern(mfc, tmp_path, sheet, template_path):
    p = sheet(HEADER + "S1,S1_SO.bam,S1_GL.bam,XY\n")
    rc, out = run(mfc, tmp_path, p, template_path, "--output-dir", "results/{sample_id}")
    assert parse((out / "S1.config.txt").read_text())["general"]["outputDir"] == "results/S1"


def test_only_selected_samples(mfc, tmp_path, sheet, template_path):
    p = sheet(HEADER + "S1,S1_SO.bam,S1_GL.bam,XY\nS2,S2_SO.bam,,XX\n")
    rc, out = run(mfc, tmp_path, p, template_path, "--sample", "S2")
    assert rc == 0
    assert sorted(f.name for f in out.iterdir()) == ["S2.config.txt"]


@pytest.mark.parametrize("body, message", [
    ("S1,S1_SO.bam,S1_GL.bam,M\n", "sex must be XX or XY"),
    ("S1,S1_SO.bam,S1_GL.bam,XY\nS1,S2_SO.bam,,XX\n", "duplicate sample_id"),
    ("S 1,S1_SO.bam,S1_GL.bam,XY\n", "invalid sample_id"),
    (",S1_SO.bam,S1_GL.bam,XY\n", "invalid sample_id"),
    ("S1,,S1_GL.bam,XY\n", "tumor_bam is empty"),
    ("S1,S1_SO.fastq.gz,,XY\n", "must be a .bam or .cram"),
])
def test_invalid_rows(mfc, tmp_path, sheet, template_path, capsys, body, message):
    p = sheet(HEADER + body)
    rc, _ = run(mfc, tmp_path, p, template_path)
    assert rc == 1
    assert message in capsys.readouterr().err


def test_missing_column(mfc, tmp_path, sheet, template_path, capsys):
    p = sheet("sample_id,tumor_bam,sex\nS1,S1_SO.bam,XY\n")
    rc, _ = run(mfc, tmp_path, p, template_path)
    assert rc == 1
    assert "normal_bam" in capsys.readouterr().err


def test_unknown_placeholder_is_an_error(mfc, tmp_path, sheet, capsys):
    tpl = tmp_path / "tpl.txt"
    tpl.write_text("[general]\nchrLenFile = {{CHR_LEN_FILE}}\nfoo = {{TYPO_NAME}}\n[sample]\nmateFile = {{TUMOR_BAM}}\n")
    p = sheet(HEADER + "S1,S1_SO.bam,,XY\n")
    rc, _ = run(mfc, tmp_path, p, tpl)
    assert rc == 1
    assert "TYPO_NAME" in capsys.readouterr().err


def test_check_files(mfc, tmp_path, sheet, template_path, capsys):
    p = sheet(HEADER + "S1,/nonexistent/S1_SO.bam,,XY\n")
    rc, _ = run(mfc, tmp_path, p, template_path, "--check-files")
    assert rc == 1
    assert "not found" in capsys.readouterr().err
