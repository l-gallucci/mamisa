"""Tests for mamisa.utils.chimera."""

import math
import pytest

from mamisa.utils.chimera import (
    compute_gc,
    windowed_gc,
    gc_statistics,
    extract_taxonomy_level,
    parse_gtdbtk_summary,
    parse_gunc_output,
    has_taxonomy_warning,
    assess_chimera_risk,
)


def test_compute_gc():
    assert compute_gc("GGCC") == 1.0
    assert compute_gc("ATAT") == 0.0
    assert compute_gc("ACGT") == 0.5
    assert math.isnan(compute_gc("NNNN"))


def test_windowed_gc_short_sequence():
    # shorter than window -> single whole-sequence value
    out = windowed_gc("ACGTACGT", window=5000, step=2500)
    assert len(out) == 1 and out[0][0] == 0


def test_windowed_gc_sliding():
    seq = "GC" * 5000          # 10 kbp, all GC
    out = windowed_gc(seq, window=1000, step=1000)
    assert len(out) >= 9
    assert all(abs(gc - 1.0) < 1e-9 for _, gc in out)


def test_gc_statistics():
    s = gc_statistics([0.4, 0.6])
    assert s["mean"] == pytest.approx(0.5)
    assert s["delta"] == pytest.approx(0.2)
    assert s["n"] == 2
    empty = gc_statistics([float("nan")])
    assert empty["n"] == 0


@pytest.mark.parametrize("level,expected", [
    ("genus", "Lactobacillus"),
    ("phylum", "Firmicutes"),
    ("species", "Unclassified"),   # empty s__
])
def test_extract_taxonomy_level(level, expected):
    cls = "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lacto;f__Lacto;g__Lactobacillus;s__"
    assert extract_taxonomy_level(cls, level) == expected


def test_parse_gtdbtk_summary(tmp_path):
    f = tmp_path / "gtdbtk.bac120.summary.tsv"
    f.write_text("user_genome\tclassification\tmsa_percent\twarnings\n"
                 "binA\td__Bacteria;p__Firmicutes\t95.0\t\n")
    out = parse_gtdbtk_summary(f)
    assert out["binA"]["classification"].startswith("d__Bacteria")
    assert out["binA"]["msa_percent"] == 95.0


def test_parse_gunc_output(tmp_path):
    f = tmp_path / "GUNC.progenomes_2.1.maxCSS_level.tsv"
    f.write_text(
        "genome\tclade_separation_score\tcontamination_portion\t"
        "n_effective_surplus_clades\tpass.GUNC\n"
        "binA\t0.82\t0.4\t2.1\tFalse\n"
        "binB\t0.10\t0.02\t0.3\tTrue\n")
    out = parse_gunc_output(tmp_path)
    assert out["binA"]["pass_gunc"] is False
    assert out["binA"]["css"] == 0.82
    assert out["binB"]["pass_gunc"] is True


def test_has_taxonomy_warning():
    assert has_taxonomy_warning({"warnings": "something"}) is True
    assert has_taxonomy_warning({"msa_percent": 5.0}) is True
    assert has_taxonomy_warning({"msa_percent": 95.0, "warnings": ""}) is False


def _risk(**kw):
    base = dict(gc_delta=0.0, gc_cv=0.0, n_contigs=1,
                contamination=None, windowed_gc_delta=None)
    base.update(kw)
    return assess_chimera_risk(**base)


def test_assess_chimera_risk_clean():
    risk, reasons = _risk()
    assert risk == "Clean" and reasons == []


def test_assess_chimera_risk_gc_high():
    risk, reasons = _risk(gc_delta=0.12, n_contigs=3)
    assert risk == "Medium"           # +3 for high inter-contig GC
    assert any("GC heterogeneity" in r for r in reasons)


def test_assess_chimera_risk_gunc_fail_adds_three():
    risk, reasons = _risk(gunc_fail=True, gunc_css=0.82)
    assert risk == "Medium"           # +3
    assert any("GUNC fail" in r for r in reasons)


def test_assess_chimera_risk_gunc_css_elevated():
    risk, reasons = _risk(gunc_fail=False, gunc_css=0.5)
    assert any("clade separation" in r.lower() for r in reasons)


def test_assess_chimera_risk_high_combo():
    # high GC (+3) + GUNC fail (+3) = 6 -> High
    risk, _ = _risk(gc_delta=0.2, n_contigs=2, gunc_fail=True, gunc_css=0.9)
    assert risk == "High"
