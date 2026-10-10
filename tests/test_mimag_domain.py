"""Tests for run-mimag-rna --split-by-domain helpers."""

from pathlib import Path

from mamisa.commands.run_mimag_rna import (
    load_domains_from_gtdbtk,
    resolve_kingdom,
    _strip_ext,
)


def _write_gtdbtk(root: Path):
    hq = root / "HQ"
    hq.mkdir(parents=True)
    (hq / "gtdbtk.bac120.summary.tsv").write_text(
        "user_genome\tclassification\n"
        "SO301_FLA_LB110_001\td__Bacteria;p__Pseudomonadota\n"
        "binB\td__Bacteria\n"
    )
    (hq / "gtdbtk.ar53.summary.tsv").write_text(
        "user_genome\tclassification\n"
        "arcBin_5\td__Archaea;p__Thermoproteota\n"
    )


def test_load_domains_by_filename(tmp_path):
    _write_gtdbtk(tmp_path)
    m = load_domains_from_gtdbtk(tmp_path)
    assert m["SO301_FLA_LB110_001"] == "bac"
    assert m["binB"] == "bac"
    assert m["arcBin_5"] == "arc"


def test_load_domains_classification_fallback(tmp_path):
    # filename matches gtdbtk.*.summary.tsv but has no bac120/ar53 token
    # -> domain taken from the classification column
    (tmp_path / "gtdbtk.other.summary.tsv").write_text(
        "user_genome\tclassification\n"
        "g1\td__Archaea;p__x\n"
        "g2\td__Bacteria;p__y\n"
    )
    m = load_domains_from_gtdbtk(tmp_path)
    assert m["g1"] == "arc"
    assert m["g2"] == "bac"


def test_resolve_kingdom_matches_with_extension(tmp_path):
    _write_gtdbtk(tmp_path)
    m = load_domains_from_gtdbtk(tmp_path)
    assert resolve_kingdom(Path("SO301_FLA_LB110_001.fa"), m, "bac") == "bac"
    assert resolve_kingdom(Path("arcBin_5.fasta"), m, "bac") == "arc"


def test_resolve_kingdom_fallback_when_missing(tmp_path):
    _write_gtdbtk(tmp_path)
    m = load_domains_from_gtdbtk(tmp_path)
    assert resolve_kingdom(Path("not_classified.fa"), m, "arc") == "arc"


def test_load_domains_empty_dir(tmp_path):
    assert load_domains_from_gtdbtk(tmp_path) == {}


def test_strip_ext():
    assert _strip_ext("bin.fa") == "bin"
    assert _strip_ext("bin.fasta.gz") == "bin"
    assert _strip_ext("bin_no_ext") == "bin_no_ext"
