"""Tests for organize-mags naming helpers (mamisa.commands.filter_checkm2)."""

from pathlib import Path

import pytest

from mamisa.commands.filter_checkm2 import (
    sanitize_token,
    extract_sample,
    taxonomy_label,
    build_output_name,
    DEFAULT_SAMPLE_REGEX,
)


@pytest.mark.parametrize("raw,clean", [
    ("Escherichia coli", "Escherichia_coli"),
    ("Firmicutes_A", "Firmicutes_A"),
    ("weird/name;x", "weird_name_x"),
    ("", "NA"),
])
def test_sanitize_token(raw, clean):
    assert sanitize_token(raw) == clean


@pytest.mark.parametrize("name,regex,sample", [
    ("SAMPLE1.bin.3", DEFAULT_SAMPLE_REGEX, "SAMPLE1"),
    ("SAMPLE1_bin_3", DEFAULT_SAMPLE_REGEX, "SAMPLE1"),
    ("S12_bin_4", r"^(S[0-9]+)_", "S12"),
])
def test_extract_sample(name, regex, sample):
    assert extract_sample(name, regex) == sample


def test_extract_sample_bad_regex_raises():
    with pytest.raises(ValueError):
        extract_sample("x", "(")


def test_taxonomy_label_uses_genus():
    cls = "d__Bacteria;p__Firmicutes;g__Lactobacillus;s__"
    assert taxonomy_label(cls, "genus") == "Lactobacillus"


def test_taxonomy_label_falls_back_up_ranks():
    # genus empty -> fall back up to phylum
    cls = "d__Bacteria;p__Bacteroidota;c__;o__;f__;g__;s__"
    assert taxonomy_label(cls, "genus") == "Bacteroidota"


def test_taxonomy_label_unclassified():
    assert taxonomy_label("", "genus") == "Unclassified"


def test_build_output_name_with_taxonomy():
    tax = {"SAMPLE1.bin.1": {"classification":
                             "d__Bacteria;p__Firmicutes;g__Lactobacillus;s__"}}
    out = build_output_name(Path("SAMPLE1.bin.1.fa"), "SAMPLE1.bin.1",
                            rename=True, tax=tax,
                            sample_regex=DEFAULT_SAMPLE_REGEX, tax_level="genus")
    assert out == "SAMPLE1__Lactobacillus__SAMPLE1.bin.1.fa"


def test_build_output_name_no_tax_hit():
    out = build_output_name(Path("S2.bin.1.fa"), "S2.bin.1",
                            rename=True, tax={},
                            sample_regex=DEFAULT_SAMPLE_REGEX, tax_level="genus")
    assert out == "S2__NoTax__S2.bin.1.fa"


def test_build_output_name_rename_disabled():
    out = build_output_name(Path("S2.bin.1.fa"), "S2.bin.1",
                            rename=False, tax={},
                            sample_regex=DEFAULT_SAMPLE_REGEX, tax_level="genus")
    assert out == "S2.bin.1.fa"
