"""Tests for the Kaiju parsers in mamisa.utils.read_taxonomy."""

from mamisa.utils.read_taxonomy import (
    parse_kaiju_output,
    parse_names_dmp,
    taxon_label,
)


KAIJU_OUT = "\t".join  # tab joiner helper


def test_parse_kaiju_output_basic(tmp_path):
    # Kaiju raw output: C/U, read_id, taxid, [extra cols ignored]
    content = (
        "C\tread1\t562\t50\t562\tACCESSION\tACGT\n"   # verbose (-v) row
        "C\tread2\t1280\n"                              # minimal 3-col row
        "U\tread3\t0\n"                                 # unclassified
        "U\tread4\t\n"                                  # unclassified, empty taxid
        "\n"                                            # blank line skipped
    )
    f = tmp_path / "kaiju.out"
    f.write_text(content)

    result = parse_kaiju_output(f)
    assert result == {"read1": 562, "read2": 1280, "read3": 0, "read4": 0}


def test_parse_kaiju_output_malformed_skipped(tmp_path):
    content = (
        "C\tread1\t562\n"
        "garbage_line_no_tabs\n"      # <3 fields -> skipped
        "C\tread2\tNOTANUMBER\n"      # non-integer taxid (classified) -> skipped
        "C\tread3\t99\n"
    )
    f = tmp_path / "kaiju.out"
    f.write_text(content)

    result = parse_kaiju_output(f)
    assert result == {"read1": 562, "read3": 99}


def test_parse_names_dmp_scientific_only(tmp_path):
    content = (
        "1\t|\troot\t|\t\t|\tscientific name\t|\n"
        "562\t|\tEscherichia coli\t|\t\t|\tscientific name\t|\n"
        "562\t|\tE. coli\t|\t\t|\tcommon name\t|\n"      # non-scientific -> ignored
        "1280\t|\tStaphylococcus aureus\t|\t\t|\tscientific name\t|\n"
    )
    f = tmp_path / "names.dmp"
    f.write_text(content)

    info = parse_names_dmp(f)
    assert info[562] == ("-", "Escherichia coli")
    assert info[1280] == ("-", "Staphylococcus aureus")
    # common-name row must not overwrite the scientific name
    assert info[562][1] == "Escherichia coli"
    # defaults present
    assert info[0] == ("-", "unclassified")


def test_taxon_label_with_names_dmp(tmp_path):
    f = tmp_path / "names.dmp"
    f.write_text("562\t|\tEscherichia coli\t|\t\t|\tscientific name\t|\n")
    info = parse_names_dmp(f)

    assert taxon_label(562, info) == "Escherichia coli [-:562]"
    assert taxon_label(0, info) == "unclassified"
    assert taxon_label(99999, info) == "taxid:99999"   # unknown -> fallback
