"""Tests for mamisa.utils.fasta."""

import gzip
import pytest

from mamisa.utils.fasta import extract_contig_id, read_fasta_streaming


@pytest.mark.parametrize("name,expected", [
    ("k141_456", "456"),            # MEGAHIT: id after k<kmer>_, not the kmer size
    ("k99_12", "12"),
    ("NODE_1_length_5000_cov_3.5", "1"),   # SPAdes: id after NODE_
    ("NODE_456_length_200", "456"),
    ("contig_123", "123"),          # trailing integer
    ("scaffold_789", "789"),
    ("stdin.part_contig_7420", "7420"),
    ("weirdname", "weirdname"),      # no integer -> stripped original
    ("  spaced  ", "spaced"),
])
def test_extract_contig_id(name, expected):
    assert extract_contig_id(name) == expected


def test_megahit_not_confused_by_kmer():
    # regression: the old implementation returned the k-mer size (141), not 456
    assert extract_contig_id("k141_456") != "141"


def _write(path, text):
    path.write_text(text)
    return path


def test_read_fasta_streaming_plain(tmp_path):
    fa = _write(tmp_path / "a.fa", ">c1 desc\nACGT\nAAAA\n>c2\nGGCC\n")
    records = list(read_fasta_streaming(fa))
    assert records == [("c1", "ACGTAAAA"), ("c2", "GGCC")]  # only first word of header


def test_read_fasta_streaming_gzip(tmp_path):
    fa = tmp_path / "b.fa.gz"
    with gzip.open(fa, "wt") as f:
        f.write(">c1\nACGTACGT\n")
    assert list(read_fasta_streaming(fa)) == [("c1", "ACGTACGT")]
