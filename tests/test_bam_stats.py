"""Tests for mamisa.utils.bam_stats (SA-tag join detection, clipped bases)."""

import math

from mamisa.utils.bam_stats import _PosAccumulator, _extract_clipped_bases


def test_extract_clipped_bases_leading_softclip():
    # 50S100M read aligned at ref pos 1000 -> soft-clip boundary at ref 1000
    seq = "A" * 150
    clipped = _extract_clipped_bases("50S100M", seq, read_pos=1000, clip_pos=1000)
    assert clipped == "A" * 50


def test_extract_clipped_bases_no_match():
    assert _extract_clipped_bases("100M", "A" * 100, 1000, 1000) == ""


def test_sa_join_dominant_partner():
    acc = _PosAccumulator(1000, this_contig="contigA")
    for i in range(10):
        sa = ["contigB"] if i < 8 else None
        acc.add_read(flag=0, tlen=0, seq="A" * 150, cigar="50S100M",
                     read_pos=1000, insert_mean=300, insert_std=50, sa_partners=sa)
    partner, support, n_clip = acc.sa_join_stats()
    assert partner == "contigB"
    assert support == 8
    assert n_clip == 10


def test_sa_join_ignores_self_contig():
    acc = _PosAccumulator(1000, this_contig="contigA")
    for _ in range(5):
        acc.add_read(flag=0, tlen=0, seq="A" * 150, cigar="50S100M",
                     read_pos=1000, insert_mean=0, insert_std=0,
                     sa_partners=["contigA"])   # self -> not a join
    partner, support, _ = acc.sa_join_stats()
    assert partner == "N/A" and support == 0


def test_to_dict_exposes_sa_columns():
    acc = _PosAccumulator(1000, this_contig="contigA")
    for _ in range(6):
        acc.add_read(flag=0, tlen=0, seq="A" * 150, cigar="50S100M",
                     read_pos=1000, insert_mean=0, insert_std=0,
                     sa_partners=["contigB"])
    d = acc.to_dict("contigA", 50000, 30.0)
    assert d["sa_partner_contig"] == "contigB"
    assert d["sa_partner_support"] == 6
    assert d["n_clipped_reads"] == 6
    assert d["sa_partner_fraction"] == 1.0


def test_secondary_reads_not_counted_for_pairs():
    acc = _PosAccumulator(1000, this_contig="contigA")
    # secondary alignment (flag 0x100) counts for depth but not primary stats
    acc.add_read(flag=0x100, tlen=0, seq="A" * 150, cigar="50S100M",
                 read_pos=1000, insert_mean=0, insert_std=0)
    assert acc.n_depth == 1
    assert acc.n_primary == 0
