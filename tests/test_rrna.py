"""Tests for MIMAG rRNA/tRNA parsing (mamisa.utils.rrna) and the HQ gate."""

from pathlib import Path

import pytest

from mamisa.utils.rrna import (
    parse_barrnap_gff,
    parse_trnascan,
    assess_mimag_rna,
    write_mimag_rna_summary,
    parse_mimag_rna_summary,
    STANDARD_AA,
    MIMAG_MIN_AA,
)
from mamisa.commands.filter_checkm2 import apply_mimag_rna_gate


BARRNAP_GFF = """\
##gff-version 3
contig1\tbarrnap:0.9\trRNA\t100\t1638\t.\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA
contig1\tbarrnap:0.9\trRNA\t2000\t2120\t.\t+\t.\tName=5S_rRNA;product=5S ribosomal RNA
contig2\tbarrnap:0.9\trRNA\t50\t2950\t.\t-\t.\tName=23S_rRNA;product=23S ribosomal RNA
contig3\tbarrnap:0.9\trRNA\t10\t400\t.\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA (partial)
"""


def test_parse_barrnap_counts(tmp_path):
    gff = tmp_path / "g.gff"
    gff.write_text(BARRNAP_GFF)
    counts = parse_barrnap_gff(gff)
    # one full 16S (the second is partial -> excluded), one 5S, one 23S
    assert counts == {'5S': 1, '16S': 1, '23S': 1}


def test_parse_barrnap_missing_file(tmp_path):
    assert parse_barrnap_gff(tmp_path / "nope.gff") == {'5S': 0, '16S': 0, '23S': 0}


def _trnascan_table(aas):
    header = "Name\ttRNA#\tBegin\tEnd\tType\tCodon\tBegin\tEnd\tScore\tNote\n"
    header += "--------\t------\t-----\t------\t----\t-----\t----\t----\t------\t----\n"
    header += "=\t=\t=\t=\t=\t=\t=\t=\t=\t=\n"
    rows = ""
    for i, aa in enumerate(aas, 1):
        rows += f"contig1\t{i}\t100\t180\t{aa}\tCTT\t0\t0\t72.3\t\n"
    return header + rows


def test_parse_trnascan_distinct_aa(tmp_path):
    out = tmp_path / "t.tsv"
    out.write_text(_trnascan_table(['Lys', 'Lys', 'Arg', 'Undet', 'Pseudo']))
    aas = parse_trnascan(out)
    assert aas == {'Lys', 'Arg'}  # dedup; Undet/Pseudo excluded


def test_parse_trnascan_pseudo_note_excluded(tmp_path):
    out = tmp_path / "t.tsv"
    table = "Name\ttRNA#\tBegin\tEnd\tType\tCodon\tB\tE\tScore\tNote\n"
    table += "-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n"
    table += "=\t=\t=\t=\t=\t=\t=\t=\t=\t=\n"
    table += "c1\t1\t1\t80\tMet\tCAT\t0\t0\t50\tpseudo\n"
    out.write_text(table)
    assert parse_trnascan(out) == set()


def test_assess_mimag_rna_full_pass():
    rrna = {'5S': 1, '16S': 1, '23S': 1}
    trna = set(list(STANDARD_AA)[:MIMAG_MIN_AA])
    v = assess_mimag_rna(rrna, trna)
    assert v['rrna_complete'] is True
    assert v['trna_complete'] is True
    assert v['mimag_rna_ok'] is True


def test_assess_mimag_rna_fails_on_missing_23s():
    rrna = {'5S': 1, '16S': 1, '23S': 0}
    trna = set(STANDARD_AA)
    v = assess_mimag_rna(rrna, trna)
    assert v['rrna_complete'] is False
    assert v['mimag_rna_ok'] is False


def test_assess_mimag_rna_fails_on_few_trna():
    rrna = {'5S': 1, '16S': 1, '23S': 1}
    trna = set(list(STANDARD_AA)[:MIMAG_MIN_AA - 1])  # 17 < 18
    v = assess_mimag_rna(rrna, trna)
    assert v['trna_complete'] is False
    assert v['mimag_rna_ok'] is False


def test_summary_round_trip(tmp_path):
    recs = [
        {'genome': 'A.fa', 'n_5S': 1, 'n_16S': 1, 'n_23S': 1, 'n_trna_aa': 20,
         'rrna_complete': True, 'trna_complete': True, 'mimag_rna_ok': True},
        {'genome': 'B.fa', 'n_5S': 0, 'n_16S': 1, 'n_23S': 1, 'n_trna_aa': 10,
         'rrna_complete': False, 'trna_complete': False, 'mimag_rna_ok': False},
    ]
    path = tmp_path / "mimag_rna_summary.tsv"
    write_mimag_rna_summary(recs, path)
    loaded = parse_mimag_rna_summary(path)
    assert loaded['A.fa']['mimag_rna_ok'] is True
    assert loaded['B.fa']['mimag_rna_ok'] is False
    assert loaded['A.fa']['n_trna_aa'] == 20


def test_gate_demotes_failing_hq():
    records = [
        {'name': 'A', 'tier': 'HQ', 'completeness': 95, 'contamination': 2},
        {'name': 'B', 'tier': 'HQ', 'completeness': 93, 'contamination': 1},
        {'name': 'C', 'tier': 'MQ', 'completeness': 80, 'contamination': 3},
    ]
    mimag = {
        'A': {'mimag_rna_ok': True},
        'B': {'mimag_rna_ok': False},   # -> demote to MQ
    }
    counts = apply_mimag_rna_gate(records, mimag)
    assert records[0]['tier'] == 'HQ'
    assert records[1]['tier'] == 'MQ'
    assert records[1]['hq_rna_demoted'] is True
    assert counts['HQ'] == 1
    assert counts['MQ'] == 2


def test_gate_missing_record_keeps_hq():
    records = [{'name': 'X', 'tier': 'HQ', 'completeness': 99, 'contamination': 0}]
    counts = apply_mimag_rna_gate(records, {})  # no record for X
    assert records[0]['tier'] == 'HQ'
    assert records[0]['mimag_rna_ok'] == ''
    assert counts['HQ'] == 1
