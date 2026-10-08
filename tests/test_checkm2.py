"""Tests for mamisa.utils.checkm2."""

import math
import pytest

from mamisa.utils.checkm2 import (
    parse_quality_value,
    assign_tier,
    parse_quality_report,
    parse_all_reports,
    NoReportsFoundError,
)


@pytest.mark.parametrize("raw,expected", [
    ("95.0", 95.0),
    (" 88.5 ", 88.5),
    ("90%", 90.0),
    (None, None),       # NaN sentinel handled below
    ("n/a", None),
])
def test_parse_quality_value(raw, expected):
    val = parse_quality_value(raw)
    if expected is None:
        assert math.isnan(val)
    else:
        assert val == expected


@pytest.mark.parametrize("comp,cont,tier", [
    (98.0, 1.0, "HQ"),
    (90.0, 5.0, "HQ"),       # boundary inclusive
    (89.9, 1.0, "MQ"),       # just under HQ completeness
    (90.0, 5.1, "MQ"),       # just over HQ contamination
    (70.0, 10.0, "MQ"),
    (55.0, 4.0, "LQ"),
    (50.0, 10.0, "LQ"),
    (49.9, 1.0, "Fail"),
    (80.0, 20.0, "Fail"),
])
def test_assign_tier(comp, cont, tier):
    assert assign_tier(comp, cont) == tier


def _report(path, rows):
    header = "Name\tCompleteness\tContamination\n"
    body = "".join(f"{n}\t{c}\t{x}\n" for n, c, x in rows)
    path.write_text(header + body)
    return path


def test_parse_quality_report(tmp_path):
    f = _report(tmp_path / "quality_report.tsv",
                [("binA", "98.5", "1.2"), ("binB", "70", "8")])
    recs = parse_quality_report(f)
    assert [r["name"] for r in recs] == ["binA", "binB"]
    assert recs[0]["completeness"] == 98.5
    assert recs[1]["contamination"] == 8.0


def test_parse_quality_report_missing_columns(tmp_path):
    bad = tmp_path / "quality_report.tsv"
    bad.write_text("Foo\tBar\n1\t2\n")
    assert parse_quality_report(bad) == []


def test_parse_all_reports_assigns_tiers(tmp_path):
    d = tmp_path / "checkm2"
    d.mkdir()
    _report(d / "quality_report.tsv",
            [("binA", "98", "1"), ("binB", "72", "9"), ("binC", "10", "1")])
    thresholds = {}
    records, counts = parse_all_reports(d, thresholds)
    tiers = {r["name"]: r["tier"] for r in records}
    assert tiers == {"binA": "HQ", "binB": "MQ", "binC": "Fail"}
    assert counts["HQ"] == 1 and counts["MQ"] == 1 and counts["Fail"] == 1


def test_parse_all_reports_raises_when_empty(tmp_path):
    empty = tmp_path / "nope"
    empty.mkdir()
    with pytest.raises(NoReportsFoundError):
        parse_all_reports(empty, {})
