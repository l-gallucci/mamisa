"""Tests for mamisa.utils.validation.check_db_path (per-tool path kind)."""

from mamisa.utils.validation import check_db_path


def test_gtdbtk_needs_directory_ok(tmp_path):
    d = tmp_path / "gtdbtk_data"
    d.mkdir()
    ok, msg = check_db_path(d, "gtdbtk")
    assert ok and msg == ""


def test_gtdbtk_rejects_file(tmp_path):
    f = tmp_path / "x.dmnd"
    f.touch()
    ok, msg = check_db_path(f, "gtdbtk")
    assert not ok and "DIRECTORY" in msg


def test_gunc_needs_file_ok(tmp_path):
    f = tmp_path / "gunc_db_progenomes2.1.dmnd"
    f.touch()
    ok, msg = check_db_path(f, "gunc")
    assert ok and msg == ""


def test_gunc_rejects_directory(tmp_path):
    d = tmp_path / "db"
    d.mkdir()
    ok, msg = check_db_path(d, "gunc")
    assert not ok and "FILE" in msg


def test_checkm2_file_ok(tmp_path):
    f = tmp_path / "uniref100.KO.1.dmnd"
    f.touch()
    ok, msg = check_db_path(f, "checkm2")
    assert ok and msg == ""


def test_non_dmnd_file_warns_but_ok(tmp_path):
    f = tmp_path / "db.bin"
    f.touch()
    ok, msg = check_db_path(f, "gunc")
    assert ok and ".dmnd" in msg      # warning, not failure


def test_missing_path(tmp_path):
    ok, msg = check_db_path(tmp_path / "nope", "gunc")
    assert not ok and "does not exist" in msg


def test_unknown_tool_passes(tmp_path):
    ok, msg = check_db_path(tmp_path, "unknown")
    assert ok and msg == ""
