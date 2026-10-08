"""Tests for setup-workflow / check-envs pure helpers."""

from mamisa.commands.setup_workflow import set_yaml_scalar
from mamisa.commands.check_envs import _parse_ver, check_env


def test_set_yaml_scalar_replaces_existing():
    text = "threads: 1\nfoo: bar\n"
    out = set_yaml_scalar(text, "threads", 40)
    assert "threads: 40" in out
    assert "foo: bar" in out            # other lines preserved


def test_set_yaml_scalar_quotes():
    out = set_yaml_scalar("genomes_dir: x\n", "genomes_dir", "/a b/c", quote=True)
    assert 'genomes_dir: "/a b/c"' in out


def test_set_yaml_scalar_appends_when_missing():
    out = set_yaml_scalar("threads: 1\n", "env_gunc", "gunc", quote=True)
    assert 'env_gunc: "gunc"' in out


def test_set_yaml_scalar_only_top_level():
    # an indented key of the same name must not be touched
    text = "block:\n  threads: 1\nthreads: 1\n"
    out = set_yaml_scalar(text, "threads", 8)
    assert "  threads: 1" in out        # nested untouched
    assert "threads: 8" in out


def test_parse_ver():
    assert _parse_ver("samtools 1.21") == (1, 21, 0)
    assert _parse_ver("diamond version 2.1.24") == (2, 1, 24)
    assert _parse_ver("no numbers here") is None


def test_check_env_unreachable(monkeypatch):
    # when conda run fails, env check returns not ok without raising
    import mamisa.commands.check_envs as ce
    monkeypatch.setattr(ce, "_conda_run", lambda env, cmd, timeout=120: (1, "no env"))
    result = ce.check_env("ghost", "gunc", repo="/repo")
    assert result["ok"] is False
    assert result["env"] == "ghost"
