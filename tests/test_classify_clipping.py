"""Tests for mamisa.commands.classify_clipping.classify_position."""

from mamisa.commands.classify_clipping import classify_position


def _stats(**kw):
    base = dict(
        near_contig_end=False,
        depth_ratio=1.0,
        discordant_fraction=0.0,
        large_insert_fraction=0.0,
        strand_fwd_fraction=0.5,
        clipped_base_entropy=2.0,
        primary_reads_in_window=50,
        sa_partner_contig="N/A",
        sa_partner_support=0,
        sa_partner_fraction="N/A",
    )
    base.update(kw)
    return base


def test_end_artefact():
    label, conf, _ = classify_position(_stats(near_contig_end=True))
    assert label == "end_artefact" and conf == "High"


def test_low_confidence_few_reads():
    label, _, _ = classify_position(_stats(primary_reads_in_window=2))
    assert label == "low_confidence"


def test_repeat_collapse_depth_spike():
    label, _, ev = classify_position(_stats(depth_ratio=4.0, clipped_base_entropy=0.5))
    assert label == "repeat_collapse"
    assert any("depth_spike" in e for e in ev)


def test_deletion_artefact_depth_drop():
    label, _, _ = classify_position(_stats(depth_ratio=0.1))
    assert label == "deletion_artefact"


def test_chimeric_join_from_sa_tag():
    label, conf, ev = classify_position(_stats(
        sa_partner_contig="contigB", sa_partner_support=8, sa_partner_fraction=0.8))
    assert label == "chimeric_join"
    assert any("sa_join_contigB" in e for e in ev)


def test_chimeric_join_requires_min_support():
    # only 2 clipped reads with SA -> below support>=3 gate, no join
    label, _, _ = classify_position(_stats(
        sa_partner_contig="contigB", sa_partner_support=2, sa_partner_fraction=1.0))
    assert label != "chimeric_join"


def test_chimeric_join_beats_other_signals():
    # strong join + some repeat signal -> join wins on priority
    label, _, _ = classify_position(_stats(
        depth_ratio=4.0,
        sa_partner_contig="contigB", sa_partner_support=10, sa_partner_fraction=0.9))
    assert label == "chimeric_join"
