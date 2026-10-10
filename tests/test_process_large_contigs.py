"""Regression tests for process-large-contigs file-name handling.

CheckM2 reports genomes by their individual-file stem, which already includes
the safe_filename hash. Extraction and assembly-rebuild must use that stem as
is, not re-apply safe_filename (which would double-hash and miss every file).
"""

from pathlib import Path

from mamisa.commands.process_large_contigs import (
    safe_filename,
    extract_hq_genomes,
    create_updated_assembly,
)


def test_safe_filename_appends_hash():
    s = safe_filename("7997BAurora_metamdbg_contig_23447")
    assert s.startswith("7997BAurora_metamdbg_contig_23447_")
    # 6-char hash suffix
    assert len(s.split("_")[-1]) == 6


def test_extract_hq_uses_stem_not_double_hash(tmp_path):
    individual = tmp_path / "02_individual"
    individual.mkdir()
    # file written exactly as split_large_contigs_individually would name it
    orig = "7997BAurora_metamdbg_contig_23447"
    stem = safe_filename(orig)                      # e.g. ..._db8c4c
    (individual / f"{stem}.fa").write_text(">c1\nACGT\n")

    # hq_list holds the CheckM2 'Name' == the stem (hash already baked in)
    out = tmp_path / "HQ_extracted"
    n = extract_hq_genomes(individual, [stem], out)

    assert n == 1
    assert (out / f"{stem}.fa").exists()


def test_create_updated_assembly_finds_kept_contig(tmp_path):
    individual = tmp_path / "02_individual"
    individual.mkdir()
    stem = safe_filename("contig_low_qual_1")
    (individual / f"{stem}.fa").write_text(">lowqual\nAAAA\n")

    regular = tmp_path / "regular.fa"
    regular.write_text(">reg1\nCCCC\n")

    out = tmp_path / "assembly_for_filtering.fa"
    create_updated_assembly(regular, individual, [stem], out)

    text = out.read_text()
    assert ">reg1" in text
    assert ">lowqual" in text   # kept large contig actually added
