"""
Misassembly data parsing utilities (anvi'o clipping files)
"""

from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Set, Tuple

from .logging import log_info, log_warning


# A clipping record: (position, clipping_coverage)
ClipRecord = Tuple[int, int]


def find_clipping_files(misassembly_dir: Path) -> List[Path]:
    """Find all *-clipping.txt files in a directory."""
    files = sorted(misassembly_dir.glob('*-clipping.txt'))
    if not files:
        log_warning(f"No *-clipping.txt files found in {misassembly_dir}")
    return files


def parse_clipping_files(clipping_files: List[Path]) -> Dict[str, List[ClipRecord]]:
    """
    Parse *-clipping.txt files produced by anvi-script-find-misassemblies.

    Expected column layout (tab-separated, one header line):
        contig_name  sample_name  position  clipping_coverage  [...]

    Returns:
        dict mapping contig_name -> sorted list of (position, coverage) tuples.
        If coverage is absent or non-numeric, it defaults to 1 so downstream
        filters can still work without crashing.
    """
    clipping_data: Dict[str, List[ClipRecord]] = defaultdict(list)

    for file_path in clipping_files:
        log_info(f"  Reading: {file_path.name}")
        with open(file_path) as f:
            f.readline()  # skip header
            for lineno, line in enumerate(f, start=2):
                if not line.strip():
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 3:
                    log_warning(
                        f"Skipping malformed line {lineno} in {file_path.name} "
                        f"(expected ≥3 fields, got {len(fields)})"
                    )
                    continue

                contig_name = fields[0]

                try:
                    position = int(fields[2])
                except ValueError:
                    log_warning(
                        f"Skipping line {lineno} in {file_path.name}: "
                        f"position '{fields[2]}' is not an integer"
                    )
                    continue

                # Field 4 (index 3) is the clipping read count
                coverage = 1
                if len(fields) >= 4:
                    try:
                        coverage = int(fields[3])
                    except ValueError:
                        pass  # non-numeric coverage → keep default 1

                clipping_data[contig_name].append((position, coverage))

    # Sort each contig's records by position
    for contig in clipping_data:
        clipping_data[contig].sort(key=lambda r: r[0])

    return dict(clipping_data)


def filter_by_coverage(
    clipping_data: Dict[str, List[ClipRecord]],
    min_coverage: int,
) -> Dict[str, List[ClipRecord]]:
    """
    Return a filtered copy keeping only clipping records where coverage >= min_coverage.
    Contigs with no records surviving the filter are dropped entirely.
    """
    if min_coverage <= 1:
        return clipping_data

    filtered: Dict[str, List[ClipRecord]] = {}
    for contig, records in clipping_data.items():
        passing = [r for r in records if r[1] >= min_coverage]
        if passing:
            filtered[contig] = passing

    n_removed = sum(len(v) for v in clipping_data.values()) - \
                sum(len(v) for v in filtered.values())
    log_info(
        f"Coverage filter (≥{min_coverage}): "
        f"kept {sum(len(v) for v in filtered.values()):,} positions, "
        f"dropped {n_removed:,} low-coverage positions"
    )
    return filtered


def load_clipping_data(
    misassembly_dir: Path,
    min_coverage: int = 0,
) -> Dict[str, List[ClipRecord]]:
    """
    Load all clipping data from a directory, optionally filtered by coverage.

    Returns dict mapping contig_name -> sorted list of (position, coverage).
    """
    files = find_clipping_files(misassembly_dir)
    if not files:
        return {}

    data = parse_clipping_files(files)
    if min_coverage > 1:
        data = filter_by_coverage(data, min_coverage)

    n_contigs = len(data)
    n_positions = sum(len(v) for v in data.values())
    log_info(f"Clipping data: {n_contigs:,} contigs, {n_positions:,} positions")
    return data


def load_clipping_positions(
    misassembly_dir: Path,
    min_coverage: int = 0,
) -> Dict[str, List[int]]:
    """
    Load clipping positions (position-only, no coverage) from a directory.
    Backward-compatible wrapper around load_clipping_data.
    """
    data = load_clipping_data(misassembly_dir, min_coverage)
    return {contig: [pos for pos, _ in records] for contig, records in data.items()}


def get_contigs_with_misassemblies(
    misassembly_dir: Path,
    min_coverage: int = 0,
) -> Set[str]:
    """Return the set of contig names that have at least one clipping position."""
    return set(load_clipping_positions(misassembly_dir, min_coverage).keys())
