"""
Misassembly data parsing utilities (anvi'o clipping files)
"""

from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Set

from .logging import log_info, log_warning


def find_clipping_files(misassembly_dir: Path) -> List[Path]:
    """Find all *-clipping.txt files in a directory."""
    files = sorted(misassembly_dir.glob('*-clipping.txt'))
    if not files:
        log_warning(f"No *-clipping.txt files found in {misassembly_dir}")
    return files


def parse_clipping_files(clipping_files: List[Path]) -> Dict[str, List[int]]:
    """
    Parse *-clipping.txt files from anvi'o.

    Returns:
        dict mapping contig_name -> sorted list of clipping positions
    """
    clipping_positions: Dict[str, List[int]] = defaultdict(list)

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
                clipping_positions[contig_name].append(position)

    for contig in clipping_positions:
        clipping_positions[contig].sort()

    return dict(clipping_positions)


def load_clipping_positions(misassembly_dir: Path) -> Dict[str, List[int]]:
    """
    Load all clipping positions from a directory.
    Returns dict mapping contig_name -> sorted list of clipping positions.
    """
    files = find_clipping_files(misassembly_dir)
    if not files:
        return {}
    positions = parse_clipping_files(files)
    log_info(f"Contigs with misassemblies: {len(positions):,}")
    return positions


def get_contigs_with_misassemblies(misassembly_dir: Path) -> Set[str]:
    """
    Return only the set of contig names that have misassemblies.
    Convenience wrapper around load_clipping_positions when positions are not needed.
    """
    return set(load_clipping_positions(misassembly_dir).keys())
