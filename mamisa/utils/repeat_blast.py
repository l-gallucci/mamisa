"""
Self-BLAST repeat detection within contigs.

Blasts each contig against itself to identify internally duplicated regions.
A region covered by self-hits is likely a collapsed tandem or dispersed repeat.

Thresholds from the anvi'o long-read assembly benchmarking workflow:
  minimum alignment length : 200 bp
  minimum identity         : 80 %

The per-contig repeat statistics and the specific repeat intervals are returned
so that classify-clipping can check whether a clipping position falls inside
a repeat region and boost the repeat_collapse score accordingly.

Reference: Meren Lab long-read assembly benchmarking workflow
  https://merenlab.org/data/benchmarking-long-read-assemblers/
"""

import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .logging import log_info, log_warning, log_error
from .validation import check_dependencies


# ──────────────────────────────────────────────────────────────────────────────
# Type aliases
# ──────────────────────────────────────────────────────────────────────────────

# 0-based half-open intervals: (start, end)
Interval = Tuple[int, int]

RepeatStats = Dict   # keys: repeat_coverage_pct, n_repeat_regions, intervals


# ──────────────────────────────────────────────────────────────────────────────
# Dependency
# ──────────────────────────────────────────────────────────────────────────────

def check_blast() -> bool:
    deps = check_dependencies(['blastn', 'makeblastdb'])
    missing = [k for k, v in deps.items() if v is None]
    if missing:
        log_error(f"BLAST+ not found ({', '.join(missing)}). "
                  "Install: conda install -c bioconda blast")
        return False
    return True


# ──────────────────────────────────────────────────────────────────────────────
# Interval helpers
# ──────────────────────────────────────────────────────────────────────────────

def _merge_intervals(intervals: List[Interval]) -> List[Interval]:
    """Merge overlapping or adjacent intervals (0-based half-open)."""
    if not intervals:
        return []
    sorted_ivs = sorted(intervals)
    merged = [sorted_ivs[0]]
    for start, end in sorted_ivs[1:]:
        if start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def _covered_bases(intervals: List[Interval]) -> int:
    return sum(end - start for start, end in intervals)


def position_in_repeat(pos: int, intervals: List[Interval]) -> bool:
    """Return True if pos (0-based) falls within any repeat interval."""
    for start, end in intervals:
        if start <= pos < end:
            return True
    return False


# ──────────────────────────────────────────────────────────────────────────────
# Self-BLAST for one contig
# ──────────────────────────────────────────────────────────────────────────────

def self_blast_contig(
    contig_name: str,
    sequence: str,
    min_len: int = 200,
    min_identity: float = 80.0,
    tmp_dir: Optional[Path] = None,
) -> RepeatStats:
    """
    BLAST a contig against itself and return repeat statistics.

    Self-hits (query = subject at identical coordinates) are excluded.
    Reciprocal hits (A→B and B→A) are deduplicated by keeping only the
    hit where query_start < subject_start.

    Returns dict with:
      repeat_coverage_pct   float  % of contig covered by repeat intervals
      n_repeat_regions      int    number of distinct merged repeat intervals
      intervals             list   merged (start, end) 0-based half-open
    """
    seq_len = len(sequence)
    empty = {'repeat_coverage_pct': 0.0, 'n_repeat_regions': 0, 'intervals': []}

    if seq_len < min_len * 2:
        return empty

    use_tmp = tmp_dir is None
    if use_tmp:
        _tmp = tempfile.mkdtemp(prefix='mamisa_selfblast_')
        tmp_path = Path(_tmp)
    else:
        tmp_path = tmp_dir
        tmp_path.mkdir(parents=True, exist_ok=True)

    contig_fa = tmp_path / f'{contig_name[:60]}_self.fa'
    db_path   = tmp_path / f'{contig_name[:60]}_selfdb'

    try:
        # Write contig
        with open(contig_fa, 'w') as f:
            f.write(f">{contig_name}\n{sequence}\n")

        # Make db
        mkdb_cmd = [
            'makeblastdb', '-in', str(contig_fa),
            '-dbtype', 'nucl', '-out', str(db_path),
        ]
        subprocess.run(mkdb_cmd, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # Run blastn against itself
        blast_cmd = [
            'blastn',
            '-query',     str(contig_fa),
            '-db',        str(db_path),
            '-dust',      'no',
            '-outfmt',    '6 qstart qend sstart send pident length',
            '-out',       '-',
            '-perc_identity', str(min_identity),
        ]
        result = subprocess.run(blast_cmd, capture_output=True,
                                text=True, check=True)

        intervals: List[Interval] = []
        for line in result.stdout.splitlines():
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            try:
                qs, qe = int(parts[0]) - 1, int(parts[1])   # 1-based → 0-based
                ss, se = int(parts[2]) - 1, int(parts[3])
                aln_len = int(parts[5])
            except ValueError:
                continue

            # Skip exact self-alignment
            if qs == ss and qe == se:
                continue
            # Skip short alignments
            if aln_len < min_len:
                continue
            # Deduplicate reciprocals: keep only q→s where qs < ss
            if qs >= ss:
                continue

            # Both the query region and the subject region are repeat instances
            intervals.append((qs, qe))
            intervals.append((ss, se))

        merged = _merge_intervals(intervals)
        covered = _covered_bases(merged)
        pct = 100.0 * covered / seq_len if seq_len > 0 else 0.0

        return {
            'repeat_coverage_pct': round(pct, 2),
            'n_repeat_regions':    len(merged),
            'intervals':           merged,
        }

    except subprocess.CalledProcessError as e:
        log_warning(f"Self-BLAST failed for {contig_name}: {e}")
        return empty
    finally:
        if use_tmp:
            import shutil
            shutil.rmtree(_tmp, ignore_errors=True)


# ──────────────────────────────────────────────────────────────────────────────
# Batch self-BLAST for multiple contigs
# ──────────────────────────────────────────────────────────────────────────────

def batch_self_blast(
    contigs: Dict[str, str],   # {contig_name: sequence}
    min_len: int = 200,
    min_identity: float = 80.0,
    threads: int = 4,
) -> Dict[str, RepeatStats]:
    """
    Run self-BLAST for each contig in `contigs`.

    Threads are used only within each blastn call (one contig at a time).
    Returns {contig_name: RepeatStats}.
    """
    if not contigs:
        return {}

    log_info(f"Self-BLAST: processing {len(contigs):,} contig(s) "
             f"(min_len={min_len} bp, min_identity={min_identity}%)…")

    results: Dict[str, RepeatStats] = {}

    with tempfile.TemporaryDirectory(prefix='mamisa_selfblast_') as tmp:
        tmp_path = Path(tmp)
        for i, (name, seq) in enumerate(contigs.items(), 1):
            stats = self_blast_contig(
                name, seq,
                min_len=min_len,
                min_identity=min_identity,
                tmp_dir=tmp_path / f'contig_{i}',
            )
            results[name] = stats
            if stats['repeat_coverage_pct'] > 0:
                log_info(
                    f"  {name[:60]}: "
                    f"{stats['repeat_coverage_pct']:.1f}% repeat coverage, "
                    f"{stats['n_repeat_regions']} region(s)"
                )

    n_with_repeats = sum(1 for s in results.values() if s['repeat_coverage_pct'] > 0)
    log_info(f"Self-BLAST done: {n_with_repeats:,}/{len(results):,} contigs have repeats")
    return results
