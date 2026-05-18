"""
Read-level taxonomy utilities for per-contig chimera detection.

Strategy
--------
1. Parse Kraken2 per-read output  → {read_id: taxid}
2. Parse Kraken2 report           → {taxid: (rank, name)} for human-readable output
3. Stream BAM via samtools        → for each mapped read: look up its taxid,
                                    accumulate per-contig and per-position counters
4. Compute per-contig profiles    → dominant taxon fraction, Shannon diversity,
                                    number of distinct taxa
5. Windowed analysis              → for long/circular contigs, slide a window across
                                    read positions and detect shifts in dominant taxon

No Python dependencies beyond the standard library.
Requires: samtools in PATH (for BAM streaming).
"""

import subprocess
import shutil
import math
from pathlib import Path
from collections import defaultdict, Counter
from typing import Dict, List, Optional, Tuple, Iterator

from .logging import log_info, log_warning, log_error
from .validation import check_dependencies


# ──────────────────────────────────────────────────────────────────────────────
# Kraken2 parsing
# ──────────────────────────────────────────────────────────────────────────────

class TaxRecord:
    __slots__ = ('classified', 'taxid', 'name', 'rank')

    def __init__(self, classified: bool, taxid: int,
                 name: str = '', rank: str = ''):
        self.classified = classified
        self.taxid = taxid
        self.name = name
        self.rank = rank


def parse_kraken2_output(kraken2_file: Path) -> Dict[str, int]:
    """
    Parse Kraken2 per-read output file (the --output file, not --report).

    Format (tab-separated):
        C/U  read_id  taxid  length  kmer_hits

    Returns {read_id: taxid}. Unclassified reads get taxid 0.
    """
    read_taxid: Dict[str, int] = {}

    with open(kraken2_file) as f:
        for lineno, line in enumerate(f, start=1):
            if not line.strip():
                continue
            fields = line.split('\t')
            if len(fields) < 3:
                log_warning(
                    f"Skipping malformed Kraken2 output line {lineno} "
                    f"(expected ≥3 fields, got {len(fields)})"
                )
                continue
            status = fields[0].strip()
            read_id = fields[1].strip()
            try:
                taxid = int(fields[2].strip())
            except ValueError:
                log_warning(f"Non-integer taxid on line {lineno}, skipping")
                continue
            if status == 'U':
                taxid = 0
            read_taxid[read_id] = taxid

    log_info(f"Loaded {len(read_taxid):,} read classifications from {kraken2_file.name}")
    return read_taxid


def parse_kraken2_report(report_file: Path) -> Dict[int, Tuple[str, str]]:
    """
    Parse a Kraken2 report file (the --report file).

    Format (tab-separated):
        %reads  n_clade  n_direct  rank_code  taxid  name(may be indented)

    Rank codes: U, R, D, P, C, O, F, G, S, S1, …

    Returns {taxid: (rank_code, name)}.
    """
    taxid_info: Dict[int, Tuple[str, str]] = {
        0: ('U', 'unclassified'),
        1: ('R', 'root'),
    }

    with open(report_file) as f:
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 6:
                continue
            rank = fields[3].strip()
            try:
                taxid = int(fields[4].strip())
            except ValueError:
                continue
            name = fields[5].strip()
            taxid_info[taxid] = (rank, name)

    log_info(f"Loaded taxonomy names for {len(taxid_info):,} taxa from {report_file.name}")
    return taxid_info


def taxon_label(taxid: int, taxid_info: Dict[int, Tuple[str, str]]) -> str:
    """Return a human-readable label for a taxid."""
    if taxid == 0:
        return 'unclassified'
    info = taxid_info.get(taxid)
    if info:
        rank, name = info
        return f"{name} [{rank}:{taxid}]"
    return f"taxid:{taxid}"


# ──────────────────────────────────────────────────────────────────────────────
# BAM streaming via samtools
# ──────────────────────────────────────────────────────────────────────────────

def check_samtools() -> bool:
    """Return True if samtools is available in PATH."""
    deps = check_dependencies(['samtools'])
    if deps['samtools'] is None:
        log_error("samtools not found in PATH — required for BAM parsing")
        log_error("Install via: conda install -c bioconda samtools")
        return False
    log_info(f"Found samtools: {deps['samtools']}")
    return True


def stream_bam(
    bam_file: Path,
    min_mapq: int = 20,
) -> Iterator[Tuple[str, str, int]]:
    """
    Stream a BAM file via ``samtools view``.

    Yields (read_id, contig_name, position) for every primary, mapped read
    that passes the MAPQ threshold.

    Flags filter: -F 2308 excludes unmapped (4), secondary (256),
    supplementary (2048) — keeping only primary alignments.
    """
    cmd = [
        'samtools', 'view',
        '-F', '2308',          # exclude unmapped, secondary, supplementary
        '-q', str(min_mapq),
        str(bam_file),
    ]

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
    )

    try:
        for line in proc.stdout:
            if line.startswith('@'):
                continue
            fields = line.split('\t')
            if len(fields) < 4:
                continue
            read_id = fields[0]
            contig = fields[2]
            if contig == '*':
                continue
            try:
                pos = int(fields[3])
            except ValueError:
                continue
            yield read_id, contig, pos
    finally:
        proc.stdout.close()
        proc.wait()
        if proc.returncode not in (0, None, -13):   # -13 = SIGPIPE (early close)
            stderr = proc.stderr.read()
            log_warning(f"samtools exited with code {proc.returncode}: {stderr[:200]}")


# ──────────────────────────────────────────────────────────────────────────────
# Per-contig aggregation
# ──────────────────────────────────────────────────────────────────────────────

class ContigReadProfile:
    """Accumulates read-taxonomy data for one contig."""

    def __init__(self, name: str, length: int = 0, window_threshold: int = 50_000):
        self.name = name
        self.length = length
        self.window_threshold = window_threshold

        self.taxid_counts: Counter = Counter()    # taxid → n_reads
        self.positions: List[Tuple[int, int]] = [] # (pos, taxid) for windowed analysis

    def add_read(self, pos: int, taxid: int):
        self.taxid_counts[taxid] += 1
        if self.length >= self.window_threshold:
            self.positions.append((pos, taxid))

    @property
    def n_reads(self) -> int:
        return sum(self.taxid_counts.values())


def compute_taxonomy_stats(
    taxid_counts: Counter,
    taxid_info: Dict[int, Tuple[str, str]],
    exclude_unclassified: bool = False,
) -> Dict:
    """
    Compute diversity statistics from a {taxid: count} counter.

    Returns a dict with:
      n_reads_total, n_classified, n_taxa (excl. unclassified),
      dominant_taxid, dominant_name, dominant_fraction,
      shannon_diversity
    """
    counts = dict(taxid_counts)
    n_total = sum(counts.values())
    n_unclassified = counts.get(0, 0)
    n_classified = n_total - n_unclassified

    # For diversity, optionally exclude unclassified
    analysis_counts = {k: v for k, v in counts.items() if k != 0} \
        if exclude_unclassified else counts

    n_analysis = sum(analysis_counts.values())
    n_taxa = len([k for k in analysis_counts if k != 0])

    dominant_taxid = 0
    dominant_fraction = 0.0
    if analysis_counts:
        dominant_taxid = max(analysis_counts, key=analysis_counts.get)
        dominant_fraction = (
            analysis_counts[dominant_taxid] / n_analysis if n_analysis > 0 else 0.0
        )

    # Shannon diversity H = -Σ p_i * ln(p_i)
    shannon = 0.0
    if n_analysis > 0:
        for count in analysis_counts.values():
            p = count / n_analysis
            if p > 0:
                shannon -= p * math.log(p)

    return {
        'n_reads_total': n_total,
        'n_classified': n_classified,
        'classified_fraction': n_classified / n_total if n_total > 0 else 0.0,
        'n_taxa': n_taxa,
        'dominant_taxid': dominant_taxid,
        'dominant_name': taxon_label(dominant_taxid, taxid_info),
        'dominant_fraction': dominant_fraction,
        'shannon_diversity': shannon,
    }


# ──────────────────────────────────────────────────────────────────────────────
# Windowed analysis for long / circular contigs
# ──────────────────────────────────────────────────────────────────────────────

def windowed_taxon_profile(
    positions: List[Tuple[int, int]],  # (mapping_pos, taxid)
    contig_length: int,
    window: int,
    step: int,
    taxid_info: Dict[int, Tuple[str, str]],
    min_reads_per_window: int = 5,
) -> List[Dict]:
    """
    Slide a window along the contig and compute the dominant taxon in each window.

    Returns a list of window dicts:
        start, end, n_reads, dominant_taxid, dominant_name,
        dominant_fraction, n_taxa, is_shift (dominant taxon changed vs prev window)
    """
    if not positions or contig_length == 0:
        return []

    windows = []
    prev_dominant = None

    for start in range(0, max(contig_length - window + 1, 1), step):
        end = min(start + window, contig_length)

        window_reads = [taxid for pos, taxid in positions if start <= pos < end]
        if len(window_reads) < min_reads_per_window:
            continue

        counter = Counter(window_reads)
        analysis = Counter({k: v for k, v in counter.items() if k != 0})
        n_analysis = sum(analysis.values())

        if n_analysis == 0:
            continue

        dom_taxid = max(analysis, key=analysis.get)
        dom_fraction = analysis[dom_taxid] / n_analysis
        n_taxa = len([k for k in analysis if k != 0])

        is_shift = (prev_dominant is not None) and (dom_taxid != prev_dominant)
        prev_dominant = dom_taxid

        windows.append({
            'start': start,
            'end': end,
            'n_reads': len(window_reads),
            'dominant_taxid': dom_taxid,
            'dominant_name': taxon_label(dom_taxid, taxid_info),
            'dominant_fraction': dom_fraction,
            'n_taxa': n_taxa,
            'is_shift': is_shift,
        })

    return windows


def detect_taxonomic_shifts(windows: List[Dict]) -> List[Dict]:
    """Return only the windows where the dominant taxon shifted from the previous window."""
    return [w for w in windows if w.get('is_shift')]


# ──────────────────────────────────────────────────────────────────────────────
# Risk assessment
# ──────────────────────────────────────────────────────────────────────────────

def assess_read_chimera_risk(
    dominant_fraction: float,
    n_taxa: int,
    n_reads: int,
    n_taxon_shifts: int,
    classified_fraction: float,
) -> Tuple[str, List[str]]:
    """
    Assess chimera risk for one contig from its read-taxonomy profile.

    Scoring:
      dominant_fraction < 0.50  → +4  (reads split roughly evenly between taxa)
      dominant_fraction < 0.65  → +3
      dominant_fraction < 0.80  → +2
      dominant_fraction < 0.90  → +1

      n_taxa > 3                → +2
      n_taxa > 1                → +1

      taxon_shifts > 2          → +3  (multiple jumps in dominant taxon along contig)
      taxon_shifts > 0          → +2  (at least one position where dominant taxon changes)

    Risk levels: High ≥ 5 | Medium ≥ 2 | Low ≥ 1 | Clean 0
    Minimum 10 reads required; below that → 'Insufficient'
    """
    if n_reads < 10:
        return 'Insufficient', [f"Too few reads ({n_reads}) for reliable classification"]

    reasons: List[str] = []
    score = 0

    if dominant_fraction < 0.50:
        score += 4
        reasons.append(
            f"Reads evenly split across taxa — dominant taxon covers only "
            f"{dominant_fraction * 100:.1f}% of classified reads"
        )
    elif dominant_fraction < 0.65:
        score += 3
        reasons.append(
            f"Low dominant taxon fraction ({dominant_fraction * 100:.1f}%)"
        )
    elif dominant_fraction < 0.80:
        score += 2
        reasons.append(
            f"Moderate dominant taxon fraction ({dominant_fraction * 100:.1f}%)"
        )
    elif dominant_fraction < 0.90:
        score += 1
        reasons.append(
            f"Slightly reduced dominant taxon fraction ({dominant_fraction * 100:.1f}%)"
        )

    if n_taxa > 3:
        score += 2
        reasons.append(f"High taxon richness per contig ({n_taxa} distinct taxa)")
    elif n_taxa > 1:
        score += 1
        reasons.append(f"Multiple taxa detected ({n_taxa})")

    if n_taxon_shifts > 2:
        score += 3
        reasons.append(
            f"Dominant taxon shifts at {n_taxon_shifts} positions along the contig — "
            "strong chimeric junction signal"
        )
    elif n_taxon_shifts > 0:
        score += 2
        reasons.append(
            f"Dominant taxon shifts at {n_taxon_shifts} position(s) along the contig"
        )

    if classified_fraction < 0.30 and n_reads >= 50:
        reasons.append(
            f"Low read classification rate ({classified_fraction * 100:.1f}%) — "
            "may indicate novel or highly divergent sequence"
        )

    if score >= 5:
        return 'High', reasons
    elif score >= 2:
        return 'Medium', reasons
    elif score >= 1:
        return 'Low', reasons
    else:
        return 'Clean', reasons
