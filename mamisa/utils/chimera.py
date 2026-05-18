"""
Chimera detection utilities for MAGs and circular contigs.

Detects potential chimeric assemblies using three orthogonal signals:
  1. GC content heterogeneity across contigs within a bin
  2. Windowed GC analysis along single circular contigs (junction detection)
  3. GTDB-Tk taxonomy consistency and confidence checks
"""

import csv
import math
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .logging import log_info, log_warning


# ---------------------------------------------------------------------------
# GC analysis
# ---------------------------------------------------------------------------

def compute_gc(sequence: str) -> float:
    """Return GC fraction (0.0–1.0). Returns NaN for sequences with no ACGT bases."""
    seq = sequence.upper()
    gc = seq.count('G') + seq.count('C')
    total = sum(seq.count(b) for b in 'ACGT')
    return gc / total if total > 0 else float('nan')


def windowed_gc(sequence: str, window: int = 5000, step: int = 2500) -> List[Tuple[int, float]]:
    """
    Compute GC content in a sliding window along a sequence.
    Returns list of (start_position, gc_fraction).
    Falls back to a single whole-sequence value when shorter than window.
    """
    seq_len = len(sequence)
    if seq_len < window:
        return [(0, compute_gc(sequence))]

    results = []
    for start in range(0, seq_len - window + 1, step):
        gc = compute_gc(sequence[start:start + window])
        if not math.isnan(gc):
            results.append((start, gc))
    return results


def gc_statistics(gc_values: List[float]) -> Dict[str, float]:
    """Summary statistics for a list of GC fractions (values in 0–1 range)."""
    valid = [v for v in gc_values if not math.isnan(v)]
    if not valid:
        return {'mean': float('nan'), 'std': float('nan'), 'cv': float('nan'),
                'min': float('nan'), 'max': float('nan'), 'delta': float('nan'), 'n': 0}

    n = len(valid)
    mean = sum(valid) / n
    variance = sum((v - mean) ** 2 for v in valid) / n if n > 1 else 0.0
    std = math.sqrt(variance)
    cv = std / mean if mean > 0 else 0.0

    return {
        'mean': mean, 'std': std, 'cv': cv,
        'min': min(valid), 'max': max(valid),
        'delta': max(valid) - min(valid),
        'n': n,
    }


def analyze_bin_gc(bin_fasta: Path, window: int = 5000, step: int = 2500) -> Dict:
    """
    Analyse GC content heterogeneity within a single bin FASTA file.

    Returns
    -------
    dict with:
      n_contigs            – number of contigs in the bin
      per_contig_gc        – {contig_name: gc_fraction}
      bin_gc_stats         – cross-contig GC summary statistics
      windowed_gc_stats    – windowed GC stats for single-contig circular candidates
      is_circular_candidate – True when n_contigs == 1 and length > 100 kbp
    """
    from .fasta import read_fasta_streaming

    contigs = list(read_fasta_streaming(bin_fasta))
    n_contigs = len(contigs)

    per_contig_gc = {name: compute_gc(seq) for name, seq in contigs}
    bin_stats = gc_statistics(list(per_contig_gc.values()))

    windowed_stats = None
    is_circular_candidate = False
    if n_contigs == 1:
        name, seq = contigs[0]
        is_circular_candidate = len(seq) > 100_000
        if is_circular_candidate:
            windows = windowed_gc(seq, window, step)
            windowed_stats = gc_statistics([gc for _, gc in windows])

    return {
        'n_contigs': n_contigs,
        'per_contig_gc': per_contig_gc,
        'bin_gc_stats': bin_stats,
        'windowed_gc_stats': windowed_stats,
        'is_circular_candidate': is_circular_candidate,
    }


# ---------------------------------------------------------------------------
# GTDB-Tk taxonomy parsing
# ---------------------------------------------------------------------------

def parse_gtdbtk_summary(summary_file: Path) -> Dict[str, Dict]:
    """
    Parse a GTDB-Tk summary TSV (bac120 or ar53).

    Returns dict: genome_name → taxonomy record.
    Handles both old (user_genome) and new (Name) column headers.
    """
    results = {}

    with open(summary_file, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            name = (row.get('user_genome') or row.get('Name') or '').strip()
            if not name:
                continue
            results[name] = {
                'classification': row.get('classification', '').strip(),
                'fastani_reference': row.get('fastani_reference', '').strip(),
                'fastani_ani': _safe_float(row.get('fastani_ani')),
                'fastani_af': _safe_float(row.get('fastani_af')),
                'msa_percent': _safe_float(row.get('msa_percent')),
                'red_value': _safe_float(row.get('red_value')),
                'warnings': row.get('warnings', '').strip(),
                'note': row.get('note', '').strip(),
            }

    log_info(f"Parsed GTDB-Tk taxonomy for {len(results):,} genomes from {summary_file.name}")
    return results


def _safe_float(value) -> Optional[float]:
    if value is None:
        return None
    s = str(value).strip()
    if s in ('', 'N/A', 'none', 'None', 'nan', 'N/a'):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def extract_taxonomy_level(classification: str, level: str = 'phylum') -> str:
    """
    Extract a specific rank from a GTDB semicolon-delimited classification string.
    e.g. 'd__Bacteria;p__Firmicutes_A;c__Bacilli;...' → 'Firmicutes_A'
    Returns 'Unclassified' when the rank is absent or empty.
    """
    prefixes = {
        'domain': 'd__', 'phylum': 'p__', 'class': 'c__',
        'order': 'o__', 'family': 'f__', 'genus': 'g__', 'species': 's__',
    }
    prefix = prefixes.get(level, 'p__')
    for part in classification.split(';'):
        part = part.strip()
        if part.startswith(prefix):
            taxon = part[len(prefix):]
            return taxon if taxon else 'Unclassified'
    return 'Unclassified'


def has_taxonomy_warning(tax_record: Dict) -> bool:
    """
    Return True when a GTDB-Tk record carries signals of unreliable placement:
      - non-empty warnings field
      - MSA percent < 10 % (poor marker gene recovery)
      - no FastANI reference match (placed by RED value only)
    """
    if tax_record.get('warnings'):
        return True
    msa = tax_record.get('msa_percent')
    if msa is not None and msa < 10.0:
        return True
    return False


# ---------------------------------------------------------------------------
# Risk assessment
# ---------------------------------------------------------------------------

def assess_chimera_risk(
    gc_delta: float,
    gc_cv: float,
    n_contigs: int,
    contamination: Optional[float],
    windowed_gc_delta: Optional[float],
    taxonomy_warning: bool = False,
) -> Tuple[str, List[str]]:
    """
    Assess chimera risk for a bin/MAG from multiple independent signals.

    Scoring:
      GC heterogeneity between contigs (multi-contig bins)
        delta > 10 %  → +3    delta > 5 %   → +1
      Windowed GC variation along a circular contig
        delta > 15 %  → +3    delta > 8 %   → +1
      CheckM2 contamination
        > 10 %        → +3    > 5 %         → +2
      GTDB-Tk placement warning
        any           → +1

    Risk levels  →  score thresholds
      High   ≥ 5
      Medium ≥ 2
      Low    ≥ 1
      Clean    0
    """
    reasons: List[str] = []
    score = 0

    if n_contigs > 1:
        if gc_delta > 0.10:
            score += 3
            reasons.append(f"High inter-contig GC heterogeneity (delta={gc_delta * 100:.1f}%)")
        elif gc_delta > 0.05:
            score += 1
            reasons.append(f"Moderate inter-contig GC heterogeneity (delta={gc_delta * 100:.1f}%)")

    if windowed_gc_delta is not None:
        if windowed_gc_delta > 0.15:
            score += 3
            reasons.append(
                f"High windowed GC variation along circular contig "
                f"(delta={windowed_gc_delta * 100:.1f}%)"
            )
        elif windowed_gc_delta > 0.08:
            score += 1
            reasons.append(
                f"Moderate windowed GC variation along circular contig "
                f"(delta={windowed_gc_delta * 100:.1f}%)"
            )

    if contamination is not None:
        if contamination > 10.0:
            score += 3
            reasons.append(f"High CheckM2 contamination ({contamination:.1f}%)")
        elif contamination > 5.0:
            score += 2
            reasons.append(f"Elevated CheckM2 contamination ({contamination:.1f}%)")

    if taxonomy_warning:
        score += 1
        reasons.append("GTDB-Tk placement warning (low confidence or poor MSA)")

    if score >= 5:
        return 'High', reasons
    elif score >= 2:
        return 'Medium', reasons
    elif score >= 1:
        return 'Low', reasons
    else:
        return 'Clean', reasons
