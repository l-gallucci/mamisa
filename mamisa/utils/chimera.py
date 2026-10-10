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
    per_contig_len = {name: len(seq) for name, seq in contigs}
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
        'per_contig_len': per_contig_len,
        'bin_gc_stats': bin_stats,
        'windowed_gc_stats': windowed_stats,
        'is_circular_candidate': is_circular_candidate,
    }


def gc_outlier_contigs(per_contig_gc: Dict[str, float],
                       per_contig_len: Dict[str, int]) -> Optional[Dict]:
    """
    Identify the two contigs that set the inter-contig GC delta (the max-GC and
    min-GC contigs), so a reviewer sees WHICH contigs drive a flag and how long
    they are. A short outlier is a weaker reason to worry than a long one.

    Returns None for single-contig bins (no inter-contig delta) or when GC is
    undefined. Otherwise {high: (name, gc, len), low: (name, gc, len)}.
    """
    valid = {n: g for n, g in per_contig_gc.items() if not math.isnan(g)}
    if len(valid) < 2:
        return None
    high_name = max(valid, key=valid.get)
    low_name = min(valid, key=valid.get)
    return {
        'high': (high_name, valid[high_name], per_contig_len.get(high_name, 0)),
        'low': (low_name, valid[low_name], per_contig_len.get(low_name, 0)),
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
            # GTDB-Tk >=2.4 renamed FastANI columns to closest_genome_* (skani).
            # Accept both so old and new summary files parse.
            results[name] = {
                'classification': row.get('classification', '').strip(),
                'closest_reference': (row.get('closest_genome_reference')
                                      or row.get('fastani_reference') or '').strip(),
                'closest_ani': _safe_float(row.get('closest_genome_ani')
                                           or row.get('fastani_ani')),
                'closest_af': _safe_float(row.get('closest_genome_af')
                                          or row.get('fastani_af')),
                'msa_percent': _safe_float(row.get('msa_percent')),
                'red_value': _safe_float(row.get('red_value')),
                'warnings': row.get('warnings', '').strip(),
                'note': row.get('note', '').strip(),
            }

    log_info(f"Parsed GTDB-Tk taxonomy for {len(results):,} genomes from {summary_file.name}")
    return results


def parse_gunc_output(gunc_dir: Path) -> Dict[str, Dict]:
    """
    Parse GUNC `GUNC.*.maxCSS_level.tsv` output(s) under a directory.

    Returns dict: genome_name → {pass_gunc, css, contamination_portion,
    n_effective_surplus_clades, taxonomic_level}.
    """
    results: Dict[str, Dict] = {}
    tsvs = sorted(gunc_dir.rglob("GUNC.*.maxCSS_level.tsv"))
    if not tsvs:
        tsvs = sorted(gunc_dir.rglob("*maxCSS_level.tsv"))
    if not tsvs:
        log_warning(f"No GUNC maxCSS_level.tsv found under {gunc_dir}")
        return results

    for tsv in tsvs:
        with open(tsv, newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                name = (row.get('genome') or row.get('Name') or '').strip()
                if not name:
                    continue
                pass_raw = (row.get('pass.GUNC') or '').strip().lower()
                results[name] = {
                    'pass_gunc': pass_raw in ('true', '1', 'yes'),
                    'css': _safe_float(row.get('clade_separation_score')),
                    'contamination_portion': _safe_float(row.get('contamination_portion')),
                    'n_effective_surplus_clades': _safe_float(row.get('n_effective_surplus_clades')),
                    'taxonomic_level': (row.get('taxonomic_level') or '').strip(),
                }
    log_info(f"Parsed GUNC results for {len(results):,} genomes from {len(tsvs)} file(s)")
    return results


def parse_gunc_contig_assignments(
    gunc_dir: Path,
    rank: str = 'phylum',
    min_genes: int = 3,
) -> Dict[str, Dict]:
    """
    Parse GUNC per-contig taxonomy files written by `gunc run
    --contig_taxonomy_output`: one `<genome>.contig_assignments.tsv` per genome
    (columns: contig, tax_level, assignment, count_of_genes_assigned).

    For each genome, at the chosen `rank`, pick each contig's dominant
    assignment (the one with the most genes, requiring at least `min_genes` so a
    one-gene stray does not invent a clade) and count how many DISTINCT clades
    the genome's contigs resolve to. Two or more is a direct per-contig
    taxonomic split — the strongest single chimera signal available here.

    Returns dict: genome_name -> {n_clades, clades (sorted list),
    n_contigs_assigned}. Empty when no files are found.
    """
    results: Dict[str, Dict] = {}
    files = sorted(gunc_dir.rglob("*.contig_assignments.tsv"))
    if not files:
        return results

    for path in files:
        genome = path.name[:-len(".contig_assignments.tsv")]
        # contig -> (best_assignment, best_gene_count) at the requested rank
        best: Dict[str, tuple] = {}
        with open(path, newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                level = (row.get('tax_level') or row.get('taxonomic_level') or '').strip()
                if level != rank:
                    continue
                contig = (row.get('contig') or '').strip()
                assignment = (row.get('assignment') or '').strip()
                if not contig or not assignment:
                    continue
                genes = _safe_float(row.get('count_of_genes_assigned'))
                genes = int(genes) if genes is not None else 0
                if genes < min_genes:
                    continue
                prev = best.get(contig)
                if prev is None or genes > prev[1]:
                    best[contig] = (assignment, genes)

        clades = sorted({a for a, _ in best.values()})
        results[genome] = {
            'n_clades': len(clades),
            'clades': clades,
            'n_contigs_assigned': len(best),
        }

    n_multi = sum(1 for r in results.values() if r['n_clades'] >= 2)
    log_info(
        f"Parsed GUNC per-contig taxonomy for {len(results):,} genomes "
        f"from {len(files)} file(s) at rank '{rank}' "
        f"({n_multi:,} with contigs spanning >=2 clades)"
    )
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
    gunc_fail: Optional[bool] = None,
    gunc_css: Optional[float] = None,
    gunc_contig_clades: Optional[int] = None,
    gc_outlier: Optional[Dict] = None,
) -> Tuple[str, List[str]]:
    """
    Assess chimera risk for a bin/MAG from multiple independent signals.

    Gene-level taxonomy (GUNC) is the strongest signal and drives the score; GC
    heterogeneity is noisy (a single short contig can set the max-min delta) so
    it corroborates rather than drives — on its own GC never reaches Medium.

    Scoring:
      GUNC per-contig taxonomy (needs gunc --contig_taxonomy_output)
        contigs span >= 2 clades → +4   (direct per-contig taxonomic split)
      GUNC gene-level taxonomic inconsistency (aggregate)
        pass.GUNC = False → +3   (genes span multiple clades)
        clade_separation_score > 0.45 → +2  (when not already failed)
      CheckM2 contamination
        > 10 %        → +3    > 5 %         → +2
      GC heterogeneity between contigs (multi-contig bins)
        delta > 10 %  → +2    delta > 5 %   → +1
      Windowed GC variation along a circular contig (weakest: prophage/GC-skew
      look the same as a chimeric junction here)
        delta > 15 %  → +1    (delta > 8 % no longer scores)
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
            score += 2
            reasons.append(f"High inter-contig GC heterogeneity (delta={gc_delta * 100:.1f}%)"
                           + _gc_outlier_note(gc_outlier))
        elif gc_delta > 0.05:
            score += 1
            reasons.append(f"Moderate inter-contig GC heterogeneity (delta={gc_delta * 100:.1f}%)"
                           + _gc_outlier_note(gc_outlier))

    if windowed_gc_delta is not None:
        if windowed_gc_delta > 0.15:
            score += 1
            reasons.append(
                f"Windowed GC variation along circular contig "
                f"(delta={windowed_gc_delta * 100:.1f}%); weak alone - "
                f"prophage/island/GC-skew look the same, corroborate with GUNC"
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

    # GUNC gene-level taxonomic consistency (orthogonal to GC/contamination)
    if gunc_fail is True:
        score += 3
        if gunc_css is not None:
            reasons.append(f"GUNC fail: genes span multiple clades (CSS={gunc_css:.2f})")
        else:
            reasons.append("GUNC fail: genes span multiple clades")
    elif gunc_css is not None and gunc_css > 0.45:
        score += 2
        reasons.append(f"Elevated GUNC clade separation score (CSS={gunc_css:.2f})")

    # GUNC per-contig taxonomy: contigs resolving to >=2 clades is a direct
    # per-contig taxonomic split, the strongest single chimera signal here.
    if gunc_contig_clades is not None and gunc_contig_clades >= 2:
        score += 4
        reasons.append(
            f"GUNC per-contig taxonomy: contigs span {gunc_contig_clades} clades"
        )

    if score >= 5:
        return 'High', reasons
    elif score >= 2:
        return 'Medium', reasons
    elif score >= 1:
        return 'Low', reasons
    else:
        return 'Clean', reasons


def _gc_outlier_note(gc_outlier: Optional[Dict]) -> str:
    """Append which contigs set the inter-contig GC delta, with their lengths,
    so a short stray outlier is recognisable as a weaker reason to worry."""
    if not gc_outlier:
        return ""
    hi_name, hi_gc, hi_len = gc_outlier['high']
    lo_name, lo_gc, lo_len = gc_outlier['low']
    return (f" [high {hi_name} {hi_gc * 100:.1f}% {hi_len:,}bp vs "
            f"low {lo_name} {lo_gc * 100:.1f}% {lo_len:,}bp]")
