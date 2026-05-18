#!/usr/bin/env python3
"""
MaMISA - classify-clipping command

Classify each anvi'o soft-clipping position as one of:
  end_artefact      — within 500 bp of a contig terminus; assembly edge effect
  repeat_collapse   — depth spike and/or low-entropy clipped bases; collapsed repeat
  deletion_artefact — local depth drop; internal deletion or coverage collapse
  chimera_candidate — discordant pairs + large inserts + optional taxonomy shift
  sv_candidate      — structural-variant signal without depth anomaly or taxon shift
  low_confidence    — insufficient or ambiguous evidence

Evidence signals used
---------------------
  depth_ratio          local_depth / contig_mean_depth
                        spike → repeat; drop → deletion
  discordant_fraction  primary reads not in a proper pair
                        elevated → inter-chromosomal or large-scale rearrangement
  large_insert_fraction primary proper-pair reads whose |TLEN| > mean + 3 σ
                        elevated → large insertion / structural variant
  strand_fwd_fraction  fraction of primary reads on the forward strand
                        extreme bias → strand-specific clipping artefact
  clipped_base_entropy Shannon entropy (bits) of the soft-clipped sequences
                        low → repetitive bases (tandem repeat boundary)
                        high → diverse sequence (genuine break or chimera)
  taxonomy_shift       True when a check-read-chimeras window overlapping this
                        position is flagged as a taxonomic shift

BAM categories
--------------
  Depth            primary + secondary (-F 2048, supplementary excluded)
  Pair statistics  primary only (FLAG & 0x100 == 0)
  Supplementary    excluded throughout
"""

import csv
import math
import argparse
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from ..utils.misassembly import load_clipping_data
from ..utils.bam_stats import (
    check_samtools,
    get_contig_mean_depths,
    estimate_insert_size,
    collect_position_stats,
)
from ..utils.validation import validate_file_exists, validate_dir_exists
from ..utils.logging import log_info, log_warning, log_error, print_header, print_section


# ──────────────────────────────────────────────────────────────────────────────
# Contig lengths from BAM header
# ──────────────────────────────────────────────────────────────────────────────

def get_contig_lengths_from_bam(bam_file: Path) -> Dict[str, int]:
    """Parse @SQ lines from BAM header to get contig lengths."""
    cmd = ['samtools', 'view', '-H', str(bam_file)]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        log_warning(f"Could not read BAM header: {e.stderr[:200]}")
        return {}

    lengths: Dict[str, int] = {}
    for line in result.stdout.splitlines():
        if not line.startswith('@SQ'):
            continue
        parts = {kv.split(':', 1)[0]: kv.split(':', 1)[1]
                 for kv in line.split('\t')[1:]
                 if ':' in kv}
        name = parts.get('SN', '')
        try:
            length = int(parts.get('LN', 0))
        except ValueError:
            length = 0
        if name and length:
            lengths[name] = length

    log_info(f"BAM header: {len(lengths):,} contigs with lengths")
    return lengths


# ──────────────────────────────────────────────────────────────────────────────
# Taxonomy shift cross-reference
# ──────────────────────────────────────────────────────────────────────────────

def load_taxonomy_shifts(
    windows_tsv: Path,
    flank: int = 5000,
) -> Dict[str, List[Tuple[int, int]]]:
    """
    Parse the chimera_read_windows.tsv from check-read-chimeras and return
    a dict mapping contig → list of (start, end) tuples for shifted windows.

    `flank` bp is added on each side so that a clipping position just outside
    a window boundary is still flagged.
    """
    shifts: Dict[str, List[Tuple[int, int]]] = {}

    with open(windows_tsv, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if (row.get('is_shift', '0') or '0').strip() not in ('1', 'True', 'true'):
                continue
            contig = (row.get('contig') or '').strip()
            try:
                start = max(0, int(row['window_start']) - flank)
                end = int(row['window_end']) + flank
            except (KeyError, ValueError):
                continue
            if contig:
                shifts.setdefault(contig, []).append((start, end))

    n_contigs = len(shifts)
    n_windows = sum(len(v) for v in shifts.values())
    log_info(f"Taxonomy shifts loaded: {n_windows:,} windows across {n_contigs:,} contigs")
    return shifts


def has_taxonomy_shift(
    contig: str,
    clip_pos: int,
    shifts: Dict[str, List[Tuple[int, int]]],
) -> bool:
    windows = shifts.get(contig, [])
    return any(start <= clip_pos <= end for start, end in windows)


# ──────────────────────────────────────────────────────────────────────────────
# Classification logic
# ──────────────────────────────────────────────────────────────────────────────

def _to_float(v) -> float:
    if v in ('N/A', None, ''):
        return float('nan')
    try:
        return float(v)
    except (TypeError, ValueError):
        return float('nan')


def classify_position(
    stats: dict,
    taxonomy_shift: bool = False,
) -> Tuple[str, str, List[str]]:
    """
    Classify one clipping position given its BAM-derived stats.

    Returns (classification, confidence, evidence_list).

    Classification labels
    ---------------------
    end_artefact      near contig terminus (<500 bp); high-confidence noise
    repeat_collapse   depth spike and/or low-entropy clipped bases
    deletion_artefact local depth drop; internal deletion or local collapse
    chimera_candidate discordant pairs + large inserts + optional taxon shift
    sv_candidate      SV signal (discordant) without depth anomaly or taxon shift
    low_confidence    ambiguous or insufficient evidence

    Confidence: High / Medium / Low
    """
    evidence: List[str] = []

    near_end = bool(stats.get('near_contig_end', False))
    depth_ratio = _to_float(stats.get('depth_ratio'))
    disc_frac = _to_float(stats.get('discordant_fraction'))
    large_ins = _to_float(stats.get('large_insert_fraction'))
    strand_fwd = _to_float(stats.get('strand_fwd_fraction'))
    entropy = _to_float(stats.get('clipped_base_entropy'))
    n_primary = int(stats.get('primary_reads_in_window', 0))

    # ── Rule 1: near contig end → artefact ─────────────────────────────────
    if near_end:
        evidence.append('near_contig_end')
        return 'end_artefact', 'High', evidence

    # ── Insufficient depth: can't say anything meaningful ──────────────────
    if n_primary < 5:
        evidence.append('low_read_count')
        return 'low_confidence', 'Low', evidence

    # ── Score each category ─────────────────────────────────────────────────
    repeat_score = 0
    deletion_score = 0
    chimera_score = 0
    sv_score = 0

    # Depth ratio
    if not math.isnan(depth_ratio):
        if depth_ratio > 3.0:
            repeat_score += 3
            evidence.append(f'depth_spike_{depth_ratio:.1f}x')
        elif depth_ratio > 2.0:
            repeat_score += 2
            evidence.append(f'depth_elevation_{depth_ratio:.1f}x')
        elif depth_ratio < 0.25:
            deletion_score += 3
            evidence.append(f'depth_drop_{depth_ratio:.2f}x')
        elif depth_ratio < 0.5:
            deletion_score += 2
            evidence.append(f'depth_low_{depth_ratio:.2f}x')

    # Discordant pairs
    if not math.isnan(disc_frac):
        if disc_frac > 0.5:
            chimera_score += 3
            sv_score += 2
            evidence.append(f'high_discordant_{disc_frac:.2f}')
        elif disc_frac > 0.3:
            chimera_score += 2
            sv_score += 2
            evidence.append(f'discordant_{disc_frac:.2f}')
        elif disc_frac > 0.15:
            sv_score += 1
            evidence.append(f'mild_discordant_{disc_frac:.2f}')

    # Large insert fraction
    if not math.isnan(large_ins):
        if large_ins > 0.3:
            chimera_score += 2
            sv_score += 2
            evidence.append(f'high_large_insert_{large_ins:.2f}')
        elif large_ins > 0.15:
            chimera_score += 1
            sv_score += 1
            evidence.append(f'large_insert_{large_ins:.2f}')

    # Clipped-base entropy (low = repetitive = repeat boundary)
    if not math.isnan(entropy):
        if entropy < 0.8:
            repeat_score += 2
            evidence.append(f'low_clip_entropy_{entropy:.2f}bits')
        elif entropy < 1.5:
            repeat_score += 1
            evidence.append(f'moderate_clip_entropy_{entropy:.2f}bits')

    # Extreme strand bias
    if not math.isnan(strand_fwd):
        if strand_fwd < 0.05 or strand_fwd > 0.95:
            chimera_score += 1
            evidence.append(f'extreme_strand_bias_{strand_fwd:.2f}')
        elif strand_fwd < 0.1 or strand_fwd > 0.9:
            evidence.append(f'strand_imbalance_{strand_fwd:.2f}')

    # Taxonomy shift (strong signal for chimera)
    if taxonomy_shift:
        chimera_score += 3
        evidence.append('taxonomy_shift')

    # ── Decision ────────────────────────────────────────────────────────────
    scores = {
        'repeat_collapse': repeat_score,
        'deletion_artefact': deletion_score,
        'chimera_candidate': chimera_score,
        'sv_candidate': sv_score,
    }

    max_score = max(scores.values())

    if max_score < 2:
        return 'low_confidence', 'Low', evidence

    # Break ties using biological priority:
    # chimera > repeat > deletion > sv  (most actionable first)
    priority = ['chimera_candidate', 'repeat_collapse', 'deletion_artefact', 'sv_candidate']
    winner = max(priority, key=lambda c: scores[c])

    if max_score >= 5:
        confidence = 'High'
    elif max_score >= 3:
        confidence = 'Medium'
    else:
        confidence = 'Low'

    return winner, confidence, evidence


# ──────────────────────────────────────────────────────────────────────────────
# Output writing
# ──────────────────────────────────────────────────────────────────────────────

_OUTPUT_FIELDS = [
    'contig', 'clip_pos', 'contig_length', 'local_depth',
    'primary_reads_in_window', 'contig_mean_depth', 'depth_ratio',
    'discordant_fraction', 'large_insert_fraction', 'strand_fwd_fraction',
    'clipped_base_entropy', 'near_contig_end', 'taxonomy_shift',
    'classification', 'confidence', 'evidence',
]


def write_report(
    results: Dict[Tuple[str, int], dict],
    taxonomy_shifts: Dict[str, List[Tuple[int, int]]],
    output_path: Path,
) -> dict:
    """Classify all positions and write the report TSV. Returns label counts."""
    label_counts: Dict[str, int] = {}

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=_OUTPUT_FIELDS, delimiter='\t')
        writer.writeheader()

        for (contig, clip_pos), stats in sorted(results.items()):
            tax_shift = has_taxonomy_shift(contig, clip_pos, taxonomy_shifts)
            label, confidence, evidence = classify_position(stats, tax_shift)
            label_counts[label] = label_counts.get(label, 0) + 1

            writer.writerow({
                **stats,
                'taxonomy_shift': tax_shift,
                'classification': label,
                'confidence': confidence,
                'evidence': '|'.join(evidence) if evidence else 'none',
            })

    return label_counts


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────

def run_classify_clipping(
    bam_file: Path,
    misassembly_dir: Path,
    output: Path,
    min_mapq: int,
    window: int,
    min_clip_coverage: int,
    taxonomy_windows: Optional[Path],
    taxonomy_flank: int,
    dry_run: bool,
) -> None:

    print_header("MaMISA - Classify Clipping Positions")

    if not check_samtools():
        sys.exit(1)

    # ── 1. Load clipping data ─────────────────────────────────────────────
    print_section("Loading clipping data")
    clipping_data = load_clipping_data(misassembly_dir, min_coverage=min_clip_coverage)
    if not clipping_data:
        log_error("No clipping positions found — nothing to classify")
        sys.exit(1)

    n_positions = sum(len(v) for v in clipping_data.values())
    log_info(f"Loaded {n_positions:,} clipping positions across {len(clipping_data):,} contigs")

    # ── 2. Contig lengths (BAM header) ────────────────────────────────────
    print_section("Reading contig metadata from BAM")
    contig_lengths = get_contig_lengths_from_bam(bam_file)

    # ── 3. Contig mean depths ─────────────────────────────────────────────
    print_section("Computing contig mean depths")
    contig_mean_depths = get_contig_mean_depths(bam_file, min_mapq)
    if not contig_mean_depths:
        log_warning("samtools coverage returned no data — depth_ratio will be N/A")

    # ── 4. Insert-size estimate ───────────────────────────────────────────
    print_section("Estimating library insert size")
    insert_mean, insert_std = estimate_insert_size(bam_file, min_mapq=min_mapq)
    if insert_mean == 0:
        log_warning("Insert-size estimation failed (single-end or unstranded library?) — "
                    "large_insert_fraction will be N/A")

    # ── 5. Per-position BAM stats ─────────────────────────────────────────
    if dry_run:
        log_info("\nDRY-RUN: skipping BAM streaming")
        log_info(f"Would process {n_positions:,} positions with window={window} bp, "
                 f"min_mapq={min_mapq}")
        log_info(f"Would write report to: {output}")
        return

    print_section("Streaming BAM for per-position evidence")
    position_stats = collect_position_stats(
        bam_file=bam_file,
        clipping_data=clipping_data,
        window=window,
        min_mapq=min_mapq,
        contig_mean_depths=contig_mean_depths,
        contig_lengths=contig_lengths,
        insert_mean=insert_mean,
        insert_std=insert_std,
    )

    if not position_stats:
        log_error("No position stats collected — check BAM file and clipping data")
        sys.exit(1)

    # ── 6. Taxonomy shift cross-reference ─────────────────────────────────
    taxonomy_shifts: Dict[str, List[Tuple[int, int]]] = {}
    if taxonomy_windows:
        print_section("Loading taxonomy shift windows")
        taxonomy_shifts = load_taxonomy_shifts(taxonomy_windows, flank=taxonomy_flank)

    # ── 7. Classify and write ─────────────────────────────────────────────
    print_section("Classifying positions")
    output.parent.mkdir(parents=True, exist_ok=True)
    label_counts = write_report(position_stats, taxonomy_shifts, output)

    # ── 8. Summary ─────────────────────────────────────────────────────────
    print_section("CLASSIFICATION SUMMARY")
    total = sum(label_counts.values())

    label_order = [
        'end_artefact',
        'repeat_collapse',
        'deletion_artefact',
        'chimera_candidate',
        'sv_candidate',
        'low_confidence',
    ]
    label_descriptions = {
        'end_artefact':      'Near contig end (assembly artefact)',
        'repeat_collapse':   'Repeat collapse / tandem repeat boundary',
        'deletion_artefact': 'Local deletion or coverage collapse',
        'chimera_candidate': 'Chimera candidate (discordant + insert anomaly)',
        'sv_candidate':      'Structural variant signal',
        'low_confidence':    'Low confidence / insufficient evidence',
    }

    for label in label_order:
        n = label_counts.get(label, 0)
        if n > 0:
            pct = 100 * n / total if total else 0
            desc = label_descriptions.get(label, label)
            print(f"  {label:<22s}  {n:>6,}  ({pct:5.1f}%)  {desc}")

    print(f"\n  Total positions classified: {total:,}")
    print(f"  Report written to: {output}")


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def register_parser(subparsers):
    parser = subparsers.add_parser(
        'classify-clipping',
        help='Classify each clipping position using BAM evidence',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Classification labels
---------------------
  end_artefact      Near contig terminus (<500 bp); likely assembly edge noise.
  repeat_collapse   Coverage spike and/or low-entropy clipped bases.
                    Boundary of a collapsed tandem or low-complexity repeat.
  deletion_artefact Local depth drop. Internal deletion, local assembly collapse,
                    or a genuine genomic deletion relative to the reference.
  chimera_candidate Discordant pairs + large insert fraction, especially when
                    paired with a taxonomy shift in the chimera windows file.
                    Prioritise for manual review or splitting.
  sv_candidate      Discordant/large-insert signal without depth anomaly or
                    taxon shift. Possible inversion, translocation, or mobile
                    element insertion.
  low_confidence    Evidence below threshold or read count too low to classify.

Taxonomy cross-reference
------------------------
  --taxonomy-windows  chimera_read_windows.tsv from check-read-chimeras.
                      Any position within --taxonomy-flank bp of a shifted
                      window gets taxonomy_shift=True and contributes +3 to
                      the chimera_candidate score.

Examples
--------
  # Basic: BAM evidence only
  mamisa classify-clipping \\
    --bam mapping.bam \\
    --misassemblies misasm_dir/ \\
    --output clipping_classified.tsv

  # With taxonomy cross-reference
  mamisa classify-clipping \\
    --bam mapping.bam \\
    --misassemblies misasm_dir/ \\
    --taxonomy-windows chimera_dir/chimera_read_windows.tsv \\
    --output clipping_classified.tsv

  # Selective: only classify high-coverage clipping sites
  mamisa classify-clipping \\
    --bam mapping.bam \\
    --misassemblies misasm_dir/ \\
    --min-clip-coverage 50 \\
    --output clipping_classified.tsv
        """
    )

    parser.add_argument('--bam', type=Path, required=True,
                        help='Sorted, indexed BAM file (reads mapped to assembly)')
    parser.add_argument('-m', '--misassemblies', type=Path, required=True,
                        help='Directory containing *-clipping.txt files from anvi\'o')
    parser.add_argument('-o', '--output', type=Path, required=True,
                        help='Output classification TSV file')

    parser.add_argument('--min-mapq', type=int, default=20,
                        help='Minimum mapping quality (default: 20)')
    parser.add_argument('--window', type=int, default=500,
                        help='Read window around each clipping position in bp (default: 500)')
    parser.add_argument('--min-clip-coverage', type=int, default=0, metavar='N',
                        help='Only classify positions with ≥N clipped reads (default: 0 = all)')

    parser.add_argument('--taxonomy-windows', type=Path, metavar='TSV',
                        help='chimera_read_windows.tsv from check-read-chimeras '
                             '(enables taxonomy_shift cross-reference)')
    parser.add_argument('--taxonomy-flank', type=int, default=5000, metavar='BP',
                        help='Extend taxonomy shift windows by this many bp on each side '
                             '(default: 5000)')

    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be done without processing the BAM')

    parser.set_defaults(func=run)
    return parser


def run(args):
    validate_file_exists(args.bam, "BAM file")
    validate_dir_exists(args.misassemblies, "Misassemblies directory")

    if args.taxonomy_windows and not args.taxonomy_windows.exists():
        log_error(f"Taxonomy windows file not found: {args.taxonomy_windows}")
        sys.exit(1)

    run_classify_clipping(
        bam_file=args.bam,
        misassembly_dir=args.misassemblies,
        output=args.output,
        min_mapq=args.min_mapq,
        window=args.window,
        min_clip_coverage=args.min_clip_coverage,
        taxonomy_windows=args.taxonomy_windows,
        taxonomy_flank=args.taxonomy_flank,
        dry_run=args.dry_run,
    )
