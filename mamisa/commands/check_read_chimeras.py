#!/usr/bin/env python3
"""
MaMISA - check-read-chimeras command

Detect chimeric contigs and bins by examining the taxonomic identity
of *reads* that map to each contig.

Problem
-------
A chimeric contig is one assembled from reads that originated in two or
more different organisms. This happens when organisms share similar k-mer
profiles, coverage, or when assembler graphs are not cleanly resolved.
The result: a single contig that is literally a mosaic of two genomes.

Standard sequence-composition checks (GC content, tetranucleotide freq.)
can miss this when the organisms have similar GC. The only reliable signal
is the taxonomy of the reads themselves: if reads mapping to one contig
come from two different organisms, the contig is chimeric.

Approach
--------
1. Parse Kraken2 per-read classifications (--kraken2-output)
2. Optionally parse Kraken2 report for taxon name resolution (--kraken2-report)
3. Stream the BAM file to link each read to its contig + position
4. Per-contig: count reads per taxon → dominant fraction, Shannon diversity
5. For long contigs (≥ --window-threshold): windowed analysis to locate
   the chimeric junction — the position where the dominant taxon shifts

Outputs
-------
  chimera_reads_report.tsv       one row per contig, overall stats
  chimera_reads_windows.tsv      one row per window (long contigs only)

Usage
-----
  # Minimal (requires BAM + Kraken2 output)
  mamisa check-read-chimeras \\
    --bam mapping.bam \\
    --kraken2-output kraken2.out \\
    --output-dir chimera_results/

  # With name resolution and window analysis
  mamisa check-read-chimeras \\
    --bam mapping.bam \\
    --kraken2-output kraken2.out \\
    --kraken2-report kraken2.report \\
    --output-dir chimera_results/ \\
    --window 10000 --window-step 5000
"""

import csv
import sys
import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import Counter, defaultdict

from ..utils.read_taxonomy import (
    parse_kraken2_output,
    parse_kraken2_report,
    check_samtools,
    stream_bam,
    ContigReadProfile,
    compute_taxonomy_stats,
    windowed_taxon_profile,
    detect_taxonomic_shifts,
    assess_read_chimera_risk,
    taxon_label,
)
from ..utils.fasta import read_fasta_streaming
from ..utils.validation import validate_file_exists
from ..utils.logging import log_info, log_warning, log_error, print_header, print_section


# ──────────────────────────────────────────────────────────────────────────────
# Core pipeline
# ──────────────────────────────────────────────────────────────────────────────

def build_contig_length_map(assembly: Optional[Path]) -> Dict[str, int]:
    """Return {contig_name: length} from a FASTA file (or empty dict)."""
    if assembly is None or not assembly.exists():
        return {}
    return {name: len(seq) for name, seq in read_fasta_streaming(assembly)}


def aggregate_reads(
    bam_file: Path,
    read_taxid: Dict[str, int],
    contig_lengths: Dict[str, int],
    min_mapq: int,
    window_threshold: int,
    min_reads: int,
) -> Dict[str, ContigReadProfile]:
    """
    Stream the BAM once and fill one ContigReadProfile per contig.

    Reads with no Kraken2 classification get taxid = -1 (not in Kraken2 output),
    distinct from taxid = 0 (explicitly unclassified by Kraken2).
    """
    profiles: Dict[str, ContigReadProfile] = {}
    n_reads_total = 0
    n_reads_no_kraken = 0

    log_info("Streaming BAM file…")

    for read_id, contig, pos in stream_bam(bam_file, min_mapq=min_mapq):
        n_reads_total += 1

        # Look up Kraken2 classification (-1 means "not in Kraken2 output")
        taxid = read_taxid.get(read_id, -1)
        if taxid == -1:
            n_reads_no_kraken += 1

        if contig not in profiles:
            length = contig_lengths.get(contig, 0)
            profiles[contig] = ContigReadProfile(
                name=contig,
                length=length,
                window_threshold=window_threshold,
            )

        profiles[contig].add_read(pos, taxid)

        if n_reads_total % 1_000_000 == 0:
            log_info(f"  Processed {n_reads_total:,} reads, "
                     f"{len(profiles):,} contigs so far…")

    log_info(f"Finished: {n_reads_total:,} reads → {len(profiles):,} contigs")
    if n_reads_no_kraken > 0:
        pct = 100 * n_reads_no_kraken / n_reads_total if n_reads_total else 0
        log_warning(
            f"{n_reads_no_kraken:,} reads ({pct:.1f}%) not found in Kraken2 output "
            "(counted as taxid=-1 = 'not classified by Kraken2')"
        )

    return profiles


def analyze_profiles(
    profiles: Dict[str, ContigReadProfile],
    taxid_info: Dict,
    window: int,
    window_step: int,
    min_reads: int,
    min_reads_per_window: int,
    exclude_unclassified: bool,
) -> Tuple[List[Dict], List[Dict]]:
    """
    Compute per-contig stats and windowed analysis.

    Returns (contig_rows, window_rows).
    """
    contig_rows: List[Dict] = []
    window_rows: List[Dict] = []

    for contig_name, profile in sorted(profiles.items()):
        if profile.n_reads < min_reads:
            continue

        stats = compute_taxonomy_stats(
            profile.taxid_counts, taxid_info,
            exclude_unclassified=exclude_unclassified,
        )

        # Windowed analysis for long contigs
        windows = []
        n_shifts = 0
        if profile.positions:
            windows = windowed_taxon_profile(
                profile.positions,
                contig_length=profile.length,
                window=window,
                step=window_step,
                taxid_info=taxid_info,
                min_reads_per_window=min_reads_per_window,
            )
            shifts = detect_taxonomic_shifts(windows)
            n_shifts = len(shifts)

        risk, reasons = assess_read_chimera_risk(
            dominant_fraction=stats['dominant_fraction'],
            n_taxa=stats['n_taxa'],
            n_reads=profile.n_reads,
            n_taxon_shifts=n_shifts,
            classified_fraction=stats['classified_fraction'],
        )

        contig_rows.append({
            'contig': contig_name,
            'contig_length': profile.length if profile.length else 'N/A',
            'n_reads_mapped': stats['n_reads_total'],
            'n_classified': stats['n_classified'],
            'classified_frac': f"{stats['classified_fraction'] * 100:.1f}",
            'n_taxa': stats['n_taxa'],
            'dominant_taxon': stats['dominant_name'],
            'dominant_frac': f"{stats['dominant_fraction'] * 100:.1f}",
            'shannon_diversity': f"{stats['shannon_diversity']:.4f}",
            'windowed_analysis': len(windows) > 0,
            'n_taxon_shifts': n_shifts,
            'chimera_risk': risk,
            'reasons': ' | '.join(reasons),
        })

        # Emit one row per window for long contigs
        for w in windows:
            window_rows.append({
                'contig': contig_name,
                'window_start': w['start'],
                'window_end': w['end'],
                'n_reads': w['n_reads'],
                'dominant_taxon': w['dominant_name'],
                'dominant_frac': f"{w['dominant_fraction'] * 100:.1f}",
                'n_taxa': w['n_taxa'],
                'is_taxon_shift': w['is_shift'],
            })

    return contig_rows, window_rows


# ──────────────────────────────────────────────────────────────────────────────
# Output helpers
# ──────────────────────────────────────────────────────────────────────────────

CONTIG_FIELDS = [
    'contig', 'contig_length', 'n_reads_mapped', 'n_classified',
    'classified_frac', 'n_taxa', 'dominant_taxon', 'dominant_frac',
    'shannon_diversity', 'windowed_analysis', 'n_taxon_shifts',
    'chimera_risk', 'reasons',
]

WINDOW_FIELDS = [
    'contig', 'window_start', 'window_end', 'n_reads',
    'dominant_taxon', 'dominant_frac', 'n_taxa', 'is_taxon_shift',
]


def write_tsv(rows: List[Dict], fields: List[str], path: Path):
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t',
                                extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    log_info(f"Written: {path} ({len(rows):,} rows)")


def print_summary(contig_rows: List[Dict]):
    risk_counts: Dict[str, int] = defaultdict(int)
    for row in contig_rows:
        risk_counts[row['chimera_risk']] += 1

    print_section("CHIMERA RISK SUMMARY (read-level)")
    total = sum(risk_counts.values())
    for level in ('High', 'Medium', 'Low', 'Clean', 'Insufficient'):
        count = risk_counts.get(level, 0)
        if count == 0:
            continue
        marker = '⚠ ' if level in ('High', 'Medium') else '  '
        print(f"  {marker}{level:12s}: {count:>6,}")

    flagged = risk_counts.get('High', 0) + risk_counts.get('Medium', 0)
    if total > 0:
        pct = 100 * flagged / total
        print(f"\n  Flagged (High + Medium): {flagged:,} / {total:,} ({pct:.1f}%)")

    # Show top 10 highest-risk contigs
    high_risk = [r for r in contig_rows if r['chimera_risk'] in ('High', 'Medium')]
    if high_risk:
        print_section("TOP FLAGGED CONTIGS")
        for row in sorted(
            high_risk,
            key=lambda r: (r['chimera_risk'] != 'High', float(row['dominant_frac']))
        )[:10]:
            print(
                f"  {row['contig'][:50]:<50s}  "
                f"risk={row['chimera_risk']:<8s}  "
                f"dominant={row['dominant_frac']}%  "
                f"taxa={row['n_taxa']}  "
                f"shifts={row['n_taxon_shifts']}"
            )


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def register_parser(subparsers):
    parser = subparsers.add_parser(
        'check-read-chimeras',
        help='Detect chimeric contigs by read-level taxonomy (Kraken2 + BAM)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
How it works
------------
Each read that maps to a contig carries a Kraken2 taxonomic assignment.
If a contig is genuine, nearly all its reads should belong to one organism.
If reads from two organisms mapped to the same contig, it is chimeric.

For long contigs (≥ --window-threshold), a sliding window tracks where
the dominant taxon changes — pinpointing the chimeric junction.

Required inputs
---------------
  --bam              BAM file (reads mapped to assembly; must be sorted)
  --kraken2-output   Kraken2 per-read output (kraken2 --output FILE)

Recommended
-----------
  --kraken2-report   Kraken2 summary report (kraken2 --report FILE)
                     Enables human-readable taxon names in output.
  --assembly         FASTA assembly (needed for accurate window boundaries
                     on circular/long contigs)

Examples
--------
  # Basic check
  mamisa check-read-chimeras \\
    --bam mapping.bam \\
    --kraken2-output kraken2.out \\
    --output-dir chimera_results/

  # Full check with name resolution and assembly
  mamisa check-read-chimeras \\
    --bam mapping.bam \\
    --kraken2-output kraken2.out \\
    --kraken2-report kraken2.report \\
    --assembly assembly.fa \\
    --output-dir chimera_results/ \\
    --window 10000 --window-step 2000

  # Bin-level check: point --bam at reads mapped to your bins FASTA
  # and use the binned FASTA as --assembly

  # Only report contigs with ≥ 50 mapped reads
  mamisa check-read-chimeras \\
    --bam mapping.bam \\
    --kraken2-output kraken2.out \\
    --output-dir chimera_results/ \\
    --min-reads 50
        """
    )

    # Required
    parser.add_argument('--bam', type=Path, required=True,
                        help='BAM file (reads mapped to assembly, sorted)')
    parser.add_argument('--kraken2-output', type=Path, required=True,
                        help='Kraken2 per-read output file (kraken2 --output)')
    parser.add_argument('-o', '--output-dir', type=Path, required=True,
                        help='Output directory for TSV reports')

    # Optional enhancements
    parser.add_argument('--kraken2-report', type=Path,
                        help='Kraken2 report file (kraken2 --report) for taxon name lookup')
    parser.add_argument('--assembly', type=Path,
                        help='Assembly FASTA (provides contig lengths for windowed analysis)')

    # Filtering
    parser.add_argument('--min-mapq', type=int, default=20,
                        help='Minimum mapping quality to include a read (default: 20)')
    parser.add_argument('--min-reads', type=int, default=10,
                        help='Minimum mapped reads per contig to analyse (default: 10)')
    parser.add_argument('--min-reads-per-window', type=int, default=5,
                        help='Minimum reads in a window to include it (default: 5)')

    # Windowed analysis
    parser.add_argument('--window-threshold', type=int, default=50_000,
                        help='Contig length above which windowed analysis is run (default: 50000 bp)')
    parser.add_argument('--window', type=int, default=10_000,
                        help='Window size for long-contig taxon profiling (default: 10000 bp)')
    parser.add_argument('--window-step', type=int, default=5_000,
                        help='Step between windows (default: 5000 bp)')

    # Taxonomy
    parser.add_argument('--exclude-unclassified', action='store_true',
                        help='Exclude unclassified reads (taxid=0) when computing '
                             'dominant taxon fraction (default: include them)')

    parser.add_argument('--dry-run', action='store_true',
                        help='Run analysis but do not write output files')

    parser.set_defaults(func=run)
    return parser


def run(args):
    """Execute check-read-chimeras."""

    validate_file_exists(args.bam, "BAM file")
    validate_file_exists(args.kraken2_output, "Kraken2 output file")

    if args.kraken2_report and not args.kraken2_report.exists():
        log_error(f"Kraken2 report not found: {args.kraken2_report}")
        sys.exit(1)

    if args.assembly and not args.assembly.exists():
        log_error(f"Assembly FASTA not found: {args.assembly}")
        sys.exit(1)

    if not check_samtools():
        sys.exit(1)

    if not args.dry_run:
        args.output_dir.mkdir(parents=True, exist_ok=True)

    print_header("MaMISA - Chimera Check (Read-Level Taxonomy)")

    # ── Step 1: Load Kraken2 data ─────────────────────────────────────────
    print_section("STEP 1: Loading Kraken2 classifications")
    read_taxid = parse_kraken2_output(args.kraken2_output)

    taxid_info: Dict = {}
    if args.kraken2_report:
        taxid_info = parse_kraken2_report(args.kraken2_report)
    else:
        log_warning(
            "No --kraken2-report provided — taxon IDs will be reported "
            "as 'taxid:NNN' without names"
        )

    # ── Step 2: Contig lengths ────────────────────────────────────────────
    print_section("STEP 2: Loading contig lengths")
    if args.assembly:
        contig_lengths = build_contig_length_map(args.assembly)
        log_info(f"Loaded lengths for {len(contig_lengths):,} contigs")
    else:
        contig_lengths = {}
        log_warning(
            "No --assembly provided — contig lengths unknown; "
            "windowed analysis will be disabled"
        )

    # ── Step 3: Stream BAM and aggregate ─────────────────────────────────
    print_section("STEP 3: Streaming BAM and aggregating reads per contig")
    profiles = aggregate_reads(
        bam_file=args.bam,
        read_taxid=read_taxid,
        contig_lengths=contig_lengths,
        min_mapq=args.min_mapq,
        window_threshold=args.window_threshold,
        min_reads=args.min_reads,
    )

    # ── Step 4: Analyse each contig ───────────────────────────────────────
    print_section("STEP 4: Analysing read taxonomy per contig")
    contig_rows, window_rows = analyze_profiles(
        profiles=profiles,
        taxid_info=taxid_info,
        window=args.window,
        window_step=args.window_step,
        min_reads=args.min_reads,
        min_reads_per_window=args.min_reads_per_window,
        exclude_unclassified=args.exclude_unclassified,
    )

    log_info(f"Analysed {len(contig_rows):,} contigs "
             f"({len(window_rows):,} window records for long contigs)")

    # ── Step 5: Write output ──────────────────────────────────────────────
    if not args.dry_run:
        print_section("STEP 5: Writing output")
        report_path = args.output_dir / "chimera_read_report.tsv"
        write_tsv(contig_rows, CONTIG_FIELDS, report_path)

        if window_rows:
            windows_path = args.output_dir / "chimera_read_windows.tsv"
            write_tsv(window_rows, WINDOW_FIELDS, windows_path)
    else:
        log_warning("DRY-RUN: output files not written")

    # ── Summary ───────────────────────────────────────────────────────────
    print_summary(contig_rows)

    log_info("\n✓ Done!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    register_parser(parser.add_subparsers(dest='command'))
    run(parser.parse_args())
