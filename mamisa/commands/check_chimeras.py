#!/usr/bin/env python3
"""
MaMISA - check-chimeras command
Detect chimeric MAGs and circular contigs using GC heterogeneity and GTDB-Tk taxonomy signals.

Chimeric MAGs ("puzzle MAGs") are a common artefact in metagenomic binning where contigs
from multiple organisms are incorrectly grouped into a single bin. Circular contigs can
also be chimeric if they were assembled from reads belonging to two different organisms
that share similar coverage and k-mer profiles.

Detection strategy
------------------
1. **GC heterogeneity** (all bins):
   - Bins whose contigs show high GC variance (delta > 10 %) are flagged.

2. **Windowed GC analysis** (single-contig circular candidates ≥ 100 kbp):
   - A sliding window scans the contig for abrupt GC shifts that may indicate
     a chimeric junction (assembly artefact).

3. **GTDB-Tk taxonomy confidence** (optional):
   - Bins with low MSA recovery (< 10 %), classification warnings, or no
     FastANI reference match are flagged as unreliably placed.

4. **CheckM2 contamination** (optional):
   - Contamination > 5 % is a strong independent chimera indicator.
"""

import csv
import sys
import argparse
from pathlib import Path
from typing import Dict, List, Optional

from ..utils.chimera import (
    analyze_bin_gc,
    parse_gtdbtk_summary,
    extract_taxonomy_level,
    has_taxonomy_warning,
    assess_chimera_risk,
)
from ..utils.checkm2 import parse_quality_report
from ..utils.validation import validate_dir_exists
from ..utils.logging import log_info, log_warning, log_error, print_header, print_section


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _find_gtdbtk_summaries(gtdbtk_dir: Path) -> List[Path]:
    """Locate GTDB-Tk summary TSV files under gtdbtk_dir."""
    candidates = [
        gtdbtk_dir / 'gtdbtk.bac120.summary.tsv',
        gtdbtk_dir / 'gtdbtk.ar53.summary.tsv',
        gtdbtk_dir / 'classify' / 'gtdbtk.bac120.summary.tsv',
        gtdbtk_dir / 'classify' / 'gtdbtk.ar53.summary.tsv',
    ]
    found = [p for p in candidates if p.exists()]
    if not found:
        found = sorted(gtdbtk_dir.rglob('gtdbtk.*.summary.tsv'))
    return found


def _load_taxonomy(gtdbtk_dir: Path) -> Dict[str, Dict]:
    summaries = _find_gtdbtk_summaries(gtdbtk_dir)
    if not summaries:
        log_warning(f"No GTDB-Tk summary files found under {gtdbtk_dir}")
        return {}
    taxonomy = {}
    for summary in summaries:
        log_info(f"  Loading taxonomy: {summary.name}")
        taxonomy.update(parse_gtdbtk_summary(summary))
    return taxonomy


def _load_contamination(checkm2_file: Path) -> Dict[str, float]:
    records = parse_quality_report(checkm2_file)
    return {rec['name']: rec['contamination'] for rec in records}


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def run_chimera_check(
    bins_dir: Path,
    output_file: Path,
    extensions: List[str],
    gtdbtk_dir: Optional[Path],
    checkm2_file: Optional[Path],
    gc_window: int,
    gc_step: int,
    taxonomy_level: str,
    dry_run: bool,
) -> Dict[str, int]:
    """
    Run chimera detection for all bins in bins_dir.

    Returns a dict of {risk_level: count}.
    """
    print_section("Loading auxiliary data")

    taxonomy: Dict[str, Dict] = {}
    if gtdbtk_dir:
        taxonomy = _load_taxonomy(gtdbtk_dir)
        if taxonomy:
            log_info(f"Taxonomy loaded for {len(taxonomy):,} genomes")

    contamination_map: Dict[str, float] = {}
    if checkm2_file:
        contamination_map = _load_contamination(checkm2_file)
        log_info(f"CheckM2 contamination loaded for {len(contamination_map):,} genomes")

    print_section("Scanning bin files")

    bin_files: List[Path] = []
    for ext in extensions:
        bin_files.extend(bins_dir.rglob(f"*.{ext}"))
    bin_files = sorted(set(bin_files))

    if not bin_files:
        log_error(f"No bin files found in {bins_dir} with extensions: {extensions}")
        return {}

    log_info(f"Found {len(bin_files):,} bin file(s)")

    print_section("Analysing chimera risk per bin")

    rows: List[Dict] = []
    risk_counts: Dict[str, int] = {'High': 0, 'Medium': 0, 'Low': 0, 'Clean': 0}

    for bin_file in bin_files:
        bin_name = bin_file.stem
        log_info(f"\n  [{bin_name}]")

        # GC analysis
        try:
            gc_data = analyze_bin_gc(bin_file, gc_window, gc_step)
        except Exception as e:
            log_warning(f"  Could not analyse {bin_file.name}: {e}")
            continue

        bin_stats = gc_data['bin_gc_stats']
        windowed_stats = gc_data.get('windowed_gc_stats')

        gc_delta = bin_stats.get('delta') or 0.0
        gc_cv = bin_stats.get('cv') or 0.0
        windowed_delta = windowed_stats.get('delta') if windowed_stats else None

        # Taxonomy
        tax_record = taxonomy.get(bin_name, {})
        classification = tax_record.get('classification', '')
        taxon = extract_taxonomy_level(classification, taxonomy_level) if classification else 'N/A'
        tax_warn = has_taxonomy_warning(tax_record) if tax_record else False

        # Contamination
        contamination = contamination_map.get(bin_name)

        # Risk score
        risk, reasons = assess_chimera_risk(
            gc_delta=gc_delta,
            gc_cv=gc_cv,
            n_contigs=gc_data['n_contigs'],
            contamination=contamination,
            windowed_gc_delta=windowed_delta,
            taxonomy_warning=tax_warn,
        )
        risk_counts[risk] = risk_counts.get(risk, 0) + 1

        log_info(
            f"  Contigs={gc_data['n_contigs']}  "
            f"GC_delta={gc_delta * 100:.1f}%  "
            f"Circular={gc_data['is_circular_candidate']}  "
            f"Contamination={'N/A' if contamination is None else f'{contamination:.1f}%'}  "
            f"Risk={risk}"
        )
        for reason in reasons:
            log_info(f"  ⚠  {reason}")

        rows.append({
            'bin': bin_name,
            'n_contigs': gc_data['n_contigs'],
            'is_circular_candidate': gc_data['is_circular_candidate'],
            'gc_mean_pct': f"{(bin_stats.get('mean') or float('nan')) * 100:.2f}",
            'gc_delta_pct': f"{gc_delta * 100:.2f}",
            'gc_cv_pct': f"{gc_cv * 100:.4f}",
            'windowed_gc_delta_pct': (
                f"{windowed_delta * 100:.2f}" if windowed_delta is not None else 'N/A'
            ),
            'checkm2_contamination': (
                f"{contamination:.2f}" if contamination is not None else 'N/A'
            ),
            'gtdbtk_taxonomy': taxon,
            'gtdbtk_warning': tax_warn,
            'chimera_risk': risk,
            'reasons': ' | '.join(reasons),
        })

    if not dry_run and rows:
        fieldnames = [
            'bin', 'n_contigs', 'is_circular_candidate',
            'gc_mean_pct', 'gc_delta_pct', 'gc_cv_pct', 'windowed_gc_delta_pct',
            'checkm2_contamination', 'gtdbtk_taxonomy', 'gtdbtk_warning',
            'chimera_risk', 'reasons',
        ]
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        log_info(f"\nChimera report written to: {output_file}")

    return risk_counts


# ---------------------------------------------------------------------------
# CLI registration
# ---------------------------------------------------------------------------

def register_parser(subparsers):
    """Register the check-chimeras subcommand."""
    parser = subparsers.add_parser(
        'check-chimeras',
        help='Detect chimeric MAGs and circular contigs using GC and taxonomy signals',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Signal interpretation
---------------------
  GC delta (multi-contig bins)
    > 10 %  High risk — contigs likely from different organisms
    > 5 %   Low/medium risk — worth inspecting
  Windowed GC delta (circular contigs ≥ 100 kbp)
    > 15 %  High risk — possible chimeric junction in circular assembly
    > 8 %   Medium risk
  CheckM2 contamination
    > 10 %  High risk
    > 5 %   Medium risk
  GTDB-Tk warning
    Any warning or MSA percent < 10 % → low-confidence placement

Examples
--------
  # Basic GC-only check
  mamisa check-chimeras \\
    --bins-dir bins/ \\
    --output chimera_report.tsv

  # Full check with GTDB-Tk + CheckM2
  mamisa check-chimeras \\
    --bins-dir bins/ \\
    --output chimera_report.tsv \\
    --gtdbtk-dir gtdbtk_results/ \\
    --checkm2-report checkm2/quality_report.tsv

  # Finer GC resolution for large circular genomes
  mamisa check-chimeras \\
    --bins-dir bins/ \\
    --output chimera_report.tsv \\
    --gc-window 2000 --gc-step 500

  # Dry-run (print to terminal, no file written)
  mamisa check-chimeras --bins-dir bins/ --dry-run
        """
    )

    parser.add_argument('--bins-dir', type=Path, required=True,
                        help='Directory containing bin/MAG FASTA files')
    parser.add_argument('-o', '--output', type=Path, default=Path('chimera_report.tsv'),
                        help='Output TSV report (default: chimera_report.tsv)')

    # Optional integration with other tools
    parser.add_argument('--gtdbtk-dir', type=Path,
                        help='GTDB-Tk output directory (adds taxonomy-based signals)')
    parser.add_argument('--checkm2-report', type=Path,
                        help='CheckM2 quality_report.tsv (adds contamination signal)')

    # GC window parameters
    parser.add_argument('--gc-window', type=int, default=5000,
                        help='Window size for circular-contig GC analysis (default: 5000 bp)')
    parser.add_argument('--gc-step', type=int, default=2500,
                        help='Step size for GC sliding window (default: 2500 bp)')

    parser.add_argument('--taxonomy-level', default='phylum',
                        choices=['domain', 'phylum', 'class', 'order', 'family',
                                 'genus', 'species'],
                        help='Taxonomic level to report in output (default: phylum)')
    parser.add_argument('--extensions', default='fa,fasta,fna',
                        help='Comma-separated bin file extensions (default: fa,fasta,fna)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Print results to terminal without writing output file')

    parser.set_defaults(func=run)
    return parser


def run(args):
    """Execute the check-chimeras command."""
    validate_dir_exists(args.bins_dir, "Bins directory")

    if args.gtdbtk_dir:
        validate_dir_exists(args.gtdbtk_dir, "GTDB-Tk directory")

    if args.checkm2_report and not args.checkm2_report.exists():
        log_error(f"CheckM2 report not found: {args.checkm2_report}")
        sys.exit(1)

    if args.gc_step >= args.gc_window:
        log_warning(
            f"--gc-step ({args.gc_step}) >= --gc-window ({args.gc_window}); "
            "windows will not overlap, which may miss narrow GC shifts"
        )

    extensions = [e.strip() for e in args.extensions.split(',')]

    print_header("MaMISA - Chimera Check")

    risk_counts = run_chimera_check(
        bins_dir=args.bins_dir,
        output_file=args.output,
        extensions=extensions,
        gtdbtk_dir=args.gtdbtk_dir,
        checkm2_file=args.checkm2_report,
        gc_window=args.gc_window,
        gc_step=args.gc_step,
        taxonomy_level=args.taxonomy_level,
        dry_run=args.dry_run,
    )

    if not risk_counts:
        log_error("No bins were analysed — check --bins-dir and --extensions")
        sys.exit(1)

    print_section("CHIMERA RISK SUMMARY")
    total = sum(risk_counts.values())
    for level in ('High', 'Medium', 'Low', 'Clean'):
        count = risk_counts.get(level, 0)
        marker = '⚠ ' if level in ('High', 'Medium') else '  '
        print(f"  {marker}{level:7s}: {count:>6,}")

    flagged = risk_counts.get('High', 0) + risk_counts.get('Medium', 0)
    if total > 0:
        pct = 100 * flagged / total
        print(f"\n  Flagged (High + Medium): {flagged:,} / {total:,} ({pct:.1f}%)")

    if args.dry_run:
        log_warning("\nDRY-RUN mode: no output file written")

    log_info("\n✓ Done!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    register_parser(parser.add_subparsers(dest='command'))
    run(parser.parse_args())
