#!/usr/bin/env python3
"""
MaMISA - filter-misassemblies command
Integrated assembly filtering with HQ genome awareness and misassembly handling.

Clipping position selection
---------------------------
anvi-script-find-misassemblies reports every position where soft-clipping
exceeds its own internal threshold (e.g. 100x coverage).  Two modes:

  --safe (default)
      Use ALL positions reported by anvi'o.  Conservative: you trust the
      tool's threshold completely.

  --min-clip-coverage N
      Only split at positions where ≥N reads are soft-clipped, providing a
      second, higher-confidence filter on top of anvi'o's own threshold.
      Higher N → fewer splits → less aggressive.
      Mutually exclusive with --safe.

HQ circular contigs
-------------------
Large, high-quality, single-sequence circular contigs (complete genomes)
that have a clipping zone cannot simply be removed or kept intact — they
need to be split at the clipping position so the two genomic pieces can
be used downstream.

  --split-hq-circular
      Force-split any HQ contig that (a) has at least one clipping position
      AND (b) is a circular candidate.  A contig is treated as circular when
      its length ≥ --hq-circular-min-length OR its FASTA header contains a
      recognised circularity annotation (circular=true, circular=Y, etc.).

Chimera awareness
-----------------
  --chimera-report PATH
      TSV produced by check-chimeras or check-read-chimeras.  Any contig
      listed as High or Medium risk (--chimera-risk-threshold) is treated as
      hq_problematic even if it has no clipping zone, preventing it from
      being silently extracted as a clean HQ genome.
"""

import re
import sys
import csv
import argparse
from pathlib import Path
from typing import Set, Dict, List, Tuple, Optional

from ..utils.fasta import read_fasta_streaming, write_fasta, extract_contig_id
from ..utils.misassembly import load_clipping_positions, load_clipping_data
from ..utils.validation import validate_file_exists, validate_dir_exists
from ..utils.logging import log_info, log_error, log_warning, print_header, print_section


# ──────────────────────────────────────────────────────────────────────────────
# Chimera report loading
# ──────────────────────────────────────────────────────────────────────────────

_RISK_ORDER = {'High': 3, 'Medium': 2, 'Low': 1, 'Clean': 0, 'Insufficient': 0}


def load_chimera_flagged(
    chimera_report: Path,
    risk_threshold: str = 'Medium',
) -> Set[str]:
    """
    Load contig/bin names whose chimera_risk is >= risk_threshold.

    Accepts output from both check-chimeras (column 'bin') and
    check-read-chimeras (column 'contig').
    """
    flagged: Set[str] = set()
    min_score = _RISK_ORDER.get(risk_threshold, 2)

    with open(chimera_report, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        name_col = 'contig' if 'contig' in (reader.fieldnames or []) else 'bin'
        for row in reader:
            name = (row.get(name_col) or '').strip()
            risk = (row.get('chimera_risk') or '').strip()
            if name and _RISK_ORDER.get(risk, 0) >= min_score:
                flagged.add(name)

    log_info(
        f"Chimera report: {len(flagged):,} contigs flagged "
        f"at {risk_threshold}+ risk from {chimera_report.name}"
    )
    return flagged


# ──────────────────────────────────────────────────────────────────────────────
# Circular contig detection
# ──────────────────────────────────────────────────────────────────────────────

# Patterns written by Flye (circular=Y), Unicycler (circular=true), etc.
_CIRCULAR_RE = re.compile(
    r'\bcircular[=:]\s*(true|yes|y|1)\b',
    re.IGNORECASE,
)


def is_circular_candidate(
    header: str,
    seq_len: int,
    min_length: int,
) -> bool:
    """
    Return True if a contig looks like a complete circular genome.

    Criteria (either is sufficient):
      1. FASTA header contains a circularity annotation (Flye, Unicycler, …)
      2. Sequence length >= min_length (length-based heuristic)
    """
    if _CIRCULAR_RE.search(header):
        return True
    return seq_len >= min_length


# ──────────────────────────────────────────────────────────────────────────────
# HQ ID loading
# ──────────────────────────────────────────────────────────────────────────────

def get_hq_contig_ids(hq_dir: Path) -> Set[str]:
    """
    Extract contig IDs from HQ genome filenames in hq_dir.
    Handles common assembler naming conventions.
    """
    hq_ids: Set[str] = set()

    if not hq_dir or not hq_dir.exists():
        log_warning(f"HQ directory not found: {hq_dir}")
        return hq_ids

    for filename in hq_dir.iterdir():
        if filename.name.startswith('.'):
            continue
        file_id = extract_contig_id(filename.stem)
        if file_id:
            hq_ids.add(file_id)

    return hq_ids


# ──────────────────────────────────────────────────────────────────────────────
# Contig splitting
# ──────────────────────────────────────────────────────────────────────────────

def split_contig(
    contig_name: str,
    sequence: str,
    split_positions: List[int],
    min_length: int,
) -> List[Tuple[str, str]]:
    """
    Split a contig at the given positions and return fragments >= min_length.
    Cuts are made at the exact clipping position (no flanking removal).
    """
    fragments = []
    cuts = [0] + split_positions + [len(sequence)]

    for i in range(len(cuts) - 1):
        frag_seq = sequence[cuts[i]:cuts[i + 1]]
        frag_len = len(frag_seq)
        if frag_len >= min_length:
            name = (
                f"{contig_name}_split{i + 1}_len{frag_len}"
                if split_positions else contig_name
            )
            fragments.append((name, frag_seq))

    return fragments


# ──────────────────────────────────────────────────────────────────────────────
# Contig classification
# ──────────────────────────────────────────────────────────────────────────────

def classify_contigs(
    assembly: Path,
    hq_ids: Set[str],
    clipping_positions: Dict[str, List[int]],
    chimera_flagged: Set[str],
    split_hq_circular: bool,
    hq_circular_min_length: int,
) -> Tuple[Dict[str, Set[str]], Dict[str, int]]:
    """
    Classify every contig in the assembly into one of five categories:

      hq_clean          HQ, no clipping, not chimeric  → safe to extract
      hq_chimeric       HQ, no clipping, but chimeric  → must split/keep
      hq_problematic    HQ + clipping zone             → split or remove
      hq_circular_clip  HQ + clipping + circular       → force-split when
                        (only populated when               --split-hq-circular
                        --split-hq-circular is set)
      nonhq_misassembly Non-HQ + clipping              → split or remove
      nonhq_clean       Non-HQ, no clipping            → keep intact

    Also returns a {contig_name: seq_len} dict for downstream use.
    """
    categories: Dict[str, Set[str]] = {
        'hq_clean': set(),
        'hq_chimeric': set(),
        'hq_problematic': set(),
        'hq_circular_clip': set(),
        'nonhq_misassembly': set(),
        'nonhq_clean': set(),
    }
    seq_lengths: Dict[str, int] = {}

    log_info("Classifying contigs…")

    for contig_name, seq in read_fasta_streaming(assembly):
        seq_len = len(seq)
        seq_lengths[contig_name] = seq_len

        contig_id = extract_contig_id(contig_name)
        is_hq = contig_id in hq_ids
        has_clip = contig_name in clipping_positions
        is_chimeric = contig_name in chimera_flagged
        is_circular = (
            split_hq_circular
            and is_hq
            and is_circular_candidate(contig_name, seq_len, hq_circular_min_length)
        )

        if is_hq and has_clip and is_circular:
            categories['hq_circular_clip'].add(contig_name)
        elif is_hq and has_clip:
            categories['hq_problematic'].add(contig_name)
        elif is_hq and is_chimeric:
            categories['hq_chimeric'].add(contig_name)
        elif is_hq:
            categories['hq_clean'].add(contig_name)
        elif has_clip:
            categories['nonhq_misassembly'].add(contig_name)
        else:
            categories['nonhq_clean'].add(contig_name)

    return categories, seq_lengths


# ──────────────────────────────────────────────────────────────────────────────
# Main processing
# ──────────────────────────────────────────────────────────────────────────────

def process_assembly(
    assembly: Path,
    hq_dir: Optional[Path],
    misassembly_dir: Path,
    output: Optional[Path],
    mode: str,
    min_length: int,
    preserve_hq_with_issues: bool,
    dry_run: bool = False,
    # clipping selection
    safe: bool = True,
    min_clip_coverage: int = 0,
    # chimera awareness
    chimera_report: Optional[Path] = None,
    chimera_risk_threshold: str = 'Medium',
    # HQ circular handling
    split_hq_circular: bool = False,
    hq_circular_min_length: int = 200_000,
) -> dict:
    """Main processing function."""

    print_header("MaMISA - Integrated Assembly Filtering")

    # ── 1. HQ IDs ─────────────────────────────────────────────────────────
    print_section("Loading HQ genome IDs")
    hq_ids = get_hq_contig_ids(hq_dir) if hq_dir else set()
    log_info(f"HQ directory: {hq_dir}")
    log_info(f"Found {len(hq_ids):,} HQ genome IDs")

    # ── 2. Clipping positions (with coverage filtering) ────────────────────
    print_section("Loading misassembly data")
    log_info(f"Misassembly directory: {misassembly_dir}")

    if safe:
        log_info("Mode: --safe (using ALL clipping positions reported by anvi'o)")
        effective_min_cov = 0
    else:
        effective_min_cov = min_clip_coverage
        log_info(
            f"Mode: selective (--min-clip-coverage {effective_min_cov}) — "
            f"only positions with ≥{effective_min_cov} clipped reads will be used"
        )

    clipping_positions = load_clipping_positions(misassembly_dir, effective_min_cov)

    if not clipping_positions and not dry_run:
        log_error(f"No clipping data found in {misassembly_dir}")
        sys.exit(1)

    # ── 3. Chimera report ─────────────────────────────────────────────────
    chimera_flagged: Set[str] = set()
    if chimera_report:
        print_section("Loading chimera report")
        chimera_flagged = load_chimera_flagged(chimera_report, chimera_risk_threshold)

    # ── 4. Classify contigs ───────────────────────────────────────────────
    print_section("Classifying contigs")
    categories, seq_lengths = classify_contigs(
        assembly=assembly,
        hq_ids=hq_ids,
        clipping_positions=clipping_positions,
        chimera_flagged=chimera_flagged,
        split_hq_circular=split_hq_circular,
        hq_circular_min_length=hq_circular_min_length,
    )

    print_section("CONTIG CLASSIFICATION")
    print(f"  HQ clean (safe to extract):      {len(categories['hq_clean']):>8,}")
    print(f"  HQ chimeric (no clipping):       {len(categories['hq_chimeric']):>8,}")
    print(f"  HQ with clipping zone:           {len(categories['hq_problematic']):>8,}")
    if split_hq_circular:
        print(f"  HQ circular + clipping (split):  {len(categories['hq_circular_clip']):>8,}")
    print(f"  Non-HQ with clipping zone:       {len(categories['nonhq_misassembly']):>8,}")
    print(f"  Non-HQ clean:                    {len(categories['nonhq_clean']):>8,}")

    # ── 5. Decision logic ─────────────────────────────────────────────────
    print_section("PROCESSING STRATEGY")

    # Always remove clean, non-chimeric HQ genomes (extracted as standalone MAGs)
    to_remove: Set[str] = categories['hq_clean'].copy()

    # Always split non-HQ misassembled contigs (in split mode)
    to_split: Set[str] = categories['nonhq_misassembly'].copy()

    to_keep_intact: Set[str] = categories['nonhq_clean'].copy()

    # HQ circular + clipping → force split (user explicitly requested this)
    if split_hq_circular and categories['hq_circular_clip']:
        to_split.update(categories['hq_circular_clip'])
        log_info(
            f"  --split-hq-circular: {len(categories['hq_circular_clip']):,} "
            f"HQ circular contig(s) with clipping will be split"
        )

    # HQ chimeric (no clipping detected, but reads say chimeric)
    if categories['hq_chimeric']:
        if preserve_hq_with_issues or split_hq_circular:
            to_split.update(categories['hq_chimeric'])
            log_info(
                f"  Chimeric HQ contigs ({len(categories['hq_chimeric']):,}): "
                f"will be split (kept in assembly)"
            )
        else:
            to_remove.update(categories['hq_chimeric'])
            log_warning(
                f"  Chimeric HQ contigs ({len(categories['hq_chimeric']):,}): "
                f"will be REMOVED — use --preserve-hq-with-issues to split them instead"
            )

    # HQ with clipping zone (not circular, or --split-hq-circular not set)
    if categories['hq_problematic']:
        if preserve_hq_with_issues:
            to_split.update(categories['hq_problematic'])
            log_info(
                f"  HQ contigs with clipping ({len(categories['hq_problematic']):,}): "
                f"will be split (--preserve-hq-with-issues)"
            )
        else:
            to_remove.update(categories['hq_problematic'])
            log_warning(
                f"  HQ contigs with clipping ({len(categories['hq_problematic']):,}): "
                f"will be REMOVED — use --preserve-hq-with-issues to split them instead"
            )

    print(f"\n  Actions:")
    print(f"     Remove:       {len(to_remove):>8,}  (HQ clean + any unpreserved HQ problematic)")
    print(f"     Split:        {len(to_split):>8,}  (misassembled / chimeric / HQ circular)")
    print(f"     Keep intact:  {len(to_keep_intact):>8,}  (non-HQ clean)")

    # ── Dry-run early exit ────────────────────────────────────────────────
    if dry_run:
        log_info("\nDRY-RUN MODE: No output file generated")

        examples = list(categories['hq_problematic'])[:5]
        if examples:
            print("\nExamples of HQ contigs with clipping zone:")
            for i, name in enumerate(examples, 1):
                n_pos = len(clipping_positions.get(name, []))
                print(f"  {i}. {name}  ({n_pos} clipping position(s))")
            if len(categories['hq_problematic']) > 5:
                print(f"  … and {len(categories['hq_problematic']) - 5} more")

        if categories['hq_chimeric']:
            print(f"\n  {len(categories['hq_chimeric']):,} HQ contigs flagged as chimeric "
                  f"(no clipping zone, but chimera report)")

        return {
            'hq_clean': len(categories['hq_clean']),
            'hq_chimeric': len(categories['hq_chimeric']),
            'hq_problematic': len(categories['hq_problematic']),
            'hq_circular_clip': len(categories['hq_circular_clip']),
            'nonhq_misassembly': len(categories['nonhq_misassembly']),
            'nonhq_clean': len(categories['nonhq_clean']),
            'to_remove': len(to_remove),
            'to_split': len(to_split),
            'to_keep': len(to_keep_intact),
        }

    # ── 6. Process assembly ───────────────────────────────────────────────
    print_section(f"Processing assembly (mode: {mode.upper()})")

    output_sequences: List[Tuple[str, str]] = []
    stats = {
        'removed': 0, 'split': 0, 'kept_intact': 0,
        'fragments_generated': 0, 'fragments_kept': 0,
        'original_bp': 0, 'final_bp': 0,
    }

    for contig_name, sequence in read_fasta_streaming(assembly):
        contig_len = len(sequence)
        stats['original_bp'] += contig_len

        if contig_name in to_remove:
            stats['removed'] += 1
            continue

        if contig_name in to_split and mode == 'split':
            split_pos = clipping_positions.get(contig_name, [])
            fragments = split_contig(contig_name, sequence, split_pos, min_length)
            stats['split'] += 1
            stats['fragments_generated'] += len(split_pos) + 1
            stats['fragments_kept'] += len(fragments)
            output_sequences.extend(fragments)
            stats['final_bp'] += sum(len(s) for _, s in fragments)

        elif contig_name in to_split and mode == 'remove':
            stats['removed'] += 1

        else:
            if contig_len >= min_length:
                output_sequences.append((contig_name, sequence))
                stats['kept_intact'] += 1
                stats['final_bp'] += contig_len
            else:
                stats['removed'] += 1

    # ── 7. Write output ───────────────────────────────────────────────────
    log_info(f"\nWriting output to: {output}")
    write_fasta(output_sequences, output)

    # ── 8. Final statistics ───────────────────────────────────────────────
    print_section("FINAL STATISTICS")
    print(f"  Contigs removed:         {stats['removed']:>10,}")
    print(f"  Contigs split:           {stats['split']:>10,}")
    print(f"  Contigs kept intact:     {stats['kept_intact']:>10,}")
    print(f"  Total output contigs:    {len(output_sequences):>10,}")
    print(f"\n  Original bases:          {stats['original_bp']:>10,} bp")
    print(f"  Final bases:             {stats['final_bp']:>10,} bp")
    pct = 100 * stats['final_bp'] / stats['original_bp'] if stats['original_bp'] else 0
    print(f"  Retained:                {pct:>9.2f}%")

    if chimera_flagged:
        n_removed = len(categories['hq_chimeric'] - to_split)
        n_split = len(categories['hq_chimeric'] & to_split)
        if n_split:
            print(f"\n  ⚠  {n_split:,} chimeric HQ contig(s) were split (kept in assembly)")
        if n_removed:
            print(f"  ⚠  {n_removed:,} chimeric HQ contig(s) were removed")

    return stats


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def register_parser(subparsers):
    parser = subparsers.add_parser(
        'filter-misassemblies',
        help='Filter assembly based on misassembly detection',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Clipping position selection
---------------------------
  --safe (default)
      Use all positions reported by anvi'o (conservative).

  --min-clip-coverage N
      Only split at positions with ≥N clipped reads (selective).
      Higher N = less aggressive. Mutually exclusive with --safe.

Examples
--------
  # Conservative: split at all anvi'o-reported positions
  mamisa filter-misassemblies \\
    --assembly assembly.fa --misassemblies misasm_dir/ \\
    --output clean.fa --safe

  # Selective: only split where ≥50 reads are clipped
  mamisa filter-misassemblies \\
    --assembly assembly.fa --misassemblies misasm_dir/ \\
    --output clean.fa --min-clip-coverage 50

  # Full pipeline: chimera-aware + force-split HQ circular genomes
  mamisa filter-misassemblies \\
    --assembly assembly.fa --misassemblies misasm_dir/ \\
    --hq-genomes hq_dir/ --output clean.fa \\
    --chimera-report chimera_results/chimera_read_report.tsv \\
    --split-hq-circular --preserve-hq-with-issues --safe
        """
    )

    parser.add_argument('-a', '--assembly', type=Path, required=True,
                        help='Input assembly FASTA file')
    parser.add_argument('-m', '--misassemblies', type=Path, required=True,
                        help='Directory with *-clipping.txt files from anvi\'o')
    parser.add_argument('-g', '--hq-genomes', type=Path,
                        help='Directory with HQ genome files (from process-large-contigs)')
    parser.add_argument('-o', '--output', type=Path,
                        help='Output filtered assembly FASTA')
    parser.add_argument('-l', '--min-length', type=int, default=2500,
                        help='Minimum contig/fragment length to keep (default: 2500 bp)')
    parser.add_argument('--mode', choices=['remove', 'split'], default='split',
                        help='How to handle misassembled non-HQ contigs (default: split)')
    parser.add_argument('--preserve-hq-with-issues', action='store_true',
                        help='Split HQ contigs with clipping zones instead of removing them')

    # ── Clipping position selection ────────────────────────────────────────
    clip_group = parser.add_mutually_exclusive_group()
    clip_group.add_argument('--safe', action='store_true', default=True,
                            help='Use ALL clipping positions reported by anvi\'o (default)')
    clip_group.add_argument('--min-clip-coverage', type=int, default=0,
                            metavar='N',
                            help='Only split at positions where ≥N reads are clipped; '
                                 'overrides --safe')

    # ── Chimera awareness ─────────────────────────────────────────────────
    parser.add_argument('--chimera-report', type=Path,
                        help='TSV from check-chimeras or check-read-chimeras; '
                             'chimeric HQ contigs are treated as problematic')
    parser.add_argument('--chimera-risk-threshold',
                        choices=['Low', 'Medium', 'High'], default='Medium',
                        help='Minimum chimera risk level to act on (default: Medium)')

    # ── HQ circular handling ───────────────────────────────────────────────
    parser.add_argument('--split-hq-circular', action='store_true',
                        help='Force-split HQ circular contigs that have a clipping zone '
                             'instead of removing them')
    parser.add_argument('--hq-circular-min-length', type=int, default=200_000,
                        metavar='BP',
                        help='Length threshold for circular candidate detection '
                             '(default: 200000 bp). Also detects circular=true/Y in headers.')

    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would happen without writing output')
    parser.add_argument('--stats', type=Path,
                        help='Save processing statistics to a TSV file')

    parser.set_defaults(func=run)
    return parser


def run(args):
    """Execute the filter-misassemblies command."""

    if not args.dry_run and not args.output:
        raise ValueError("--output is required unless using --dry-run")

    validate_file_exists(args.assembly, "Assembly file")
    validate_dir_exists(args.misassemblies, "Misassemblies directory")

    if args.chimera_report and not args.chimera_report.exists():
        log_error(f"Chimera report not found: {args.chimera_report}")
        sys.exit(1)

    # --min-clip-coverage overrides --safe
    use_safe = args.safe and args.min_clip_coverage == 0

    stats = process_assembly(
        assembly=args.assembly,
        hq_dir=args.hq_genomes,
        misassembly_dir=args.misassemblies,
        output=args.output,
        mode=args.mode,
        min_length=args.min_length,
        preserve_hq_with_issues=args.preserve_hq_with_issues,
        dry_run=args.dry_run,
        safe=use_safe,
        min_clip_coverage=args.min_clip_coverage,
        chimera_report=args.chimera_report,
        chimera_risk_threshold=args.chimera_risk_threshold,
        split_hq_circular=args.split_hq_circular,
        hq_circular_min_length=args.hq_circular_min_length,
    )

    if args.stats and not args.dry_run:
        with open(args.stats, 'w') as f:
            f.write("metric\tvalue\n")
            for key, value in stats.items():
                f.write(f"{key}\t{value}\n")
        log_info(f"\nStatistics saved to: {args.stats}")

    log_info("\n✓ Done!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    register_parser(parser.add_subparsers(dest='command'))
    run(parser.parse_args())
