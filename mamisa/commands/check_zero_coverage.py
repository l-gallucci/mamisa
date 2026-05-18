#!/usr/bin/env python3
"""
MaMISA - check-zero-coverage command

Identify assembly regions that have no read coverage and validate whether
those regions are actually supported by source reads via BLAST and (optionally)
k-mer composition.

Pipeline
--------
  1. Find zero-coverage intervals
       --zero-cov-dir   use anvi'o *-zero_cov.txt files (fast, already computed)
       --bam            compute directly from BAM with samtools depth -a
  2. Extract sequences for intervals >= --min-length
  3. Build a nucleotide BLAST database from the long reads (--reads)
  4. blastn -dust no (following the anvi'o long-read benchmarking workflow)
       Regions with zero hits = assembler-invented sequence
  5. (optional) Meryl 21-mer support fraction (--meryl)
       Fraction of region k-mers present in source reads.

Output TSV columns
------------------
  contig  start  end  length  n_blast_hits  kmer_support_pct  verdict
  verdict: supported | partial | unsupported

Reference: https://merenlab.org/data/benchmarking-long-read-assemblers/
"""

import csv
import sys
import argparse
from pathlib import Path
from typing import Optional

from ..utils.zero_coverage import (
    check_blast, check_meryl,
    parse_zero_cov_files, compute_zero_cov_from_bam,
    extract_region_sequences, _write_region_fasta,
    prepare_reads_fasta, make_blast_db,
    blast_regions_vs_reads, run_meryl_kmer_support,
    assign_verdict,
)
from ..utils.validation import validate_file_exists
from ..utils.logging import log_info, log_warning, log_error, print_header, print_section


# ──────────────────────────────────────────────────────────────────────────────
# Core function
# ──────────────────────────────────────────────────────────────────────────────

def run_check_zero_coverage(
    assembly: Path,
    reads: Path,
    output: Path,
    bam: Optional[Path],
    zero_cov_dir: Optional[Path],
    min_length: int,
    min_mapq: int,
    evalue: float,
    threads: int,
    max_target_seqs: int,
    use_meryl: bool,
    meryl_k: int,
    keep_db: bool,
    dry_run: bool,
) -> None:

    print_header("MaMISA - Check Zero-Coverage Regions")

    # ── 1. Dependency checks ──────────────────────────────────────────────
    if not check_blast():
        sys.exit(1)
    if use_meryl and not check_meryl():
        log_warning("Disabling Meryl k-mer check (meryl not found)")
        use_meryl = False

    # ── 2. Find zero-coverage regions ────────────────────────────────────
    print_section("Finding zero-coverage regions")

    if zero_cov_dir:
        log_info(f"Source: anvi'o *-zero_cov.txt files in {zero_cov_dir}")
        regions = parse_zero_cov_files(zero_cov_dir)
    else:
        log_info(f"Source: samtools depth -a on {bam.name}")
        regions = compute_zero_cov_from_bam(bam, min_mapq=min_mapq,
                                            min_length=min_length)

    # Apply length filter
    regions = [(c, s, e) for c, s, e in regions if e - s >= min_length]
    log_info(f"Regions ≥{min_length} bp: {len(regions):,}")

    if not regions:
        log_warning("No zero-coverage regions found — nothing to validate")
        _write_empty(output)
        return

    # ── 3. Extract sequences ──────────────────────────────────────────────
    print_section("Extracting region sequences from assembly")
    region_seqs = extract_region_sequences(assembly, regions, min_length)

    if not region_seqs:
        log_error("No sequences could be extracted — check assembly FASTA")
        sys.exit(1)

    if dry_run:
        log_info(f"\nDRY-RUN: would BLAST {len(region_seqs):,} regions against reads")
        log_info(f"  Reads file: {reads}")
        log_info(f"  Output:     {output}")
        _print_region_summary(regions)
        return

    # ── 4. Prepare work directory ─────────────────────────────────────────
    work_dir = output.parent / f'{output.stem}_blast_tmp'
    work_dir.mkdir(parents=True, exist_ok=True)

    query_fa = work_dir / 'query_regions.fasta'
    _write_region_fasta(region_seqs, query_fa)

    # ── 5. Build BLAST database from reads ────────────────────────────────
    print_section("Building BLAST database from reads")
    reads_fa = prepare_reads_fasta(reads, work_dir)
    db_path  = work_dir / 'reads_blast_db'

    if not make_blast_db(reads_fa, db_path):
        sys.exit(1)

    # ── 6. BLAST query regions against reads ──────────────────────────────
    print_section("BLASTing zero-coverage regions against reads")
    log_info(f"  evalue={evalue}, max_target_seqs={max_target_seqs}, "
             f"threads={threads}, dust=no")
    hits = blast_regions_vs_reads(
        query_fa, db_path,
        threads=threads,
        evalue=evalue,
        max_target_seqs=max_target_seqs,
    )

    # ── 7. Meryl k-mer support (optional) ─────────────────────────────────
    kmer_support = {}
    if use_meryl:
        print_section("Meryl k-mer support validation")
        meryl_dir = work_dir / 'meryl'
        kmer_support = run_meryl_kmer_support(
            region_seqs, reads_fa,
            k=meryl_k, threads=threads, work_dir=meryl_dir,
        )

    # ── 8. Write output ───────────────────────────────────────────────────
    print_section("Writing report")
    output.parent.mkdir(parents=True, exist_ok=True)

    verdicts = {'supported': 0, 'partial': 0, 'unsupported': 0}

    with open(output, 'w', newline='') as f:
        fieldnames = ['contig', 'start', 'end', 'length',
                      'n_blast_hits', 'kmer_support_pct', 'verdict']
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        for rid, (contig, start, end, _) in sorted(region_seqs.items()):
            length       = end - start
            n_hits       = hits.get(rid, 0)
            kmer_frac    = kmer_support.get(rid, float('nan'))
            import math
            kmer_pct_str = (f"{100 * kmer_frac:.1f}"
                            if not math.isnan(kmer_frac) else 'N/A')
            verdict = assign_verdict(n_hits, kmer_frac, use_meryl)
            verdicts[verdict] = verdicts.get(verdict, 0) + 1

            writer.writerow({
                'contig':           contig,
                'start':            start,
                'end':              end,
                'length':           length,
                'n_blast_hits':     n_hits,
                'kmer_support_pct': kmer_pct_str,
                'verdict':          verdict,
            })

    # ── 9. Clean up ───────────────────────────────────────────────────────
    if not keep_db:
        import shutil
        shutil.rmtree(work_dir, ignore_errors=True)
    else:
        log_info(f"BLAST database kept at: {work_dir}")

    # ── 10. Summary ───────────────────────────────────────────────────────
    print_section("SUMMARY")
    total = sum(verdicts.values())
    print(f"  {'unsupported':<14}  {verdicts.get('unsupported', 0):>6,}  "
          f"no read backing — likely assembler artefact")
    print(f"  {'partial':<14}  {verdicts.get('partial', 0):>6,}  "
          f"weak read support — inspect manually")
    print(f"  {'supported':<14}  {verdicts.get('supported', 0):>6,}  "
          f"read-supported zero-cov (coverage gap, mapping artefact)")
    print(f"\n  Total regions: {total:,}")
    print(f"  Report:        {output}")

    n_bad = verdicts.get('unsupported', 0)
    if n_bad > 0:
        log_warning(
            f"{n_bad:,} unsupported regions detected — these sequences have "
            f"no read backing and may represent assembler hallucinations. "
            f"Consider masking or removing them before downstream analysis."
        )


def _write_empty(output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, 'w') as f:
        f.write("contig\tstart\tend\tlength\tn_blast_hits\tkmer_support_pct\tverdict\n")


def _print_region_summary(regions) -> None:
    lengths = [e - s for _, s, e in regions]
    if lengths:
        print(f"\n  Regions: {len(lengths):,}")
        print(f"  Min length:  {min(lengths):,} bp")
        print(f"  Max length:  {max(lengths):,} bp")
        print(f"  Mean length: {sum(lengths) // len(lengths):,} bp")
        print(f"  Total bases: {sum(lengths):,} bp")


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def register_parser(subparsers):
    parser = subparsers.add_parser(
        'check-zero-coverage',
        help='Identify assembler-invented sequence via BLAST and k-mer validation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
What this detects
-----------------
Assembly regions with zero read coverage can be genuine (e.g. highly variable
or very low-abundance loci) or artefacts generated by the assembler with no
source-read support. BLAST against the original reads distinguishes them:

  unsupported   No read in the dataset covers this sequence.
                Likely an assembler hallucination (chimeric join, repeat
                resolution error, or graph traversal artefact).
  partial       Some hits found but weak. Warrants manual inspection.
  supported     Read-covered region that the mapper failed to align.
                Usually a mapping artefact; sequence itself is real.

Zero-coverage source (choose one)
----------------------------------
  --zero-cov-dir   Use *-zero_cov.txt files from anvi-script-find-misassemblies.
                   Faster — the intervals are already computed.
  --bam            Compute directly from the BAM with samtools depth -a.
                   Use when anvi'o output is not available.

Examples
--------
  # Using anvi'o output
  mamisa check-zero-coverage \\
    --assembly assembly.fa \\
    --reads long_reads.fq.gz \\
    --zero-cov-dir misassemblies/ \\
    --output zero_cov_report.tsv

  # Using BAM directly, with Meryl k-mer cross-check
  mamisa check-zero-coverage \\
    --assembly assembly.fa \\
    --reads long_reads.fa \\
    --bam mapping.bam \\
    --meryl \\
    --output zero_cov_report.tsv
        """
    )

    parser.add_argument('-a', '--assembly', type=Path, required=True,
                        help='Assembly FASTA file')
    parser.add_argument('--reads', type=Path, required=True,
                        help='Long-read source file (FASTA or FASTQ, optionally gzip)')
    parser.add_argument('-o', '--output', type=Path, required=True,
                        help='Output TSV report')

    # Zero-cov source (mutually exclusive)
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument('--zero-cov-dir', type=Path, metavar='DIR',
                     help='Directory with *-zero_cov.txt files from anvi\'o')
    src.add_argument('--bam', type=Path,
                     help='Sorted, indexed BAM file (computes zero-cov on the fly)')

    parser.add_argument('--min-length', type=int, default=500, metavar='BP',
                        help='Minimum zero-coverage region length to analyse (default: 500)')
    parser.add_argument('--min-mapq', type=int, default=20,
                        help='Min mapping quality for depth computation (--bam only; default: 20)')
    parser.add_argument('--evalue', type=float, default=1e-5,
                        help='BLAST e-value cutoff (default: 1e-5)')
    parser.add_argument('--max-target-seqs', type=int, default=1,
                        help='BLAST max target seqs per query (default: 1)')
    parser.add_argument('--threads', type=int, default=4,
                        help='Threads for BLAST and Meryl (default: 4)')

    parser.add_argument('--meryl', action='store_true',
                        help='Enable Meryl k-mer support validation (requires meryl in PATH)')
    parser.add_argument('--meryl-k', type=int, default=21,
                        help='k-mer size for Meryl (default: 21)')
    parser.add_argument('--keep-db', action='store_true',
                        help='Keep BLAST database files after run')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be done without running BLAST')

    parser.set_defaults(func=run)
    return parser


def run(args):
    validate_file_exists(args.assembly, "Assembly")
    validate_file_exists(args.reads, "Reads file")

    if args.bam:
        validate_file_exists(args.bam, "BAM file")
    if args.zero_cov_dir and not args.zero_cov_dir.exists():
        log_error(f"Zero-coverage directory not found: {args.zero_cov_dir}")
        sys.exit(1)

    run_check_zero_coverage(
        assembly      = args.assembly,
        reads         = args.reads,
        output        = args.output,
        bam           = args.bam,
        zero_cov_dir  = args.zero_cov_dir,
        min_length    = args.min_length,
        min_mapq      = args.min_mapq,
        evalue        = args.evalue,
        threads       = args.threads,
        max_target_seqs = args.max_target_seqs,
        use_meryl     = args.meryl,
        meryl_k       = args.meryl_k,
        keep_db       = args.keep_db,
        dry_run       = args.dry_run,
    )
