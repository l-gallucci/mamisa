#!/usr/bin/env python3
"""
MaMISA - process-large-contigs command
Extract and QC large contigs, then intelligently filter based on quality and misassemblies
"""

import re
import sys
import argparse
import subprocess
import shutil
from pathlib import Path
from typing import Set, Dict, List, Tuple

from ..utils.fasta import read_fasta_streaming, write_fasta
from ..utils.misassembly import get_contigs_with_misassemblies
from ..utils.checkm2 import parse_quality_report
from ..utils.validation import validate_file_exists, validate_dir_exists, check_dependencies
from ..utils.logging import log_info, log_error, log_warning, print_header, print_section


def safe_filename(name: str) -> str:
    """
    Sanitize a contig name for use as a filename.
    Replaces any character that is not alphanumeric, dot, hyphen, or underscore.
    """
    return re.sub(r'[^\w.\-]', '_', name)


def extract_large_contigs(assembly: Path, max_length: int, output_dir: Path) -> Tuple[int, int]:
    """
    Separate contigs into large (> max_length) and regular (<= max_length).

    Returns:
        (n_large, n_regular)
    """
    large_contigs = []
    regular_contigs = []

    for name, seq in read_fasta_streaming(assembly):
        if len(seq) > max_length:
            large_contigs.append((name, seq))
        else:
            regular_contigs.append((name, seq))

    large_file = output_dir / "large_contigs.fa"
    regular_file = output_dir / "assembly_regular.fa"

    if large_contigs:
        write_fasta(large_contigs, large_file)
        log_info(f"Extracted {len(large_contigs):,} large contigs to: {large_file}")
    else:
        log_warning("No large contigs found")

    write_fasta(regular_contigs, regular_file)
    log_info(f"Regular contigs: {len(regular_contigs):,} → {regular_file}")

    return len(large_contigs), len(regular_contigs)


def split_large_contigs_individually(large_contigs_file: Path, output_dir: Path) -> int:
    """
    Write each large contig to its own FASTA file for CheckM2 input.

    Returns:
        Number of files created
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    count = 0

    for name, seq in read_fasta_streaming(large_contigs_file):
        out_file = output_dir / f"{safe_filename(name)}.fa"
        write_fasta([(name, seq)], out_file)
        count += 1

    log_info(f"Split into {count:,} individual files")
    return count


def run_checkm2(genome_dir: Path, output_dir: Path, threads: int = 1,
                checkm2_env: str = None) -> bool:
    """
    Run CheckM2 on a directory of genome files.
    Streams CheckM2 stdout/stderr live so the user can monitor progress.

    Returns:
        True if CheckM2 completed successfully
    """
    deps = check_dependencies(['checkm2'])

    if deps['checkm2'] is None:
        log_error("CheckM2 not found in PATH")
        if checkm2_env:
            log_error(f"Please activate the conda environment: conda activate {checkm2_env}")
        return False

    log_info(f"Found CheckM2: {deps['checkm2']}")

    cmd = [
        'checkm2', 'predict',
        '--threads', str(threads),
        '--input', str(genome_dir),
        '--output-directory', str(output_dir),
        '-x', 'fa',
        '--force',
    ]

    log_info(f"Running: {' '.join(cmd)}")

    try:
        # Do NOT use capture_output=True — let stdout/stderr stream to the terminal
        # so the user can see CheckM2 progress in real time.
        subprocess.run(cmd, check=True)
        log_info("CheckM2 completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        log_error(f"CheckM2 failed with exit code {e.returncode}")
        return False


def make_filtering_decisions(quality_data: Dict[str, Dict],
                             misasm_contigs: Set[str],
                             min_completeness: float,
                             max_contamination: float) -> Dict[str, List[str]]:
    """
    Decide what to do with each large contig based on quality and misassembly status.

    Categories:
      - extract_hq:       HQ and clean → extract as genome
      - keep_hq_misasm:   HQ but misassembled → keep in assembly for splitting
      - keep_low_quality: Below HQ threshold → keep in assembly
    """
    decisions: Dict[str, List[str]] = {
        'extract_hq': [],
        'keep_hq_misasm': [],
        'keep_low_quality': [],
    }

    print_section("LARGE CONTIGS ANALYSIS")

    for contig_name, qc_data in quality_data.items():
        comp = qc_data['completeness']
        cont = qc_data['contamination']
        passes_qc = comp >= min_completeness and cont <= max_contamination
        has_misasm = contig_name in misasm_contigs

        # Use log_info so output goes to stderr (not stdout)
        log_info(f"\n{contig_name}:")
        log_info(f"  Completeness:  {comp:>6.1f}%")
        log_info(f"  Contamination: {cont:>6.1f}%")
        log_info(f"  Misassemblies: {'YES' if has_misasm else 'NO'}")

        if passes_qc and not has_misasm:
            decisions['extract_hq'].append(contig_name)
            log_info("  → EXTRACT (HQ clean genome)")
        elif passes_qc and has_misasm:
            decisions['keep_hq_misasm'].append(contig_name)
            log_info("  → KEEP (HQ with misassemblies - will be split)")
        else:
            decisions['keep_low_quality'].append(contig_name)
            log_info("  → KEEP (low quality)")

    print_section("DECISION SUMMARY")
    log_info(f"  HQ to extract:               {len(decisions['extract_hq']):>6,}")
    log_info(f"  HQ with misasm (keep):       {len(decisions['keep_hq_misasm']):>6,}")
    log_info(f"  Low quality (keep):          {len(decisions['keep_low_quality']):>6,}")

    return decisions


def save_decisions(decisions: Dict[str, List[str]], quality_data: Dict[str, Dict],
                   misasm_contigs: Set[str], output_file: Path):
    """Save per-contig filtering decisions to a TSV file."""

    with open(output_file, 'w') as f:
        f.write("contig\tdecision\tcompleteness\tcontamination\thas_misassembly\n")

        for contig in decisions['extract_hq']:
            qc = quality_data[contig]
            f.write(f"{contig}\tEXTRACT_HQ\t{qc['completeness']:.2f}"
                    f"\t{qc['contamination']:.2f}\tFalse\n")

        for contig in decisions['keep_hq_misasm']:
            qc = quality_data[contig]
            f.write(f"{contig}\tKEEP_HQ_MISASM\t{qc['completeness']:.2f}"
                    f"\t{qc['contamination']:.2f}\tTrue\n")

        for contig in decisions['keep_low_quality']:
            qc = quality_data[contig]
            has_misasm = contig in misasm_contigs
            f.write(f"{contig}\tKEEP_LOW_QUAL\t{qc['completeness']:.2f}"
                    f"\t{qc['contamination']:.2f}\t{has_misasm}\n")

    log_info(f"Decisions saved to: {output_file}")


def extract_hq_genomes(large_contigs_individual_dir: Path, hq_list: List[str],
                       output_dir: Path) -> int:
    """Copy HQ genome files to a dedicated output directory."""
    output_dir.mkdir(parents=True, exist_ok=True)
    extracted = 0

    for contig_name in hq_list:
        source = large_contigs_individual_dir / f"{safe_filename(contig_name)}.fa"
        dest = output_dir / source.name

        if source.exists():
            shutil.copy2(source, dest)
            extracted += 1
        else:
            log_warning(f"Source file not found for {contig_name}: {source}")

    log_info(f"Extracted {extracted:,} HQ genomes to: {output_dir}")
    return extracted


def create_updated_assembly(regular_assembly: Path, large_contigs_individual_dir: Path,
                            contigs_to_keep: List[str], output_file: Path):
    """Create updated assembly: all regular contigs + selected large contigs."""
    sequences = list(read_fasta_streaming(regular_assembly))
    log_info(f"Added {len(sequences):,} regular contigs")

    added = 0
    for contig_name in contigs_to_keep:
        contig_file = large_contigs_individual_dir / f"{safe_filename(contig_name)}.fa"
        if contig_file.exists():
            for name, seq in read_fasta_streaming(contig_file):
                sequences.append((name, seq))
                added += 1
        else:
            log_warning(f"Contig file not found: {contig_file}")

    log_info(f"Added {added:,} large contigs to keep")
    write_fasta(sequences, output_file)
    log_info(f"Updated assembly: {len(sequences):,} contigs → {output_file}")


def register_parser(subparsers):
    """Register this command's parser"""
    parser = subparsers.add_parser(
        'process-large-contigs',
        help='Extract, QC, and filter large contigs intelligently',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Complete workflow
  mamisa process-large-contigs \\
    --assembly assembly.fa \\
    --misassemblies misasm_dir/ \\
    --output-dir filtered/ \\
    --max-length 300000 \\
    --min-completeness 50 \\
    --max-contamination 10 \\
    --threads 40

  # Skip CheckM2 if already run
  mamisa process-large-contigs \\
    --assembly assembly.fa \\
    --misassemblies misasm_dir/ \\
    --output-dir filtered/ \\
    --skip-checkm2 \\
    --checkm2-results filtered/03_checkm2/
        """
    )

    parser.add_argument('-a', '--assembly', type=Path, required=True,
                        help='Input assembly FASTA file')
    parser.add_argument('-m', '--misassemblies', type=Path, required=True,
                        help="Directory with anvi'o *-clipping.txt files")
    parser.add_argument('-o', '--output-dir', type=Path, required=True,
                        help='Output directory')

    parser.add_argument('--max-length', type=int, default=300000,
                        help='Contigs longer than this are considered large (default: 300000 bp)')
    parser.add_argument('--min-completeness', type=float, default=50.0,
                        help='Minimum completeness to classify a contig as HQ (default: 50%%)')
    parser.add_argument('--max-contamination', type=float, default=10.0,
                        help='Maximum contamination to classify a contig as HQ (default: 10%%)')

    parser.add_argument('--skip-checkm2', action='store_true',
                        help='Skip CheckM2 step (requires --checkm2-results)')
    parser.add_argument('--checkm2-results', type=Path,
                        help='Existing CheckM2 output directory (used with --skip-checkm2)')
    parser.add_argument('--threads', type=int, default=1,
                        help='Number of threads for CheckM2 (default: 1)')
    parser.add_argument('--checkm2-env', default='checkm2',
                        help='CheckM2 conda environment name (default: checkm2)')

    parser.set_defaults(func=run)
    return parser


def run(args):
    """Execute the process-large-contigs command"""

    validate_file_exists(args.assembly, "Assembly file")
    validate_dir_exists(args.misassemblies, "Misassemblies directory")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print_header("MaMISA - Process Large Contigs")

    # Step 1: Extract large contigs
    print_section("STEP 1: Extracting Large Contigs")
    extract_dir = args.output_dir / "01_extracted"
    extract_dir.mkdir(exist_ok=True)

    n_large, n_regular = extract_large_contigs(args.assembly, args.max_length, extract_dir)

    if n_large == 0:
        log_warning("No large contigs found — nothing to process")
        log_info("Assembly can be used as-is for binning")
        return

    # Step 2: Split large contigs into individual files
    print_section("STEP 2: Splitting Large Contigs into Individual Files")
    individual_dir = args.output_dir / "02_individual"
    large_contigs_file = extract_dir / "large_contigs.fa"
    split_large_contigs_individually(large_contigs_file, individual_dir)

    # Step 3: Run CheckM2
    print_section("STEP 3: Quality Assessment (CheckM2)")
    checkm2_dir = args.output_dir / "03_checkm2"

    if args.skip_checkm2:
        if not args.checkm2_results:
            log_error("--skip-checkm2 requires --checkm2-results")
            sys.exit(1)
        checkm2_dir = args.checkm2_results
        log_info(f"Using existing CheckM2 results: {checkm2_dir}")
    else:
        log_info(f"Running CheckM2 with {args.threads} thread(s)...")
        if not run_checkm2(individual_dir, checkm2_dir, args.threads, args.checkm2_env):
            log_error("CheckM2 failed — aborting")
            sys.exit(1)

    # Step 4: Parse CheckM2 results
    print_section("STEP 4: Parsing CheckM2 Results")
    quality_file = checkm2_dir / "quality_report.tsv"

    if not quality_file.exists():
        log_error(f"CheckM2 quality report not found: {quality_file}")
        sys.exit(1)

    raw_records = parse_quality_report(quality_file)
    if not raw_records:
        log_error("No quality data found in CheckM2 report")
        sys.exit(1)

    quality_data = {rec['name']: rec for rec in raw_records}
    log_info(f"Parsed quality data for {len(quality_data):,} contigs")

    # Step 5: Parse misassembly data
    print_section("STEP 5: Loading Misassembly Data")
    misasm_contigs = get_contigs_with_misassemblies(args.misassemblies)

    # Step 6: Make filtering decisions
    decisions = make_filtering_decisions(
        quality_data, misasm_contigs,
        args.min_completeness, args.max_contamination,
    )

    # Step 7: Save decisions
    print_section("STEP 7: Saving Results")
    decisions_file = args.output_dir / "filtering_decisions.tsv"
    save_decisions(decisions, quality_data, misasm_contigs, decisions_file)

    # Step 8: Extract HQ genomes
    hq_dir = args.output_dir / "HQ_extracted"
    n_extracted = extract_hq_genomes(individual_dir, decisions['extract_hq'], hq_dir)

    # Step 9: Create updated assembly
    contigs_to_keep = decisions['keep_hq_misasm'] + decisions['keep_low_quality']
    updated_assembly = args.output_dir / "assembly_for_filtering.fa"
    create_updated_assembly(
        extract_dir / "assembly_regular.fa",
        individual_dir,
        contigs_to_keep,
        updated_assembly,
    )

    # Final summary
    print_section("RESULTS")
    log_info(f"  Large contigs processed:     {n_large:>10,}")
    log_info(f"  HQ genomes extracted:        {n_extracted:>10,}")
    log_info(f"  Contigs in updated assembly: {n_regular + len(contigs_to_keep):>10,}")
    log_info(f"\n  Decisions:         {decisions_file}")
    log_info(f"  HQ genomes:        {hq_dir}/")
    log_info(f"  Updated assembly:  {updated_assembly}")

    log_info("\nNext step — run filter-misassemblies on the updated assembly:")
    log_info(f"  mamisa filter-misassemblies \\")
    log_info(f"    --assembly {updated_assembly} \\")
    log_info(f"    --misassemblies {args.misassemblies} \\")
    log_info(f"    --hq-genomes {hq_dir} \\")
    log_info(f"    --output final_clean_assembly.fa \\")
    log_info(f"    --mode split \\")
    log_info(f"    --preserve-hq-with-issues")

    log_info("\n✓ Done!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    register_parser(parser.add_subparsers(dest='command'))
    args = parser.parse_args()
    run(args)
