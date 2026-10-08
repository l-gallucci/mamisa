#!/usr/bin/env python3
"""
MaMISA - run-checkm2 command
Wrapper for running CheckM2 completeness/contamination prediction
"""

import sys
import shlex
import subprocess
import argparse
from pathlib import Path
from typing import Optional

from ..utils.validation import validate_dir_exists, check_dependencies, check_db_path
from ..utils.logging import log_info, log_error, log_warning, print_header, print_section


def check_checkm2_available() -> bool:
    """Check if checkm2 is available in PATH"""
    deps = check_dependencies(['checkm2'])

    if deps['checkm2'] is None:
        log_error("checkm2 not found in PATH")
        log_error("Please install CheckM2: https://github.com/chklovski/CheckM2")
        return False

    log_info(f"Found checkm2: {deps['checkm2']}")
    return True


def get_checkm2_version() -> Optional[str]:
    """Get CheckM2 version"""
    try:
        result = subprocess.run(['checkm2', '--version'],
                                capture_output=True, text=True, check=True)
        return result.stdout.strip()
    except Exception:
        return None


def run_checkm2_predict(genome_dir: Path, output_dir: Path,
                        extension: str, threads: int,
                        database: Optional[Path] = None,
                        force: bool = False,
                        extra_args: list = None) -> int:
    """
    Run CheckM2 predict.

    Returns:
        Exit code from checkm2
    """

    cmd = [
        'checkm2', 'predict',
        '--input', str(genome_dir),
        '--output-directory', str(output_dir),
        '--threads', str(threads),
        '-x', extension,
    ]

    if database:
        cmd.extend(['--database_path', str(database)])
    if force:
        cmd.append('--force')
    if extra_args:
        cmd.extend(extra_args)

    log_info(f"\nRunning command:")
    log_info(f"  {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, check=True)
        return result.returncode
    except subprocess.CalledProcessError as e:
        log_error(f"CheckM2 failed with exit code {e.returncode}")
        return e.returncode
    except KeyboardInterrupt:
        log_warning("\nInterrupted by user")
        return 130


def process_tier_directory(base_dir: Path, tier: str, output_base: Path,
                           extension: str, threads: int,
                           database: Optional[Path] = None,
                           force: bool = False,
                           extra_args: list = None) -> dict:
    """Process a single tier subdirectory (HQ/MQ/LQ)."""

    tier_genome_dir = base_dir / tier
    tier_output_dir = output_base / tier

    if not tier_genome_dir.exists():
        log_warning(f"Tier directory not found: {tier_genome_dir}")
        return {'status': 'skipped', 'reason': 'directory_not_found'}

    genome_files = list(tier_genome_dir.glob(f"*.{extension}"))
    n_genomes = len(genome_files)

    if n_genomes == 0:
        log_warning(f"No genomes found in {tier_genome_dir} with extension .{extension}")
        return {'status': 'skipped', 'reason': 'no_genomes', 'n_genomes': 0}

    log_info(f"\nProcessing tier: {tier}")
    log_info(f"  Genomes: {n_genomes:,}")
    log_info(f"  Input: {tier_genome_dir}")
    log_info(f"  Output: {tier_output_dir}")

    tier_output_dir.mkdir(parents=True, exist_ok=True)

    exit_code = run_checkm2_predict(
        genome_dir=tier_genome_dir,
        output_dir=tier_output_dir,
        extension=extension,
        threads=threads,
        database=database,
        force=force,
        extra_args=extra_args,
    )

    return {
        'status': 'completed' if exit_code == 0 else 'failed',
        'exit_code': exit_code,
        'n_genomes': n_genomes,
    }


def register_parser(subparsers):
    """Register this command's parser"""
    parser = subparsers.add_parser(
        'run-checkm2',
        help='Run CheckM2 completeness/contamination prediction',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run on a single genome directory
  mamisa run-checkm2 \\
    --genome-dir genomes/ \\
    --output checkm2_results/ \\
    --extension fa \\
    --threads 40

  # Run on all tiers in a Selected directory
  mamisa run-checkm2 \\
    --selected-dir filtered/Selected/ \\
    --output checkm2_results/ \\
    --extension fa \\
    --threads 40

  # With an explicit database path
  mamisa run-checkm2 \\
    --genome-dir genomes/ \\
    --output checkm2_results/ \\
    --database /path/to/CheckM2_database/uniref100.KO.1.dmnd
        """
    )

    # Input options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--selected-dir', type=Path,
                             help='Directory with tier subdirectories (HQ/MQ/LQ)')
    input_group.add_argument('--genome-dir', type=Path,
                             help='Single directory with genome files')

    parser.add_argument('-o', '--output', type=Path, required=True,
                        help='Output directory for CheckM2 results')

    # CheckM2 options
    parser.add_argument('--extension', default='fa',
                        help='Genome file extension (default: fa)')
    parser.add_argument('--threads', type=int, default=1,
                        help='Number of threads to use (default: 1)')
    parser.add_argument('--database', type=Path,
                        help='Path to CheckM2 diamond database (optional)')
    parser.add_argument('--force', action='store_true',
                        help='Overwrite existing CheckM2 output')

    # Tier selection (for --selected-dir mode)
    parser.add_argument('--tiers', default='HQ,MQ,LQ',
                        help='Comma-separated tiers to process (default: HQ,MQ,LQ)')

    # Extra arguments to pass to checkm2
    parser.add_argument('--checkm2-args', default='',
                        help='Additional arguments to pass to checkm2 (quoted string)')

    parser.set_defaults(func=run)
    return parser


def run(args):
    """Execute the run-checkm2 command"""

    print_header("MaMISA - Run CheckM2")

    print_section("Checking Dependencies")
    if not check_checkm2_available():
        sys.exit(1)

    version = get_checkm2_version()
    if version:
        log_info(f"CheckM2 version: {version}")

    extra_args = shlex.split(args.checkm2_args) if args.checkm2_args else []

    if args.database:
        ok, msg = check_db_path(args.database, 'checkm2')
        if not ok:
            log_error(msg)
            sys.exit(1)
        if msg:
            log_warning(msg)
        log_info(f"Using database: {args.database}")

    print_section("Processing Genomes")

    results = {}

    if args.selected_dir:
        validate_dir_exists(args.selected_dir, "Selected directory")

        tiers = [t.strip() for t in args.tiers.split(',')]
        log_info(f"Processing tiers: {', '.join(tiers)}")

        for tier in tiers:
            results[tier] = process_tier_directory(
                base_dir=args.selected_dir,
                tier=tier,
                output_base=args.output,
                extension=args.extension,
                threads=args.threads,
                database=args.database,
                force=args.force,
                extra_args=extra_args,
            )

    else:
        validate_dir_exists(args.genome_dir, "Genome directory")

        genome_files = list(args.genome_dir.glob(f"*.{args.extension}"))
        n_genomes = len(genome_files)

        if n_genomes == 0:
            log_error(f"No genomes found with extension .{args.extension}")
            sys.exit(1)

        log_info(f"Found {n_genomes:,} genome files")

        args.output.mkdir(parents=True, exist_ok=True)

        exit_code = run_checkm2_predict(
            genome_dir=args.genome_dir,
            output_dir=args.output,
            extension=args.extension,
            threads=args.threads,
            database=args.database,
            force=args.force,
            extra_args=extra_args,
        )

        results['single'] = {
            'status': 'completed' if exit_code == 0 else 'failed',
            'exit_code': exit_code,
            'n_genomes': n_genomes,
        }

    print_section("SUMMARY")

    total_genomes = 0
    failed = 0

    for tier, result in results.items():
        status = result['status']
        n_genomes = result.get('n_genomes', 0)
        total_genomes += n_genomes

        if status == 'completed':
            print(f"  {tier:3s}: ✓ Completed ({n_genomes:,} genomes)")
        elif status == 'failed':
            failed += 1
            print(f"  {tier:3s}: ✗ Failed ({n_genomes:,} genomes)")
        elif status == 'skipped':
            reason = result.get('reason', 'unknown')
            print(f"  {tier:3s}: - Skipped ({reason})")

    print(f"\n  Total genomes processed: {total_genomes:,}")

    if failed > 0:
        log_error(f"\n{failed} tier(s) failed")
        sys.exit(1)

    log_info("\n✓ Done!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    register_parser(parser.add_subparsers(dest='command'))
    args = parser.parse_args()
    run(args)
