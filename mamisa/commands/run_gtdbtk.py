#!/usr/bin/env python3
"""
MaMISA - run-gtdbtk command
Wrapper for running GTDB-Tk taxonomy classification
"""

import sys
import shlex
import subprocess
import argparse
import shutil
from pathlib import Path
from typing import Optional

from ..utils.validation import validate_dir_exists, check_dependencies
from ..utils.logging import log_info, log_error, log_warning, print_header, print_section


def check_gtdbtk_available() -> bool:
    """Check if gtdbtk is available in PATH"""
    deps = check_dependencies(['gtdbtk'])
    
    if deps['gtdbtk'] is None:
        log_error("gtdbtk not found in PATH")
        log_error("Please install GTDB-Tk: https://ecogenomics.github.io/GTDBTk/")
        return False
    
    log_info(f"Found gtdbtk: {deps['gtdbtk']}")
    return True


def get_gtdbtk_version() -> Optional[str]:
    """Get GTDB-Tk version"""
    try:
        result = subprocess.run(['gtdbtk', '--version'],
                                capture_output=True, text=True,
                                stdin=subprocess.DEVNULL, timeout=120)
        return (result.stdout or result.stderr).strip() or None
    except Exception:
        return None


def run_gtdbtk_classify(genome_dir: Path, output_dir: Path,
                       extension: str, cpus: int,
                       place_species: bool = False,
                       extra_args: list = None) -> int:
    """
    Run GTDB-Tk classify_wf

    Returns:
        Exit code from gtdbtk
    """

    cmd = [
        'gtdbtk',
        'classify_wf',
        '--genome_dir', str(genome_dir),
        '--out_dir', str(output_dir),
        '--cpus', str(cpus),
        '--extension', extension
    ]

    # --place_species: force pplacer tree placement even when skani classifies.
    # (Replaces the deprecated --skip_ani_screen; GTDB-Tk >=2.7.)
    if place_species:
        cmd.append('--place_species')

    if extra_args:
        cmd.extend(extra_args)
    
    log_info(f"\nRunning command:")
    log_info(f"  {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, stdin=subprocess.DEVNULL)
        return result.returncode
    except subprocess.CalledProcessError as e:
        log_error(f"GTDB-Tk failed with exit code {e.returncode}")
        return e.returncode
    except KeyboardInterrupt:
        log_warning("\nInterrupted by user")
        return 130


def process_tier_directory(base_dir: Path, tier: str, output_base: Path,
                          extension: str, cpus: int,
                          place_species: bool = False,
                          extra_args: list = None) -> dict:
    """
    Process a single tier directory
    """
    
    tier_genome_dir = base_dir / tier
    tier_output_dir = output_base / tier
    
    if not tier_genome_dir.exists():
        log_warning(f"Tier directory not found: {tier_genome_dir}")
        return {'status': 'skipped', 'reason': 'directory_not_found'}
    
    # Count genomes
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
    
    exit_code = run_gtdbtk_classify(
        genome_dir=tier_genome_dir,
        output_dir=tier_output_dir,
        extension=extension,
        cpus=cpus,
        place_species=place_species,
        extra_args=extra_args
    )
    
    return {
        'status': 'completed' if exit_code == 0 else 'failed',
        'exit_code': exit_code,
        'n_genomes': n_genomes
    }


def register_parser(subparsers):
    """Register this command's parser"""
    parser = subparsers.add_parser(
        'run-gtdbtk',
        help='Run GTDB-Tk taxonomy classification',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run on all tiers in Selected directory
  mamisa run-gtdbtk \\
    --selected-dir filtered/Selected/ \\
    --output gtdbtk_results/ \\
    --extension fa \\
    --cpus 40
  
  # Run on specific tier only
  mamisa run-gtdbtk \\
    --genome-dir filtered/Selected/HQ/ \\
    --output gtdbtk_results/HQ/ \\
    --extension fa \\
    --cpus 40
  
  # Force species placement in the pplacer tree (GTDB-Tk >=2.7)
  mamisa run-gtdbtk \\
    --selected-dir filtered/Selected/ \\
    --output gtdbtk_results/ \\
    --extension fa \\
    --cpus 40 \\
    --place-species

Note: GTDB-Tk removed Mash (--mash_db) in v2.5.0 and now screens with skani.
The old --mash-db flag is accepted but ignored.
        """
    )
    
    # Input options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--selected-dir', type=Path,
                            help='Directory with tier subdirectories (HQ/MQ/LQ)')
    input_group.add_argument('--genome-dir', type=Path,
                            help='Single directory with genome files')
    
    parser.add_argument('-o', '--output', type=Path, required=True,
                        help='Output directory for GTDB-Tk results')
    
    # GTDB-Tk options
    parser.add_argument('--extension', default='fa',
                        help='Genome file extension (default: fa)')
    parser.add_argument('--cpus', type=int, default=1,
                        help='Number of CPUs to use (default: 1)')
    parser.add_argument('--place-species', action='store_true',
                        help='Pass --place_species: place genomes in the pplacer '
                             'tree even when skani classifies them (GTDB-Tk >=2.7)')
    parser.add_argument('--mash-db', type=Path,
                        help='DEPRECATED and ignored: GTDB-Tk removed Mash in '
                             'v2.5.0 (now uses skani)')
    
    # Tier selection (for --selected-dir mode)
    parser.add_argument('--tiers', default='HQ,MQ,LQ',
                        help='Comma-separated tiers to process (default: HQ,MQ,LQ)')
    
    # Extra arguments to pass to gtdbtk
    parser.add_argument('--gtdbtk-args', default='',
                        help='Additional arguments to pass to gtdbtk (quoted string)')
    
    parser.set_defaults(func=run)
    return parser


def run(args):
    """Execute the run-gtdbtk command"""
    
    print_header("MaMISA - Run GTDB-Tk")
    
    # Check if gtdbtk is available
    print_section("Checking Dependencies")
    if not check_gtdbtk_available():
        sys.exit(1)
    
    version = get_gtdbtk_version()
    if version:
        log_info(f"GTDB-Tk version: {version}")
    
    # Parse extra arguments — use shlex.split so quoted values with spaces are handled
    extra_args = []
    if args.gtdbtk_args:
        extra_args = shlex.split(args.gtdbtk_args)
    
    # --mash-db is obsolete: GTDB-Tk dropped Mash in v2.5.0 (skani now)
    if args.mash_db:
        log_warning("--mash-db is ignored: GTDB-Tk removed the Mash screen in "
                    "v2.5.0 and now uses skani. Drop this flag.")
    
    # Process genomes
    print_section("Processing Genomes")
    
    results = {}
    
    if args.selected_dir:
        # Process multiple tiers
        validate_dir_exists(args.selected_dir, "Selected directory")
        
        tiers = [t.strip() for t in args.tiers.split(',')]
        log_info(f"Processing tiers: {', '.join(tiers)}")
        
        for tier in tiers:
            result = process_tier_directory(
                base_dir=args.selected_dir,
                tier=tier,
                output_base=args.output,
                extension=args.extension,
                cpus=args.cpus,
                place_species=args.place_species,
                extra_args=extra_args
            )
            results[tier] = result
    
    else:
        # Process single directory
        validate_dir_exists(args.genome_dir, "Genome directory")
        
        genome_files = list(args.genome_dir.glob(f"*.{args.extension}"))
        n_genomes = len(genome_files)
        
        if n_genomes == 0:
            log_error(f"No genomes found with extension .{args.extension}")
            sys.exit(1)
        
        log_info(f"Found {n_genomes:,} genome files")
        
        args.output.mkdir(parents=True, exist_ok=True)
        
        exit_code = run_gtdbtk_classify(
            genome_dir=args.genome_dir,
            output_dir=args.output,
            extension=args.extension,
            cpus=args.cpus,
            place_species=args.place_species,
            extra_args=extra_args
        )

        results['single'] = {
            'status': 'completed' if exit_code == 0 else 'failed',
            'exit_code': exit_code,
            'n_genomes': n_genomes
        }
    
    # Print summary
    print_section("SUMMARY")
    
    total_genomes = 0
    successful = 0
    failed = 0
    skipped = 0
    
    for tier, result in results.items():
        status = result['status']
        n_genomes = result.get('n_genomes', 0)
        total_genomes += n_genomes
        
        if status == 'completed':
            successful += 1
            print(f"  {tier:3s}: ✓ Completed ({n_genomes:,} genomes)")
        elif status == 'failed':
            failed += 1
            print(f"  {tier:3s}: ✗ Failed ({n_genomes:,} genomes)")
        elif status == 'skipped':
            skipped += 1
            reason = result.get('reason', 'unknown')
            print(f"  {tier:3s}: - Skipped ({reason})")
    
    print(f"\n  Total genomes processed: {total_genomes:,}")
    
    if failed > 0:
        log_error(f"\n{failed} tier(s) failed")
        sys.exit(1)
    
    log_info("\n✓ Done!")


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    register_parser(parser.add_subparsers(dest='command'))
    args = parser.parse_args()
    run(args)
