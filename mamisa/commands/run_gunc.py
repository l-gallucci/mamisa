#!/usr/bin/env python3
"""
MaMISA - run-gunc command

Wrapper for GUNC (Genome UNClutterer) chimerism/contamination detection.

GUNC detects chimeric and contaminated genomes by testing taxonomic consistency
of genes across the genome — orthogonal to MaMISA's GC- and read-based signals.
This fills the gap where genes from different clades are joined into one MAG.

NOTE ON DEPENDENCIES
--------------------
GUNC 1.1.1 is an old (~2021) pure-python package whose conda recipe declares
`pandas>=2.0.0` with NO upper bound, so a naive install pulls pandas 3.x and a
bleeding-edge Python that break it at runtime. A verified-working environment:

  mamba create -n gunc -c conda-forge -c bioconda \\
      "gunc=1.1.1" "python=3.10" "pandas>=2,<3" "numpy<2" diamond prodigal

GUNC also needs its reference database (~13 GB), set via --db-file or the
GUNC_DB environment variable (`gunc download_db`).

This command runs an env-health check before launching and reports clearly if
the installed GUNC is broken by an incompatible pandas/numpy/python.
"""

import os
import sys
import shlex
import subprocess
import argparse
from pathlib import Path
from typing import Optional

from ..utils.validation import validate_dir_exists, check_dependencies, check_db_path
from ..utils.logging import log_info, log_error, log_warning, print_header, print_section


RECOMMENDED_ENV = (
    'mamba create -n gunc -c conda-forge -c bioconda '
    '"gunc=1.1.1" "python=3.10" "pandas>=2,<3" "numpy<2" diamond prodigal'
)


def check_gunc_env() -> bool:
    """
    Verify gunc + its backends are present AND that the installed gunc actually
    runs (i.e. its pandas/numpy/python stack is not broken). Returns True if OK.
    """
    deps = check_dependencies(['gunc', 'diamond', 'prodigal'])
    missing = [name for name, path in deps.items() if path is None]
    if missing:
        log_error(f"Missing dependencies: {', '.join(missing)}")
        log_error("Install a known-good GUNC environment with:")
        log_error(f"  {RECOMMENDED_ENV}")
        return False

    for name, path in deps.items():
        log_info(f"Found {name}: {path}")

    # Env-health: `gunc --version` imports pandas/numpy. A broken stack surfaces
    # here as a traceback / non-zero exit rather than a version string.
    try:
        result = subprocess.run(['gunc', '--version'],
                                capture_output=True, text=True,
                                stdin=subprocess.DEVNULL, timeout=60)
    except Exception as e:
        log_error(f"Could not run 'gunc --version': {e}")
        return False

    combined = (result.stdout + result.stderr)
    broken_markers = ('Traceback', 'ImportError', 'ModuleNotFoundError',
                      'AttributeError', 'cannot import name')
    if result.returncode != 0 or any(m in combined for m in broken_markers):
        log_error("GUNC is installed but fails to run — its Python/pandas/numpy "
                  "stack is likely incompatible (GUNC is old).")
        log_error("Captured output:")
        for line in combined.strip().splitlines()[:8]:
            log_error(f"  {line}")
        log_error("Rebuild the environment with pinned versions:")
        log_error(f"  {RECOMMENDED_ENV}")
        return False

    version = result.stdout.strip() or result.stderr.strip()
    log_info(f"GUNC version: {version}")
    return True


def resolve_db(db_file: Optional[Path]) -> Optional[Path]:
    """Resolve the GUNC database from --db-file or the GUNC_DB env var."""
    if db_file:
        return db_file
    env_db = os.environ.get('GUNC_DB')
    if env_db:
        return Path(env_db)
    return None


def run_gunc(input_dir: Path, output_dir: Path, file_suffix: str,
             threads: int, db_file: Optional[Path],
             extra_args: list = None) -> int:
    """Run `gunc run` on a directory of genomes. Returns exit code."""

    cmd = [
        'gunc', 'run',
        '--input_dir', str(input_dir),
        '--out_dir', str(output_dir),
        '--file_suffix', file_suffix,
        '--threads', str(threads),
    ]
    if db_file:
        cmd.extend(['--db_file', str(db_file)])
    if extra_args:
        cmd.extend(extra_args)

    log_info("\nRunning command:")
    log_info(f"  {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, check=True, stdin=subprocess.DEVNULL)
        return result.returncode
    except subprocess.CalledProcessError as e:
        log_error(f"GUNC failed with exit code {e.returncode}")
        return e.returncode
    except KeyboardInterrupt:
        log_warning("\nInterrupted by user")
        return 130


def process_tier_directory(base_dir: Path, tier: str, output_base: Path,
                           file_suffix: str, threads: int,
                           db_file: Optional[Path],
                           extra_args: list = None) -> dict:
    """Run GUNC on a single tier subdirectory (HQ/MQ/LQ)."""
    tier_input = base_dir / tier
    tier_output = output_base / tier

    if not tier_input.exists():
        log_warning(f"Tier directory not found: {tier_input}")
        return {'status': 'skipped', 'reason': 'directory_not_found'}

    genomes = list(tier_input.glob(f"*{file_suffix}"))
    if not genomes:
        log_warning(f"No genomes in {tier_input} with suffix {file_suffix}")
        return {'status': 'skipped', 'reason': 'no_genomes', 'n_genomes': 0}

    log_info(f"\nProcessing tier: {tier}  ({len(genomes):,} genomes)")
    tier_output.mkdir(parents=True, exist_ok=True)

    exit_code = run_gunc(tier_input, tier_output, file_suffix, threads,
                         db_file, extra_args)
    return {
        'status': 'completed' if exit_code == 0 else 'failed',
        'exit_code': exit_code,
        'n_genomes': len(genomes),
    }


def register_parser(subparsers):
    parser = subparsers.add_parser(
        'run-gunc',
        help='Run GUNC gene-level chimerism/contamination detection',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
GUNC detects chimeras by taxonomic inconsistency of genes across a genome,
complementing check-chimeras (GC/read signals). Feed the output back in with
`mamisa check-chimeras --gunc-dir <out>`.

Dependencies (GUNC is old — pin versions):
  {RECOMMENDED_ENV}
  gunc download_db -db gtdb_214 ./gunc_db/   # GUNC>=1.1.1, from Zenodo; or set GUNC_DB
  # DB choices: progenomes_2.1 (default), progenomes_3, gtdb_95, gtdb_214
  # gtdb_214 recommended for environmental MAGs

Examples:
  # Single genome directory
  mamisa run-gunc \\
    --genome-dir genomes/ --output gunc_out/ \\
    --file-suffix .fa --threads 20 --db-file gunc_db/gunc_db_gtdb214.dmnd

  # All tiers from organize-mags, DB from $GUNC_DB
  mamisa run-gunc \\
    --selected-dir filtered/Selected/ --output gunc_out/ \\
    --file-suffix .fa --threads 20
        """
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--selected-dir', type=Path,
                             help='Directory with tier subdirectories (HQ/MQ/LQ)')
    input_group.add_argument('--genome-dir', type=Path,
                             help='Single directory with genome files')

    parser.add_argument('-o', '--output', type=Path, required=True,
                        help='Output directory for GUNC results')
    parser.add_argument('--file-suffix', default='.fa',
                        help='Genome file suffix (default: .fa)')
    parser.add_argument('--threads', type=int, default=1,
                        help='Number of threads (default: 1)')
    parser.add_argument('--db-file', type=Path,
                        help='GUNC diamond database (.dmnd); else uses $GUNC_DB')
    parser.add_argument('--tiers', default='HQ,MQ,LQ',
                        help='Comma-separated tiers to process (default: HQ,MQ,LQ)')
    parser.add_argument('--gunc-args', default='',
                        help='Additional arguments to pass to gunc run (quoted)')

    parser.set_defaults(func=run)
    return parser


def run(args):
    print_header("MaMISA - Run GUNC")

    print_section("Checking Dependencies")
    if not check_gunc_env():
        sys.exit(1)

    db_file = resolve_db(args.db_file)
    if db_file is None:
        log_warning("No --db-file and GUNC_DB not set; GUNC will error unless its "
                    "default DB location is configured.")
    else:
        ok, msg = check_db_path(db_file, 'gunc')
        if not ok:
            log_error(msg)
            sys.exit(1)
        if msg:
            log_warning(msg)
        log_info(f"Using GUNC database: {db_file}")

    extra_args = shlex.split(args.gunc_args) if args.gunc_args else []

    print_section("Running GUNC")
    results = {}

    if args.selected_dir:
        validate_dir_exists(args.selected_dir, "Selected directory")
        tiers = [t.strip() for t in args.tiers.split(',')]
        for tier in tiers:
            results[tier] = process_tier_directory(
                args.selected_dir, tier, args.output,
                args.file_suffix, args.threads, db_file, extra_args)
    else:
        validate_dir_exists(args.genome_dir, "Genome directory")
        genomes = list(args.genome_dir.glob(f"*{args.file_suffix}"))
        if not genomes:
            log_error(f"No genomes with suffix {args.file_suffix}")
            sys.exit(1)
        log_info(f"Found {len(genomes):,} genomes")
        args.output.mkdir(parents=True, exist_ok=True)
        exit_code = run_gunc(args.genome_dir, args.output, args.file_suffix,
                             args.threads, db_file, extra_args)
        results['single'] = {
            'status': 'completed' if exit_code == 0 else 'failed',
            'exit_code': exit_code, 'n_genomes': len(genomes),
        }

    print_section("SUMMARY")
    failed = 0
    for tier, result in results.items():
        status = result['status']
        n = result.get('n_genomes', 0)
        if status == 'completed':
            print(f"  {tier:3s}: ✓ Completed ({n:,} genomes)")
        elif status == 'failed':
            failed += 1
            print(f"  {tier:3s}: ✗ Failed ({n:,} genomes)")
        else:
            print(f"  {tier:3s}: - Skipped ({result.get('reason','unknown')})")

    if failed:
        log_error(f"\n{failed} run(s) failed")
        sys.exit(1)
    log_info("\n✓ Done! Feed results to check-chimeras with --gunc-dir")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    register_parser(parser.add_subparsers(dest='command'))
    args = parser.parse_args()
    run(args)
