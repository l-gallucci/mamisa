#!/usr/bin/env python3
"""
MaMISA - setup-workflow command

One command to prepare the Snakemake workflow:
  1. points each workflow/envs/*.yaml at this local checkout
     (replaces the `- mamisa` pip line with `- -e <repo>`),
  2. writes the paths/threads into workflow/config.yaml, and
  3. optionally builds every per-rule conda env (`--install`).

Pure file edits otherwise — no YAML dependency, comments in config.yaml are kept.
"""

import re
import sys
import shutil
import argparse
import subprocess
from pathlib import Path

from .. import PACKAGE_ROOT
from ..utils.validation import check_db_path
from ..utils.logging import log_info, log_warning, log_error, print_header, print_section


def default_repo() -> Path:
    """The repo root = parent of the installed `mamisa` package directory."""
    return PACKAGE_ROOT.parent


def set_yaml_scalar(text: str, key: str, value, quote: bool = False) -> str:
    """
    Replace (or append) a top-level `key: value` line in a YAML text, preserving
    comments and layout. Only touches a line that starts at column 0.
    """
    rendered = f'"{value}"' if quote else str(value)
    pattern = re.compile(rf'^{re.escape(key)}:.*$', re.MULTILINE)
    if pattern.search(text):
        return pattern.sub(f'{key}: {rendered}', text)
    # append if missing
    sep = '' if text.endswith('\n') else '\n'
    return f'{text}{sep}{key}: {rendered}\n'


def patch_env_files(envs_dir: Path, repo: Path, dry_run: bool) -> int:
    """Replace the `- mamisa` pip entry with `- -e <repo>` in every env yaml."""
    env_files = sorted(envs_dir.glob('*.yaml'))
    if not env_files:
        log_warning(f"No env files found in {envs_dir}")
        return 0

    # matches a list item that is exactly `mamisa` or an existing `-e <path>`
    entry_re = re.compile(r'^(?P<indent>\s*-\s+)(mamisa|-e\s+\S+)\s*$', re.MULTILINE)
    patched = 0
    for env_file in env_files:
        text = env_file.read_text()
        new_text, n = entry_re.subn(lambda m: f"{m.group('indent')}-e {repo}", text)
        if n == 0:
            log_warning(f"  {env_file.name}: no `- mamisa` pip entry found, left as-is")
            continue
        if dry_run:
            log_info(f"  Would patch {env_file.name}: pip -> -e {repo}")
        else:
            env_file.write_text(new_text)
            log_info(f"  Patched {env_file.name}: pip -> -e {repo}")
        patched += 1
    return patched


def register_parser(subparsers):
    parser = subparsers.add_parser(
        'setup-workflow',
        help='Configure the Snakemake workflow (env yamls + config.yaml) and optionally build envs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Point the workflow at this checkout and fill config paths
  mamisa setup-workflow \\
    --genomes-dir bins/ --genome-ext fa --threads 40 \\
    --gtdbtk-data /data/gtdbtk_r220 \\
    --gunc-db /data/gunc/gunc_db_progenomes2.1.dmnd

  # Same, then build every per-rule conda env (needs snakemake + mamba/conda)
  mamisa setup-workflow --genomes-dir bins/ --install
"""
    )

    parser.add_argument('--repo', type=Path, default=None,
                        help='Path to the MaMISA checkout (default: auto-detected)')
    parser.add_argument('--workflow-dir', type=Path, default=None,
                        help='Workflow directory (default: <repo>/workflow)')

    # config values (only those given are written)
    parser.add_argument('--genomes-dir', type=Path)
    parser.add_argument('--genome-ext', default=None)
    parser.add_argument('--outdir', type=Path, default=None)
    parser.add_argument('--threads', type=int, default=None)
    parser.add_argument('--gtdbtk-data', type=Path, default=None,
                        help='GTDB-Tk data DIRECTORY (GTDBTK_DATA_PATH)')
    parser.add_argument('--gunc-db', type=Path, default=None,
                        help='GUNC database FILE (.dmnd)')
    parser.add_argument('--checkm2-db', type=Path, default=None,
                        help='CheckM2 database FILE (.dmnd)')
    parser.add_argument('--sample-regex', default=None)
    parser.add_argument('--tax-level', default=None,
                        choices=['domain', 'phylum', 'class', 'order',
                                 'family', 'genus', 'species'])

    # existing conda env names (enables the named-env run mode)
    parser.add_argument('--mamisa-env', default=None,
                        help='Existing env with the light CLI tools + MaMISA')
    parser.add_argument('--checkm2-env', default=None, help='Existing env with CheckM2')
    parser.add_argument('--gtdbtk-env', default=None, help='Existing env with GTDB-Tk')
    parser.add_argument('--gunc-env', default=None, help='Existing env with GUNC')

    parser.add_argument('--install', action='store_true',
                        help='Build all per-rule conda envs via snakemake --conda-create-envs-only')
    parser.add_argument('--no-check', action='store_true',
                        help='Skip the env sanity check when env names are given')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would change without writing')

    parser.set_defaults(func=run)
    return parser


def run(args):
    print_header("MaMISA - Setup Workflow")

    repo = (args.repo or default_repo()).resolve()
    workflow_dir = (args.workflow_dir or (repo / 'workflow')).resolve()
    envs_dir = workflow_dir / 'envs'
    config_file = workflow_dir / 'config.yaml'

    if not workflow_dir.is_dir():
        log_error(f"Workflow directory not found: {workflow_dir}")
        log_error("Pass --repo /path/to/mamisa or --workflow-dir explicitly.")
        sys.exit(1)
    if not config_file.exists():
        log_error(f"config.yaml not found: {config_file}")
        sys.exit(1)

    log_info(f"Repo:        {repo}")
    log_info(f"Workflow:    {workflow_dir}")

    # 1. patch env files to install this checkout
    print_section("Patching env files")
    patch_env_files(envs_dir, repo, args.dry_run)

    # 2. write config values
    print_section("Updating config.yaml")
    text = config_file.read_text()
    updates = []
    if args.genomes_dir is not None:
        text = set_yaml_scalar(text, 'genomes_dir', args.genomes_dir, quote=True); updates.append('genomes_dir')
    if args.genome_ext is not None:
        text = set_yaml_scalar(text, 'genome_ext', args.genome_ext, quote=True); updates.append('genome_ext')
    if args.outdir is not None:
        text = set_yaml_scalar(text, 'outdir', args.outdir, quote=True); updates.append('outdir')
    if args.threads is not None:
        text = set_yaml_scalar(text, 'threads', args.threads); updates.append('threads')
    # database paths — validate KIND (dir vs file) per tool before recording
    for flag_val, cfg_key, tool in (
        (args.gtdbtk_data, 'gtdbtk_data', 'gtdbtk'),
        (args.gunc_db, 'gunc_db', 'gunc'),
        (args.checkm2_db, 'checkm2_db', 'checkm2'),
    ):
        if flag_val is None:
            continue
        ok, msg = check_db_path(flag_val, tool)
        if not ok:
            log_error(f"  {msg}")
            sys.exit(1)
        if msg:
            log_warning(f"  {msg}")
        text = set_yaml_scalar(text, cfg_key, flag_val, quote=True)
        updates.append(cfg_key)

    if args.sample_regex is not None:
        text = set_yaml_scalar(text, 'sample_regex', args.sample_regex, quote=True); updates.append('sample_regex')
    if args.tax_level is not None:
        text = set_yaml_scalar(text, 'tax_level', args.tax_level, quote=True); updates.append('tax_level')

    # existing env names — recording them switches the workflow to named-env mode
    env_map = {
        'env_mamisa': args.mamisa_env,
        'env_checkm2': args.checkm2_env,
        'env_gtdbtk': args.gtdbtk_env,
        'env_gunc': args.gunc_env,
    }
    any_env = any(v for v in env_map.values())
    for key, val in env_map.items():
        if val:
            text = set_yaml_scalar(text, key, val, quote=True); updates.append(key)
    if any_env:
        text = set_yaml_scalar(text, 'use_named_envs', 'true'); updates.append('use_named_envs')

    if updates:
        if args.dry_run:
            log_info(f"  Would set: {', '.join(updates)}")
        else:
            config_file.write_text(text)
            log_info(f"  Set: {', '.join(updates)}")
    else:
        log_info("  No config values passed; left unchanged")

    # 3. sanity-check any env names the user supplied
    provided_envs = {
        'mamisa': args.mamisa_env,
        'checkm2': args.checkm2_env,
        'gtdbtk': args.gtdbtk_env,
        'gunc': args.gunc_env,
    }
    provided_envs = {r: e for r, e in provided_envs.items() if e}
    if provided_envs and not args.no_check and not args.dry_run:
        print_section("Sanity-checking supplied envs")
        if shutil.which('conda') is None:
            log_warning("conda not found; skipping env sanity check")
        else:
            from .check_envs import run_checks
            if not run_checks(provided_envs, repo):
                log_error("\nOne or more envs are not compatible — see fix hints above.")
                log_error("Re-run after fixing, or pass --no-check to skip this gate.")
                sys.exit(1)

    # 4. optionally build the conda envs
    if args.install:
        print_section("Building per-rule conda envs")
        if args.dry_run:
            log_info("  Would run: snakemake --use-conda --conda-create-envs-only ...")
        elif shutil.which('snakemake') is None:
            log_error("snakemake not found in PATH.")
            log_error("Install it: conda install -c conda-forge -c bioconda snakemake-minimal")
            sys.exit(1)
        else:
            cmd = ['snakemake', '--use-conda', '--conda-create-envs-only',
                   '--cores', str(args.threads or 1),
                   '-s', str(workflow_dir / 'Snakefile'),
                   '--configfile', str(config_file)]
            log_info(f"  {' '.join(cmd)}")
            rc = subprocess.run(cmd).returncode
            if rc != 0:
                log_error(f"Env creation failed (exit {rc})")
                sys.exit(rc)

    log_info("\n✓ Done! Next: mamisa fetch-databases ... then "
             "snakemake --use-conda --cores N -s workflow/Snakefile --configfile workflow/config.yaml")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    register_parser(parser.add_subparsers(dest='command'))
    args = parser.parse_args()
    run(args)
