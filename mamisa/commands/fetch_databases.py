#!/usr/bin/env python3
"""
MaMISA - fetch-databases command

One command to get the reference databases the workflow needs:
  - if you ALREADY have a database, point to it with --gtdbtk-data / --gunc-db /
    --checkm2-db and it is validated and written into workflow/config.yaml;
  - otherwise pass --download and the matching tool's own downloader is run
    (inside its conda env) into --db-dir, then the path is written to config.

Databases handled: GTDB-Tk data, GUNC DB (~13 GB), CheckM2 DIAMOND DB.
"""

import os
import sys
import shutil
import argparse
import subprocess
from pathlib import Path

from .. import PACKAGE_ROOT
from .setup_workflow import set_yaml_scalar, default_repo
from ..utils.validation import check_db_path
from ..utils.logging import log_info, log_warning, log_error, print_header, print_section


def conda_run(env: str, cmd: list, extra_env: dict = None) -> int:
    """Run a command inside a named conda env via `conda run`."""
    full = ['conda', 'run', '-n', env] + cmd
    log_info(f"  {' '.join(str(c) for c in full)}")
    run_env = dict(os.environ)
    if extra_env:
        run_env.update({k: str(v) for k, v in extra_env.items()})
    try:
        return subprocess.run(full, env=run_env).returncode
    except FileNotFoundError:
        log_error("conda not found in PATH")
        return 127


def register_parser(subparsers):
    parser = subparsers.add_parser(
        'fetch-databases',
        help='Download reference databases, or register existing ones into config.yaml',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # I already have the databases — just record their locations
  mamisa fetch-databases \\
    --gtdbtk-data /data/gtdbtk_r220 \\
    --gunc-db /data/gunc/gunc_db_progenomes2.1.dmnd \\
    --checkm2-db /data/checkm2/uniref100.KO.1.dmnd

  # Download the ones I'm missing into /data/mamisa_dbs
  mamisa fetch-databases --download --db-dir /data/mamisa_dbs --tools gunc,checkm2

Downloads run inside the tool's conda env (override names with --*-env).
"""
    )

    parser.add_argument('--config', type=Path, default=None,
                        help='config.yaml to update (default: <repo>/workflow/config.yaml)')
    parser.add_argument('--db-dir', type=Path, default=Path('databases'),
                        help='Directory for downloaded databases (default: ./databases)')
    parser.add_argument('--tools', default='gtdbtk,gunc,checkm2',
                        help='Which databases to handle (default: gtdbtk,gunc,checkm2)')
    parser.add_argument('--download', action='store_true',
                        help='Download any database not supplied as an existing path')

    # existing database locations (path KIND is checked per tool)
    parser.add_argument('--gtdbtk-data', type=Path,
                        help='Existing GTDB-Tk data DIRECTORY (GTDBTK_DATA_PATH)')
    parser.add_argument('--gunc-db', type=Path,
                        help='Existing GUNC database FILE (.dmnd)')
    parser.add_argument('--gunc-db-source',
                        choices=['progenomes_2.1', 'progenomes_3',
                                 'gtdb_95', 'gtdb_214'],
                        default='progenomes_2.1',
                        help='Which GUNC reference DB to download (GUNC >=1.1.1 '
                             'names; default: progenomes_2.1). For environmental '
                             'MAGs prefer gtdb_214 (GTDB r214): broader taxa, '
                             'fewer false chimeric calls on novel lineages')
    parser.add_argument('--checkm2-db', type=Path,
                        help='Existing CheckM2 database FILE (.dmnd)')

    # conda env names used for downloads
    parser.add_argument('--gtdbtk-env', default='gtdbtk')
    parser.add_argument('--gunc-env', default='gunc')
    parser.add_argument('--checkm2-env', default='checkm2')

    parser.add_argument('--dry-run', action='store_true')

    parser.set_defaults(func=run)
    return parser


def _write_config(config_file: Path, key: str, value, dry_run: bool):
    if config_file is None or not config_file.exists():
        log_warning(f"  config not found; not recording {key}={value}")
        return
    if dry_run:
        log_info(f"  Would set {key} = {value} in {config_file.name}")
        return
    text = set_yaml_scalar(config_file.read_text(), key, value, quote=True)
    config_file.write_text(text)
    log_info(f"  Recorded {key} in {config_file.name}")


def handle_gtdbtk(args, config_file):
    print_section("GTDB-Tk data")
    if args.gtdbtk_data:
        ok, msg = check_db_path(args.gtdbtk_data, 'gtdbtk')
        if not ok:
            log_error(f"  {msg}")
            return False
        if msg:
            log_warning(f"  {msg}")
        log_info(f"  Using existing: {args.gtdbtk_data}")
        _write_config(config_file, 'gtdbtk_data', args.gtdbtk_data.resolve(), args.dry_run)
        return True
    if not args.download:
        log_info("  No path given and --download not set; skipped")
        return True
    target = (args.db_dir / 'gtdbtk').resolve()
    log_info(f"  Downloading GTDB-Tk data into {target} (large; this takes a while)")
    if args.dry_run:
        log_info(f"  Would run: conda run -n {args.gtdbtk_env} download-db.sh "
                 f"(GTDBTK_DATA_PATH={target})")
        return True
    target.mkdir(parents=True, exist_ok=True)
    rc = conda_run(args.gtdbtk_env, ['download-db.sh'], extra_env={'GTDBTK_DATA_PATH': target})
    if rc != 0:
        log_error("  GTDB-Tk download failed. download-db.sh ships with gtdbtk; "
                  "check the env or download the release tarball manually.")
        return False
    _write_config(config_file, 'gtdbtk_data', target, args.dry_run)
    return True


def handle_gunc(args, config_file):
    print_section("GUNC database")
    if args.gunc_db:
        ok, msg = check_db_path(args.gunc_db, 'gunc')
        if not ok:
            log_error(f"  {msg}")
            return False
        if msg:
            log_warning(f"  {msg}")
        log_info(f"  Using existing: {args.gunc_db}")
        _write_config(config_file, 'gunc_db', args.gunc_db.resolve(), args.dry_run)
        return True
    if not args.download:
        log_info("  No path given and --download not set; skipped")
        return True
    source = getattr(args, 'gunc_db_source', 'progenomes_2.1')
    target = (args.db_dir / 'gunc').resolve()
    log_info(f"  Downloading GUNC DB ({source}) into {target} (~13-16 GB)")
    # GUNC >=1.1.1 fetches from Zenodo (old EMBL webserver is a fallback).
    dl_cmd = ['gunc', 'download_db', '-db', source, str(target)]
    if args.dry_run:
        log_info(f"  Would run: conda run -n {args.gunc_env} {' '.join(dl_cmd)}")
        return True
    target.mkdir(parents=True, exist_ok=True)
    rc = conda_run(args.gunc_env, dl_cmd)
    if rc != 0:
        log_error("  GUNC download failed.")
        return False
    # Prefer the .dmnd matching the requested source if several are present
    # (filenames like gunc_db_gtdb214.dmnd / gunc_db_progenomes2.1.dmnd).
    dmnds = sorted(target.glob('*.dmnd'))
    if dmnds:
        tag = source.replace('_', '').lower()
        chosen = next((d for d in dmnds if tag in d.name.lower()), dmnds[0])
        _write_config(config_file, 'gunc_db', chosen.resolve(), args.dry_run)
    else:
        log_warning(f"  No .dmnd found under {target}; set gunc_db manually")
    return True


def handle_checkm2(args, config_file):
    print_section("CheckM2 database")
    if args.checkm2_db:
        ok, msg = check_db_path(args.checkm2_db, 'checkm2')
        if not ok:
            log_error(f"  {msg}")
            return False
        if msg:
            log_warning(f"  {msg}")
        log_info(f"  Using existing: {args.checkm2_db}")
        _write_config(config_file, 'checkm2_db', args.checkm2_db.resolve(), args.dry_run)
        return True
    if not args.download:
        log_info("  No path given and --download not set; skipped")
        return True
    target = (args.db_dir / 'checkm2').resolve()
    log_info(f"  Downloading CheckM2 DB into {target}")
    if args.dry_run:
        log_info(f"  Would run: conda run -n {args.checkm2_env} checkm2 database "
                 f"--download --path {target}")
        return True
    target.mkdir(parents=True, exist_ok=True)
    rc = conda_run(args.checkm2_env, ['checkm2', 'database', '--download', '--path', str(target)])
    if rc != 0:
        log_error("  CheckM2 download failed.")
        return False
    dmnds = sorted(target.rglob('*.dmnd'))
    if dmnds:
        _write_config(config_file, 'checkm2_db', dmnds[0].resolve(), args.dry_run)
    return True


def run(args):
    print_header("MaMISA - Fetch Databases")

    config_file = args.config
    if config_file is None:
        config_file = default_repo() / 'workflow' / 'config.yaml'
    if not config_file.exists():
        log_warning(f"config.yaml not found at {config_file}; paths will not be recorded")
        config_file = None
    else:
        log_info(f"Config: {config_file}")

    if args.download and shutil.which('conda') is None:
        log_error("conda not found in PATH — required for --download")
        sys.exit(1)

    tools = {t.strip().lower() for t in args.tools.split(',') if t.strip()}
    handlers = {
        'gtdbtk': handle_gtdbtk,
        'gunc': handle_gunc,
        'checkm2': handle_checkm2,
    }

    ok = True
    for name in ('gtdbtk', 'gunc', 'checkm2'):
        if name in tools:
            ok = handlers[name](args, config_file) and ok

    print_section("SUMMARY")
    if args.dry_run:
        log_warning("DRY-RUN: nothing downloaded or written")
    if not ok:
        log_error("One or more databases failed; see messages above")
        sys.exit(1)
    log_info("✓ Done!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    register_parser(parser.add_subparsers(dest='command'))
    args = parser.parse_args()
    run(args)
