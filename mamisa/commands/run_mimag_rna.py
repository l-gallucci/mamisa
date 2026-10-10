#!/usr/bin/env python3
"""
MaMISA - run-mimag-rna command

Wrapper for the MIMAG rRNA/tRNA high-quality criterion. For each genome it runs:
  - barrnap       -> 5S / 16S / 23S rRNA genes
  - tRNAscan-SE   -> tRNA genes (counted as distinct standard amino acids)
and writes a single `mimag_rna_summary.tsv`. That summary can be fed back into
`organize-mags --mimag-rna-dir <out>` so the HQ tier enforces full MIMAG, not
just completeness/contamination.

CheckM2 does NOT evaluate rRNA/tRNA, so MaMISA's HQ tier is otherwise only an
"HQ by completeness/contamination" label. This command closes that gap.

Dependencies (one small conda env):
  mamba create -n mimag-rna -c conda-forge -c bioconda barrnap trnascan-se
"""

import sys
import shlex
import subprocess
import argparse
from pathlib import Path
from typing import Optional, List

from ..utils.validation import validate_dir_exists, check_dependencies
from ..utils.logging import log_info, log_error, log_warning, print_header, print_section
from ..utils.rrna import (
    parse_barrnap_gff,
    parse_trnascan,
    assess_mimag_rna,
    write_mimag_rna_summary,
)


RECOMMENDED_ENV = (
    'mamba create -n mimag-rna -c conda-forge -c bioconda barrnap trnascan-se'
)


def check_mimag_env() -> bool:
    """Verify barrnap + tRNAscan-SE are present."""
    deps = check_dependencies(['barrnap', 'tRNAscan-SE'])
    missing = [name for name, path in deps.items() if path is None]
    if missing:
        log_error(f"Missing dependencies: {', '.join(missing)}")
        log_error("Install with:")
        log_error(f"  {RECOMMENDED_ENV}")
        return False
    for name, path in deps.items():
        log_info(f"Found {name}: {path}")
    return True


def run_barrnap(genome: Path, kingdom: str, gff_out: Path,
                threads: int) -> bool:
    """Run barrnap on one genome, writing a GFF. Returns True on success."""
    cmd = ['barrnap', '--kingdom', kingdom, '--threads', str(threads),
           str(genome)]
    try:
        with open(gff_out, 'w') as out:
            subprocess.run(cmd, check=True, stdout=out,
                           stderr=subprocess.DEVNULL, stdin=subprocess.DEVNULL)
        return True
    except subprocess.CalledProcessError as e:
        log_error(f"barrnap failed on {genome.name} (exit {e.returncode})")
        return False


def run_trnascan(genome: Path, domain_flag: str, out_file: Path) -> bool:
    """Run tRNAscan-SE on one genome. Returns True on success."""
    # -B bacteria / -A archaea / -G general; -o tabular output; -q quiet
    cmd = ['tRNAscan-SE', domain_flag, '-q', '-o', str(out_file), str(genome)]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, stdin=subprocess.DEVNULL)
        return True
    except subprocess.CalledProcessError as e:
        log_error(f"tRNAscan-SE failed on {genome.name} (exit {e.returncode})")
        return False


_KINGDOM = {'bac': ('bac', '-B'), 'arc': ('arc', '-A'), 'euk': ('euk', '-E')}


def process_genome(genome: Path, work_dir: Path, kingdom: str,
                   threads: int) -> Optional[dict]:
    """Run both tools on one genome and return its MIMAG rRNA/tRNA record."""
    barrnap_kingdom, trnascan_flag = _KINGDOM.get(kingdom, ('bac', '-B'))

    gff_out = work_dir / f"{genome.stem}.barrnap.gff"
    trna_out = work_dir / f"{genome.stem}.trnascan.tsv"

    if not run_barrnap(genome, barrnap_kingdom, gff_out, threads):
        return None
    if not run_trnascan(genome, trnascan_flag, trna_out):
        return None

    rrna_counts = parse_barrnap_gff(gff_out)
    trna_aas = parse_trnascan(trna_out)
    verdict = assess_mimag_rna(rrna_counts, trna_aas)
    verdict['genome'] = genome.name
    return verdict


def collect_genomes(input_dir: Path, extension: str) -> List[Path]:
    return sorted(input_dir.glob(f"*.{extension.lstrip('.')}"))


_FASTA_EXTS = ('.fa', '.fasta', '.fna', '.fa.gz', '.fasta.gz', '.fna.gz')


def _strip_ext(name: str) -> str:
    low = name.lower()
    for ext in _FASTA_EXTS:
        if low.endswith(ext):
            return name[:-len(ext)]
    return name


def load_domains_from_gtdbtk(gtdbtk_dir: Path) -> dict:
    """
    Map genome id -> 'bac' / 'arc' using GTDB-Tk summaries.

    GTDB-Tk writes one summary per domain (gtdbtk.bac120.summary.tsv,
    gtdbtk.ar53/ar122.summary.tsv), so the domain is taken from the filename;
    the classification column (d__Archaea / d__Bacteria) is the fallback.
    Each genome is keyed by both its full name and extension-stripped stem.
    """
    import csv
    mapping: dict = {}
    summaries = sorted(gtdbtk_dir.rglob("gtdbtk.*.summary.tsv"))
    if not summaries:
        log_warning(f"No gtdbtk.*.summary.tsv found under {gtdbtk_dir}")
        return mapping
    for summary in summaries:
        fname = summary.name.lower()
        if 'bac120' in fname:
            file_domain = 'bac'
        elif 'ar53' in fname or 'ar122' in fname:
            file_domain = 'arc'
        else:
            file_domain = None
        try:
            with open(summary) as f:
                for row in csv.DictReader(f, delimiter='\t'):
                    g = (row.get('user_genome') or '').strip()
                    if not g:
                        continue
                    domain = file_domain
                    if domain is None:
                        cls = row.get('classification') or ''
                        domain = 'arc' if 'd__Archaea' in cls else 'bac'
                    mapping[g] = domain
                    mapping[_strip_ext(g)] = domain
        except OSError as e:
            log_warning(f"Could not read {summary}: {e}")
    n = len({k for k in mapping})
    log_info(f"Loaded domains for {n:,} genome keys from {len(summaries)} "
             f"GTDB-Tk summary file(s)")
    return mapping


def resolve_kingdom(genome: Path, domains: dict, default: str) -> str:
    """Resolve barrnap/tRNAscan kingdom for a genome from the domain map."""
    for key in (genome.name, genome.stem, _strip_ext(genome.name)):
        if key in domains:
            return domains[key]
    log_warning(f"  {genome.name}: domain not found in GTDB-Tk; "
                f"using default --kingdom {default}")
    return default


def register_parser(subparsers):
    parser = subparsers.add_parser(
        'run-mimag-rna',
        help='Check MIMAG rRNA (5S/16S/23S) + tRNA criterion via barrnap + tRNAscan-SE',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
Runs barrnap (rRNA) + tRNAscan-SE (tRNA) per genome and writes
mimag_rna_summary.tsv. Feed it to organize-mags to enforce full MIMAG HQ:

  mamisa organize-mags ... --mimag-rna-dir <this output>

Dependencies (one small env):
  {RECOMMENDED_ENV}

Examples:
  # Single genome directory (bacteria)
  mamisa run-mimag-rna \\
    --genome-dir genomes/ --output mimag_rna/ \\
    --extension fa --kingdom bac --threads 8

  # All HQ/MQ/LQ tiers from organize-mags
  mamisa run-mimag-rna \\
    --selected-dir filtered/Selected/ --output mimag_rna/ \\
    --extension fa

  # Mixed bacteria + archaea: pick bac/arc per genome from GTDB-Tk
  mamisa run-mimag-rna \\
    --selected-dir filtered/Selected/ --output mimag_rna/ \\
    --extension fa --tiers HQ \\
    --split-by-domain --gtdbtk-dir taxonomy/
        """
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--selected-dir', type=Path,
                             help='Directory with tier subdirectories (HQ/MQ/LQ)')
    input_group.add_argument('--genome-dir', type=Path,
                             help='Single directory with genome files')

    parser.add_argument('-o', '--output', type=Path, required=True,
                        help='Output directory for rRNA/tRNA results')
    parser.add_argument('--extension', default='fa',
                        help='Genome file extension (default: fa)')
    parser.add_argument('--kingdom', default='bac',
                        choices=['bac', 'arc', 'euk'],
                        help='Kingdom for barrnap/tRNAscan-SE (default: bac). '
                             'With --split-by-domain this is only the fallback '
                             'for genomes missing from GTDB-Tk')
    parser.add_argument('--split-by-domain', action='store_true',
                        help='Pick bac/arc per genome from GTDB-Tk results '
                             '(requires --gtdbtk-dir). Handles mixed '
                             'bacteria+archaea bin sets in one run')
    parser.add_argument('--gtdbtk-dir', type=Path,
                        help='GTDB-Tk output dir (summaries) for --split-by-domain')
    parser.add_argument('--threads', type=int, default=1,
                        help='Threads for barrnap (default: 1)')
    parser.add_argument('--tiers', default='HQ,MQ,LQ',
                        help='Comma-separated tiers to process (default: HQ,MQ,LQ)')
    parser.add_argument('--keep-intermediate', action='store_true',
                        help='Keep per-genome barrnap/tRNAscan output files')

    parser.set_defaults(func=run)
    return parser


def _process_dir(input_dir: Path, output_dir: Path, extension: str,
                 kingdom_for, threads: int) -> List[dict]:
    """
    Process every genome in a directory; return list of records.

    `kingdom_for` is a callable genome_path -> 'bac'/'arc'/'euk', so the
    kingdom can be fixed or resolved per genome (e.g. from GTDB-Tk).
    """
    genomes = collect_genomes(input_dir, extension)
    if not genomes:
        log_warning(f"No genomes in {input_dir} with extension .{extension}")
        return []

    work_dir = output_dir / "intermediate"
    work_dir.mkdir(parents=True, exist_ok=True)

    records = []
    for i, genome in enumerate(genomes, 1):
        kingdom = kingdom_for(genome)
        log_info(f"  [{i}/{len(genomes)}] {genome.name} [{kingdom}]")
        rec = process_genome(genome, work_dir, kingdom, threads)
        if rec:
            rec['kingdom'] = kingdom
            records.append(rec)
    return records


def run(args):
    print_header("MaMISA - MIMAG rRNA/tRNA Check")

    print_section("Checking Dependencies")
    if not check_mimag_env():
        sys.exit(1)

    # Build the per-genome kingdom resolver (fixed, or split by GTDB-Tk domain).
    if args.split_by_domain:
        if not args.gtdbtk_dir:
            log_error("--split-by-domain requires --gtdbtk-dir")
            sys.exit(1)
        validate_dir_exists(args.gtdbtk_dir, "GTDB-Tk directory")
        print_section("Loading domains from GTDB-Tk")
        domains = load_domains_from_gtdbtk(args.gtdbtk_dir)
        if not domains:
            log_error("No domains loaded from GTDB-Tk; cannot split by domain")
            sys.exit(1)
        kingdom_for = lambda g: resolve_kingdom(g, domains, args.kingdom)
    else:
        kingdom_for = lambda g: args.kingdom

    args.output.mkdir(parents=True, exist_ok=True)

    print_section("Running barrnap + tRNAscan-SE")
    all_records: List[dict] = []

    if args.selected_dir:
        validate_dir_exists(args.selected_dir, "Selected directory")
        tiers = [t.strip() for t in args.tiers.split(',')]
        for tier in tiers:
            tier_dir = args.selected_dir / tier
            if not tier_dir.exists():
                log_warning(f"Tier directory not found: {tier_dir}")
                continue
            log_info(f"\nTier: {tier}")
            tier_out = args.output / tier
            recs = _process_dir(tier_dir, tier_out, args.extension,
                                kingdom_for, args.threads)
            for r in recs:
                r['tier'] = tier
            all_records.extend(recs)
    else:
        validate_dir_exists(args.genome_dir, "Genome directory")
        all_records = _process_dir(args.genome_dir, args.output,
                                   args.extension, kingdom_for, args.threads)

    if not all_records:
        log_error("No genomes processed")
        sys.exit(1)

    summary_path = args.output / "mimag_rna_summary.tsv"
    write_mimag_rna_summary(all_records, summary_path)
    log_info(f"\nMIMAG rRNA/tRNA summary: {summary_path}")

    if not args.keep_intermediate:
        import shutil
        for sub in args.output.rglob("intermediate"):
            shutil.rmtree(sub, ignore_errors=True)

    print_section("SUMMARY")
    n_ok = sum(1 for r in all_records if r['mimag_rna_ok'])
    n_rrna = sum(1 for r in all_records if r['rrna_complete'])
    n_trna = sum(1 for r in all_records if r['trna_complete'])
    total = len(all_records)
    print(f"  Genomes checked:              {total:>8,}")
    print(f"  rRNA complete (5S+16S+23S):   {n_rrna:>8,}")
    print(f"  tRNA complete (>=18 aa):      {n_trna:>8,}")
    print(f"  Full MIMAG rRNA/tRNA pass:    {n_ok:>8,}")
    if args.split_by_domain:
        n_bac = sum(1 for r in all_records if r.get('kingdom') == 'bac')
        n_arc = sum(1 for r in all_records if r.get('kingdom') == 'arc')
        print(f"  By domain:                    bac={n_bac:,}  arc={n_arc:,}")

    log_info("\nDone. Feed to organize-mags with --mimag-rna-dir")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    register_parser(parser.add_subparsers(dest='command'))
    args = parser.parse_args()
    run(args)
