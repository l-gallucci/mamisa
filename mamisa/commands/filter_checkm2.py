#!/usr/bin/env python3
"""
MaMISA - filter-checkm2 command
Filter genomes based on CheckM2 quality reports
"""

import sys
import argparse
import csv
import shutil
from pathlib import Path
from typing import Dict, List

from ..utils.checkm2 import parse_all_reports
from ..utils.validation import validate_dir_exists, validate_file_exists
from ..utils.logging import log_info, log_error, log_warning, print_header, print_section


def load_name_map(map_file: Path) -> Dict[str, str]:
    """Load genome name mapping from a two-column TSV file."""
    name_map = {}

    with open(map_file) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if row and not row[0].startswith('#') and len(row) >= 2:
                name_map[row[0].strip()] = row[1].strip()

    log_info(f"Loaded {len(name_map):,} name mappings")
    return name_map


def normalize_name(name: str, strip_prefix: str = "", strip_suffix: str = "",
                   add_prefix: str = "", add_suffix: str = "") -> str:
    """Apply prefix/suffix transformations to a genome name."""
    result = name
    if strip_prefix and result.startswith(strip_prefix):
        result = result[len(strip_prefix):]
    if strip_suffix and result.endswith(strip_suffix):
        result = result[:-len(strip_suffix)]
    return f"{add_prefix}{result}{add_suffix}"


def find_genome_file(genomes_dir: Path, basename: str, extensions: List[str]) -> Path:
    """Find a genome file matching basename with any of the given extensions."""
    # Try direct match (basename already has an extension)
    for ext in extensions:
        if basename.endswith(f".{ext}"):
            candidates = list(genomes_dir.rglob(basename))
            if candidates:
                return candidates[0].resolve()

    # Try appending each extension
    for ext in extensions:
        candidates = list(genomes_dir.rglob(f"{basename}.{ext}"))
        if candidates:
            if len(candidates) > 1:
                log_warning(f"Multiple matches for {basename}.{ext}, using first: {candidates[0]}")
            return candidates[0].resolve()

    return None


def link_or_copy_genomes(records: List[dict], genomes_dir: Path, output_dir: Path,
                         selected_tiers: List[str], extensions: List[str],
                         name_map: Dict[str, str], normalize_params: dict,
                         mode: str = "symlink", dry_run: bool = False) -> Dict:
    """Link or copy genome files into tier-specific subdirectories."""

    stats = {
        'total_genomes': len(records),
        'selected_tiers': {tier: 0 for tier in selected_tiers},
        'found': 0, 'not_found': 0, 'linked_or_copied': 0,
    }

    if not dry_run:
        for tier in selected_tiers:
            (output_dir / "Selected" / tier).mkdir(parents=True, exist_ok=True)

    for record in records:
        tier = record['tier']
        if tier not in selected_tiers:
            continue

        stats['selected_tiers'][tier] += 1
        name = record['name']

        basename = name_map.get(name) or normalize_name(
            name,
            normalize_params.get('strip_prefix', ''),
            normalize_params.get('strip_suffix', ''),
            normalize_params.get('add_prefix', ''),
            normalize_params.get('add_suffix', ''),
        )

        genome_file = find_genome_file(genomes_dir, basename, extensions)

        if not genome_file:
            log_warning(f"Genome not found: {basename}")
            stats['not_found'] += 1
            continue

        stats['found'] += 1
        dest = output_dir / "Selected" / tier / genome_file.name

        if dry_run:
            log_info(f"Would {mode}: {genome_file.name} -> {tier}/")
            stats['linked_or_copied'] += 1
            continue

        try:
            if mode == "symlink":
                dest.unlink(missing_ok=True)
                dest.symlink_to(genome_file)
                log_info(f"Symlinked: {genome_file.name} -> {tier}/")
            else:
                shutil.copy2(genome_file, dest)
                log_info(f"Copied: {genome_file.name} -> {tier}/")
            stats['linked_or_copied'] += 1
        except Exception as e:
            log_error(f"Failed to {mode} {genome_file.name}: {e}")

    return stats


def register_parser(subparsers):
    """Register this command's parser"""
    parser = subparsers.add_parser(
        'filter-checkm2',
        help='Filter genomes based on CheckM2 quality reports',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Filter with default MiMAG thresholds
  mamisa filter-checkm2 \\
    --checkm2-root checkm2_results/ \\
    --genomes-dir genomes/ \\
    --output filtered/

  # Select only HQ and MQ genomes with custom thresholds
  mamisa filter-checkm2 \\
    --checkm2-root checkm2_results/ \\
    --genomes-dir genomes/ \\
    --output filtered/ \\
    --tiers HQ,MQ \\
    --hq-comp-min 95 \\
    --hq-cont-max 3

  # Copy instead of symlink
  mamisa filter-checkm2 \\
    --checkm2-root checkm2_results/ \\
    --genomes-dir genomes/ \\
    --output filtered/ \\
    --copy
        """
    )

    parser.add_argument('--checkm2-root', type=Path, required=True,
                        help='Root directory with CheckM2 results')
    parser.add_argument('--genomes-dir', type=Path, required=True,
                        help='Directory containing genome files')
    parser.add_argument('-o', '--output', type=Path, required=True,
                        help='Output directory')

    # Quality thresholds
    parser.add_argument('--hq-comp-min', type=float, default=90.0,
                        help='HQ minimum completeness (default: 90)')
    parser.add_argument('--hq-cont-max', type=float, default=5.0,
                        help='HQ maximum contamination (default: 5)')
    parser.add_argument('--mq-comp-min', type=float, default=70.0,
                        help='MQ minimum completeness (default: 70)')
    parser.add_argument('--mq-cont-max', type=float, default=10.0,
                        help='MQ maximum contamination (default: 10)')
    parser.add_argument('--lq-comp-min', type=float, default=50.0,
                        help='LQ minimum completeness (default: 50)')
    parser.add_argument('--lq-cont-max', type=float, default=10.0,
                        help='LQ maximum contamination (default: 10)')

    parser.add_argument('--tiers', default='HQ,MQ,LQ',
                        help='Comma-separated tiers to select (default: HQ,MQ,LQ)')
    parser.add_argument('--extensions', default='fa,fasta,fna,fa.gz,fasta.gz,fna.gz',
                        help='Comma-separated genome file extensions')

    # Symlink vs copy — both map to args.copy; default is symlink (copy=False)
    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument('--symlink', dest='copy', action='store_false',
                            help='Create symlinks (default)')
    mode_group.add_argument('--copy', dest='copy', action='store_true',
                            help='Copy files instead of symlinking')
    parser.set_defaults(copy=False)

    # Name normalization
    parser.add_argument('--name-prefix', default='', help='Add prefix to genome names')
    parser.add_argument('--name-suffix', default='', help='Add suffix to genome names')
    parser.add_argument('--strip-prefix', default='', help='Remove prefix from genome names')
    parser.add_argument('--strip-suffix', default='', help='Remove suffix from genome names')
    parser.add_argument('--name-map', type=Path,
                        help='TSV file mapping original names to search names')

    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be done without doing it')

    parser.set_defaults(func=run)
    return parser


def run(args):
    """Execute the filter-checkm2 command"""

    validate_dir_exists(args.checkm2_root, "CheckM2 root directory")
    validate_dir_exists(args.genomes_dir, "Genomes directory")

    if not args.dry_run:
        args.output.mkdir(parents=True, exist_ok=True)

    name_map = {}
    if args.name_map:
        validate_file_exists(args.name_map, "Name map file")
        name_map = load_name_map(args.name_map)

    print_header("MaMISA - Filter CheckM2 Results")

    thresholds = {
        'hq_comp': args.hq_comp_min, 'hq_cont': args.hq_cont_max,
        'mq_comp': args.mq_comp_min, 'mq_cont': args.mq_cont_max,
        'lq_comp': args.lq_comp_min, 'lq_cont': args.lq_cont_max,
    }

    print_section("Parsing CheckM2 Reports")
    records, tier_counts = parse_all_reports(args.checkm2_root, thresholds)

    print_section("Quality Tier Summary")
    for tier in ['HQ', 'MQ', 'LQ', 'Fail']:
        print(f"  {tier:4s}: {tier_counts.get(tier, 0):>8,} genomes")
    print(f"  {'Total':4s}: {len(records):>8,} genomes")

    merged_file = args.output / "merged_quality.tsv"
    if not args.dry_run:
        with open(merged_file, 'w', newline='') as f:
            writer = csv.DictWriter(
                f,
                fieldnames=['report_path', 'name', 'completeness', 'contamination', 'tier'],
                delimiter='\t',
            )
            writer.writeheader()
            writer.writerows(records)
        log_info(f"\nMerged quality report: {merged_file}")

    print_section("Organizing Genomes by Tier")

    selected_tiers = [t.strip() for t in args.tiers.split(',')]
    extensions = [e.strip() for e in args.extensions.split(',')]
    mode = "copy" if args.copy else "symlink"

    normalize_params = {
        'strip_prefix': args.strip_prefix,
        'strip_suffix': args.strip_suffix,
        'add_prefix': args.name_prefix,
        'add_suffix': args.name_suffix,
    }

    stats = link_or_copy_genomes(
        records=records,
        genomes_dir=args.genomes_dir,
        output_dir=args.output,
        selected_tiers=selected_tiers,
        extensions=extensions,
        name_map=name_map,
        normalize_params=normalize_params,
        mode=mode,
        dry_run=args.dry_run,
    )

    print_section("RESULTS")
    print(f"  Total genomes in reports:    {stats['total_genomes']:>8,}")
    for tier in selected_tiers:
        print(f"  {tier} genomes selected:        {stats['selected_tiers'].get(tier, 0):>8,}")
    print(f"  Genome files found:          {stats['found']:>8,}")
    print(f"  Genome files not found:      {stats['not_found']:>8,}")
    print(f"  Files {mode}ed:             {stats['linked_or_copied']:>8,}")

    if args.dry_run:
        log_warning("\nDRY-RUN mode: No files were actually created")

    log_info("\n✓ Done!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    register_parser(parser.add_subparsers(dest='command'))
    args = parser.parse_args()
    run(args)
