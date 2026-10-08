#!/usr/bin/env python3
"""
MaMISA - organize-mags command (formerly filter-checkm2)

Parse CheckM2 quality reports, assign MiMAG quality tiers, split genomes into
HQ/MQ/LQ subdirectories, and optionally rename each genome with its sample of
origin and GTDB-Tk taxonomy.

The command is registered as `organize-mags`; `filter-checkm2` is kept as a
deprecated alias for backward compatibility.
"""

import sys
import re
import argparse
import csv
import shutil
from pathlib import Path
from typing import Dict, List, Optional

from ..utils.checkm2 import parse_all_reports, NoReportsFoundError
from ..utils.chimera import (
    parse_gtdbtk_summary,
    extract_taxonomy_level,
)
from ..utils.validation import validate_dir_exists, validate_file_exists
from ..utils.logging import log_info, log_error, log_warning, print_header, print_section


# Default: sample is the leading token before the first '.' or '_'
#   SAMPLE1.bin.3  -> SAMPLE1
#   SAMPLE1_bin_3  -> SAMPLE1
DEFAULT_SAMPLE_REGEX = r'^([^._]+)'

# Order used when falling back to a higher rank for the taxonomy label
_RANK_FALLBACK = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'domain']


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


def find_genome_file(genomes_dir: Path, basename: str, extensions: List[str]) -> Optional[Path]:
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


# ---------------------------------------------------------------------------
# Naming helpers
# ---------------------------------------------------------------------------

def sanitize_token(token: str) -> str:
    """Make a token safe for use inside a filename."""
    token = token.strip()
    # collapse whitespace and path/field separators to underscores
    token = re.sub(r'[\s;/\\]+', '_', token)
    # drop anything else that is not filename-friendly
    token = re.sub(r'[^A-Za-z0-9_.+-]', '', token)
    return token or 'NA'


def extract_sample(name: str, sample_regex: str) -> str:
    """Extract the sample of origin from a genome name via a capturing regex."""
    try:
        m = re.search(sample_regex, name)
    except re.error as e:
        raise ValueError(f"Invalid --sample-regex {sample_regex!r}: {e}")
    if m and m.groups():
        return sanitize_token(m.group(1))
    if m and m.group(0):
        return sanitize_token(m.group(0))
    return 'NA'


def taxonomy_label(classification: str, level: str) -> str:
    """
    Build a taxonomy token for the chosen rank, falling back to higher ranks
    when the requested rank is unclassified/empty.
    """
    if not classification:
        return 'Unclassified'

    taxon = extract_taxonomy_level(classification, level)
    if taxon and taxon != 'Unclassified':
        return sanitize_token(taxon)

    # fall back up the ranks starting just above the requested one
    try:
        start = _RANK_FALLBACK.index(level) + 1
    except ValueError:
        start = 0
    for higher in _RANK_FALLBACK[start:]:
        taxon = extract_taxonomy_level(classification, higher)
        if taxon and taxon != 'Unclassified':
            return sanitize_token(f"{taxon}")
    return 'Unclassified'


def load_taxonomy(gtdbtk_dir: Path) -> Dict[str, Dict]:
    """
    Load and merge all GTDB-Tk summary TSVs (bac120 + ar53) found under a
    directory. Returns {genome_name: taxonomy_record}.
    """
    summaries = sorted(gtdbtk_dir.rglob("gtdbtk.*.summary.tsv"))
    if not summaries:
        # also accept generic *summary.tsv
        summaries = sorted(gtdbtk_dir.rglob("*summary.tsv"))

    if not summaries:
        log_warning(f"No GTDB-Tk summary TSV found under {gtdbtk_dir}")
        return {}

    merged: Dict[str, Dict] = {}
    for summary in summaries:
        merged.update(parse_gtdbtk_summary(summary))
    log_info(f"Loaded taxonomy for {len(merged):,} genomes from {len(summaries)} summary file(s)")
    return merged


def lookup_taxonomy(tax: Dict[str, Dict], record_name: str,
                    genome_stem: str) -> Optional[str]:
    """Resolve a classification string for a genome by trying several keys."""
    for key in (record_name, genome_stem, record_name.rsplit('.', 1)[0]):
        if key in tax:
            return tax[key].get('classification', '')
    return None


def build_output_name(genome_file: Path, record_name: str,
                      rename: bool, tax: Dict[str, Dict],
                      sample_regex: str, tax_level: str) -> str:
    """Compute the destination filename for a genome."""
    if not rename:
        return genome_file.name

    stem = genome_file.name
    # strip known compound extensions for the stem used in the new name
    for ext in ('.fa.gz', '.fasta.gz', '.fna.gz', '.fa', '.fasta', '.fna'):
        if stem.endswith(ext):
            suffix = ext
            stem = stem[:-len(ext)]
            break
    else:
        suffix = genome_file.suffix
        stem = genome_file.stem

    sample = extract_sample(record_name, sample_regex)

    taxon = 'NoTax'
    if tax:
        classification = lookup_taxonomy(tax, record_name, stem)
        if classification is not None:
            taxon = taxonomy_label(classification, tax_level)

    return f"{sanitize_token(sample)}__{taxon}__{sanitize_token(stem)}{suffix}"


def link_or_copy_genomes(records: List[dict], genomes_dir: Path, output_dir: Path,
                         selected_tiers: List[str], extensions: List[str],
                         name_map: Dict[str, str], normalize_params: dict,
                         mode: str = "symlink", dry_run: bool = False,
                         rename: bool = False, tax: Optional[Dict] = None,
                         sample_regex: str = DEFAULT_SAMPLE_REGEX,
                         tax_level: str = 'genus') -> Dict:
    """Link or copy genome files into tier-specific subdirectories."""

    tax = tax or {}
    stats = {
        'total_genomes': len(records),
        'selected_tiers': {tier: 0 for tier in selected_tiers},
        'found': 0, 'not_found': 0, 'linked_or_copied': 0, 'renamed': 0,
    }

    if not dry_run:
        for tier in selected_tiers:
            (output_dir / "Selected" / tier).mkdir(parents=True, exist_ok=True)

    # map written destination -> source, to detect collisions after renaming
    seen_dest: Dict[Path, Path] = {}

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

        out_name = build_output_name(
            genome_file, name, rename, tax, sample_regex, tax_level)
        if rename and out_name != genome_file.name:
            stats['renamed'] += 1

        dest = output_dir / "Selected" / tier / out_name

        if dest in seen_dest:
            log_warning(f"Name collision for {dest.name} "
                        f"({seen_dest[dest].name} vs {genome_file.name}); skipping second")
            continue
        seen_dest[dest] = genome_file

        if dry_run:
            log_info(f"Would {mode}: {genome_file.name} -> {tier}/{out_name}")
            stats['linked_or_copied'] += 1
            continue

        try:
            if mode == "symlink":
                dest.unlink(missing_ok=True)
                dest.symlink_to(genome_file)
                log_info(f"Symlinked: {genome_file.name} -> {tier}/{out_name}")
            else:
                shutil.copy2(genome_file, dest)
                log_info(f"Copied: {genome_file.name} -> {tier}/{out_name}")
            stats['linked_or_copied'] += 1
        except Exception as e:
            log_error(f"Failed to {mode} {genome_file.name}: {e}")

    return stats


def register_parser(subparsers):
    """Register this command's parser"""
    parser = subparsers.add_parser(
        'organize-mags',
        aliases=['filter-checkm2'],
        help='Split genomes into HQ/MQ/LQ tiers, optionally rename by sample + taxonomy',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Split with default MiMAG thresholds (symlinks)
  mamisa organize-mags \\
    --checkm2-root checkm2_results/ \\
    --genomes-dir genomes/ \\
    --output filtered/

  # Rename each genome <sample>__<taxon>__<orig> using GTDB-Tk taxonomy
  mamisa organize-mags \\
    --checkm2-root checkm2_results/ \\
    --genomes-dir genomes/ \\
    --gtdbtk-dir gtdbtk_results/ \\
    --output filtered/ \\
    --rename --copy --tax-level genus

  # Custom sample parsing (capture group 1 = sample id)
  mamisa organize-mags \\
    --checkm2-root checkm2_results/ --genomes-dir genomes/ \\
    --output filtered/ --rename \\
    --sample-regex '^(S[0-9]+)_'
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

    # Taxonomy-aware renaming
    rename_group = parser.add_argument_group('renaming (sample + taxonomy)')
    rename_group.add_argument('--rename', action='store_true',
                              help='Rename outputs as <sample>__<taxon>__<orig>')
    rename_group.add_argument('--gtdbtk-dir', type=Path,
                              help='Directory with GTDB-Tk summary TSVs (for --rename taxonomy)')
    rename_group.add_argument('--sample-regex', default=DEFAULT_SAMPLE_REGEX,
                              help=f"Regex with a capture group for the sample id "
                                   f"(default: {DEFAULT_SAMPLE_REGEX!r} = token before first . or _)")
    rename_group.add_argument('--tax-level', default='genus',
                              choices=['domain', 'phylum', 'class', 'order',
                                       'family', 'genus', 'species'],
                              help='Taxonomic rank used in the name (default: genus)')

    # Name normalization (for locating files on disk)
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
    """Execute the organize-mags command"""

    # Deprecation notice when invoked via the old alias
    if getattr(args, 'command', '') == 'filter-checkm2':
        log_warning("'filter-checkm2' is deprecated; use 'organize-mags' instead.")

    validate_dir_exists(args.checkm2_root, "CheckM2 root directory")
    validate_dir_exists(args.genomes_dir, "Genomes directory")

    if not args.dry_run:
        args.output.mkdir(parents=True, exist_ok=True)

    name_map = {}
    if args.name_map:
        validate_file_exists(args.name_map, "Name map file")
        name_map = load_name_map(args.name_map)

    print_header("MaMISA - Organize MAGs by Quality Tier")

    # Load taxonomy if requested
    tax = {}
    if args.rename and args.gtdbtk_dir:
        validate_dir_exists(args.gtdbtk_dir, "GTDB-Tk directory")
        print_section("Loading GTDB-Tk Taxonomy")
        tax = load_taxonomy(args.gtdbtk_dir)
    elif args.rename and not args.gtdbtk_dir:
        log_warning("--rename without --gtdbtk-dir: taxonomy token will be 'NoTax'")

    thresholds = {
        'hq_comp': args.hq_comp_min, 'hq_cont': args.hq_cont_max,
        'mq_comp': args.mq_comp_min, 'mq_cont': args.mq_cont_max,
        'lq_comp': args.lq_comp_min, 'lq_cont': args.lq_cont_max,
    }

    print_section("Parsing CheckM2 Reports")
    try:
        records, tier_counts = parse_all_reports(args.checkm2_root, thresholds)
    except NoReportsFoundError as e:
        log_error(str(e))
        sys.exit(1)

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

    if args.rename and mode == "symlink":
        log_warning("Renaming with symlinks: links point to original files under "
                    "their new names (use --copy for standalone renamed files).")

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
        rename=args.rename,
        tax=tax,
        sample_regex=args.sample_regex,
        tax_level=args.tax_level,
    )

    print_section("RESULTS")
    print(f"  Total genomes in reports:    {stats['total_genomes']:>8,}")
    for tier in selected_tiers:
        print(f"  {tier} genomes selected:        {stats['selected_tiers'].get(tier, 0):>8,}")
    print(f"  Genome files found:          {stats['found']:>8,}")
    print(f"  Genome files not found:      {stats['not_found']:>8,}")
    print(f"  Files {mode}ed:             {stats['linked_or_copied']:>8,}")
    if args.rename:
        print(f"  Files renamed:               {stats['renamed']:>8,}")

    if args.dry_run:
        log_warning("\nDRY-RUN mode: No files were actually created")

    log_info("\n✓ Done!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    register_parser(parser.add_subparsers(dest='command'))
    args = parser.parse_args()
    run(args)
