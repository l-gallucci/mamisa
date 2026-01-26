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
from typing import Dict, List, Tuple
from collections import defaultdict

from ..utils.validation import validate_dir_exists, validate_file_exists
from ..utils.logging import log_info, log_error, log_warning, print_header, print_section


def load_name_map(map_file: Path) -> Dict[str, str]:
    """Load genome name mapping from TSV file"""
    name_map = {}
    
    with open(map_file) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if row and not row[0].startswith('#') and len(row) >= 2:
                name_map[row[0].strip()] = row[1].strip()
    
    log_info(f"Loaded {len(name_map):,} name mappings")
    return name_map


def clean_percent(value: str) -> float:
    """Convert percentage string to float"""
    if value is None:
        return float('nan')
    value = value.strip().replace('%', '')
    try:
        return float(value)
    except:
        return float('nan')


def normalize_name(name: str, strip_prefix: str = "", strip_suffix: str = "",
                   add_prefix: str = "", add_suffix: str = "") -> str:
    """Normalize genome name with prefix/suffix operations"""
    result = name
    
    if strip_prefix and result.startswith(strip_prefix):
        result = result[len(strip_prefix):]
    
    if strip_suffix and result.endswith(strip_suffix):
        result = result[:-len(strip_suffix)]
    
    return f"{add_prefix}{result}{add_suffix}"


def assign_tier(completeness: float, contamination: float,
                hq_comp: float, hq_cont: float,
                mq_comp: float, mq_cont: float,
                lq_comp: float, lq_cont: float) -> str:
    """
    Assign quality tier based on MiMAG standards
    
    Returns: HQ, MQ, LQ, or Fail
    """
    if completeness >= hq_comp and contamination <= hq_cont:
        return "HQ"
    if completeness >= mq_comp and contamination <= mq_cont:
        return "MQ"
    if completeness >= lq_comp and contamination <= lq_cont:
        return "LQ"
    return "Fail"


def find_quality_reports(checkm2_root: Path) -> List[Path]:
    """Recursively find all quality_report.tsv files"""
    reports = sorted(checkm2_root.rglob("quality_report.tsv"))
    return reports


def parse_checkm2_reports(checkm2_root: Path, thresholds: dict) -> Tuple[List, Dict]:
    """
    Parse all CheckM2 quality reports
    
    Returns:
        Tuple of (all_records, tier_counts)
    """
    reports = find_quality_reports(checkm2_root)
    
    if not reports:
        log_error(f"No quality_report.tsv files found in {checkm2_root}")
        sys.exit(1)
    
    log_info(f"Found {len(reports):,} CheckM2 report(s)")
    
    all_records = []
    tier_counts = defaultdict(int)
    
    for report_path in reports:
        log_info(f"  Processing: {report_path}")
        
        with open(report_path, newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            
            required_cols = {'Name', 'Completeness', 'Contamination'}
            if not required_cols.issubset(reader.fieldnames or set()):
                log_warning(f"Missing required columns in {report_path}, skipping")
                continue
            
            for row in reader:
                name = (row.get('Name') or '').strip()
                if not name:
                    continue
                
                comp = clean_percent(row.get('Completeness'))
                cont = clean_percent(row.get('Contamination'))
                
                tier = assign_tier(
                    comp, cont,
                    thresholds['hq_comp'], thresholds['hq_cont'],
                    thresholds['mq_comp'], thresholds['mq_cont'],
                    thresholds['lq_comp'], thresholds['lq_cont']
                )
                
                all_records.append({
                    'report_path': str(report_path),
                    'name': name,
                    'completeness': comp,
                    'contamination': cont,
                    'tier': tier
                })
                
                tier_counts[tier] += 1
    
    return all_records, dict(tier_counts)


def find_genome_file(genomes_dir: Path, basename: str, extensions: List[str]) -> Path:
    """Find genome file matching basename and extensions"""
    
    # Check if basename already has an extension
    for ext in extensions:
        if basename.endswith(f".{ext}"):
            # Try direct match
            candidates = list(genomes_dir.rglob(basename))
            if candidates:
                return candidates[0].resolve()
    
    # Try adding extensions
    for ext in extensions:
        pattern = f"{basename}.{ext}"
        candidates = list(genomes_dir.rglob(pattern))
        if candidates:
            if len(candidates) > 1:
                log_warning(f"Multiple matches for {pattern}, using first: {candidates[0]}")
            return candidates[0].resolve()
    
    return None


def link_or_copy_genomes(records: List[dict], genomes_dir: Path, output_dir: Path,
                        selected_tiers: List[str], extensions: List[str],
                        name_map: Dict[str, str], normalize_params: dict,
                        mode: str = "symlink", dry_run: bool = False) -> Dict:
    """
    Link or copy genome files to tier-specific directories
    """
    
    stats = {
        'total_genomes': len(records),
        'selected_tiers': {},
        'found': 0,
        'not_found': 0,
        'linked_or_copied': 0
    }
    
    for tier in selected_tiers:
        stats['selected_tiers'][tier] = 0
    
    # Create output directories
    if not dry_run:
        for tier in selected_tiers:
            tier_dir = output_dir / "Selected" / tier
            tier_dir.mkdir(parents=True, exist_ok=True)
    
    # Process each genome
    for record in records:
        tier = record['tier']
        
        # Skip if not in selected tiers
        if tier not in selected_tiers:
            continue
        
        stats['selected_tiers'][tier] += 1
        
        name = record['name']
        
        # Normalize name
        if name in name_map:
            basename = name_map[name]
        else:
            basename = normalize_name(
                name,
                normalize_params.get('strip_prefix', ''),
                normalize_params.get('strip_suffix', ''),
                normalize_params.get('add_prefix', ''),
                normalize_params.get('add_suffix', '')
            )
        
        # Find genome file
        genome_file = find_genome_file(genomes_dir, basename, extensions)
        
        if not genome_file:
            log_warning(f"Genome not found: {basename}")
            stats['not_found'] += 1
            continue
        
        stats['found'] += 1
        
        # Destination path
        dest = output_dir / "Selected" / tier / genome_file.name
        
        if dry_run:
            log_info(f"Would {mode}: {genome_file.name} -> {tier}/")
            stats['linked_or_copied'] += 1
            continue
        
        # Link or copy
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
    
    # Input/Output
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
    
    # Tier selection
    parser.add_argument('--tiers', default='HQ,MQ,LQ',
                        help='Comma-separated tiers to select (default: HQ,MQ,LQ)')
    
    # File handling
    parser.add_argument('--extensions', default='fa,fasta,fna,fa.gz,fasta.gz,fna.gz',
                        help='Comma-separated genome file extensions')
    
    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument('--symlink', action='store_true', default=True,
                            help='Create symlinks (default)')
    mode_group.add_argument('--copy', action='store_true',
                            help='Copy files instead of symlinking')
    
    # Name normalization
    parser.add_argument('--name-prefix', default='',
                        help='Add prefix to genome names')
    parser.add_argument('--name-suffix', default='',
                        help='Add suffix to genome names')
    parser.add_argument('--strip-prefix', default='',
                        help='Remove prefix from genome names')
    parser.add_argument('--strip-suffix', default='',
                        help='Remove suffix from genome names')
    parser.add_argument('--name-map', type=Path,
                        help='TSV file mapping original names to search names')
    
    # Other options
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be done without doing it')
    
    parser.set_defaults(func=run)
    return parser


def run(args):
    """Execute the filter-checkm2 command"""
    import argparse
    
    # Validation
    validate_dir_exists(args.checkm2_root, "CheckM2 root directory")
    validate_dir_exists(args.genomes_dir, "Genomes directory")
    
    if not args.dry_run:
        args.output.mkdir(parents=True, exist_ok=True)
    
    # Load name map if provided
    name_map = {}
    if args.name_map:
        validate_file_exists(args.name_map, "Name map file")
        name_map = load_name_map(args.name_map)
    
    print_header("MaMISA - Filter CheckM2 Results")
    
    # Parse quality thresholds
    thresholds = {
        'hq_comp': args.hq_comp_min,
        'hq_cont': args.hq_cont_max,
        'mq_comp': args.mq_comp_min,
        'mq_cont': args.mq_cont_max,
        'lq_comp': args.lq_comp_min,
        'lq_cont': args.lq_cont_max
    }
    
    # Parse CheckM2 reports
    print_section("Parsing CheckM2 Reports")
    records, tier_counts = parse_checkm2_reports(args.checkm2_root, thresholds)
    
    # Print tier summary
    print_section("Quality Tier Summary")
    for tier in ['HQ', 'MQ', 'LQ', 'Fail']:
        count = tier_counts.get(tier, 0)
        print(f"  {tier:4s}: {count:>8,} genomes")
    print(f"  {'Total':4s}: {len(records):>8,} genomes")
    
    # Save merged quality report
    merged_file = args.output / "merged_quality.tsv"
    if not args.dry_run:
        with open(merged_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['report_path', 'name', 'completeness', 
                                                   'contamination', 'tier'], 
                                   delimiter='\t')
            writer.writeheader()
            writer.writerows(records)
        log_info(f"\nMerged quality report: {merged_file}")
    
    # Link or copy genomes
    print_section("Organizing Genomes by Tier")
    
    selected_tiers = [t.strip() for t in args.tiers.split(',')]
    extensions = [e.strip() for e in args.extensions.split(',')]
    mode = "copy" if args.copy else "symlink"
    
    normalize_params = {
        'strip_prefix': args.strip_prefix,
        'strip_suffix': args.strip_suffix,
        'add_prefix': args.name_prefix,
        'add_suffix': args.name_suffix
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
        dry_run=args.dry_run
    )
    
    # Print final statistics
    print_section("RESULTS")
    print(f"  Total genomes in reports:    {stats['total_genomes']:>8,}")
    for tier in selected_tiers:
        count = stats['selected_tiers'].get(tier, 0)
        print(f"  {tier} genomes selected:        {count:>8,}")
    print(f"  Genome files found:          {stats['found']:>8,}")
    print(f"  Genome files not found:      {stats['not_found']:>8,}")
    print(f"  Files {mode}ed:             {stats['linked_or_copied']:>8,}")
    
    if args.dry_run:
        log_warning("\nDRY-RUN mode: No files were actually created")
    
    log_info("\n✓ Done!")


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    register_parser(parser.add_subparsers(dest='command'))
    args = parser.parse_args()
    run(args)
