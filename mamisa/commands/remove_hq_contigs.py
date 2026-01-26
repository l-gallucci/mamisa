#!/usr/bin/env python3
"""
MaMISA - remove-hq-contigs command
Remove high-quality genome contigs from assembly
"""

import sys
import argparse
from pathlib import Path
from typing import Set
import re

from ..utils.fasta import read_fasta_streaming, write_fasta
from ..utils.validation import validate_file_exists, validate_dir_exists
from ..utils.logging import log_info, log_error, log_warning, print_header, print_section


def extract_contig_id(text: str) -> str:
    """
    Extract contig identifier in multiple formats:
    - contig_123 → 123
    - stdin.part_contig_7420 → 7420
    - NODE_456 → 456
    Returns the first number found or full text if no number
    """
    match = re.search(r'(\d+)', text)
    if match:
        return match.group(1)
    return text.strip()


def load_hq_ids_from_dir(hq_dir: Path) -> Set[str]:
    """
    Extract contig IDs from HQ genome directory
    Scans all files and extracts numeric IDs from filenames
    """
    hq_ids = set()
    
    log_info(f"Scanning HQ directory: {hq_dir}")
    
    count = 0
    for filepath in hq_dir.rglob('*'):
        if filepath.is_file() and not filepath.name.startswith('.'):
            file_id = extract_contig_id(filepath.stem)
            if file_id:
                hq_ids.add(file_id)
                count += 1
    
    log_info(f"Found {len(hq_ids):,} unique HQ genome IDs from {count:,} files")
    return hq_ids


def load_hq_ids_from_list(list_file: Path) -> Set[str]:
    """
    Load HQ contig IDs from a text file (one per line)
    """
    hq_ids = set()
    
    log_info(f"Loading HQ IDs from list: {list_file}")
    
    with open(list_file) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                contig_id = extract_contig_id(line)
                if contig_id:
                    hq_ids.add(contig_id)
    
    log_info(f"Loaded {len(hq_ids):,} HQ genome IDs")
    return hq_ids


def process_assembly(assembly: Path, hq_ids: Set[str], output: Path,
                     min_length: int = 0, dry_run: bool = False) -> dict:
    """
    Remove HQ contigs from assembly
    """
    
    print_header("MaMISA - Remove HQ Contigs")
    
    print_section("Configuration")
    log_info(f"Assembly: {assembly}")
    log_info(f"HQ contigs to remove: {len(hq_ids):,}")
    log_info(f"Minimum length filter: {min_length:,} bp")
    if not dry_run:
        log_info(f"Output: {output}")
    
    # Statistics
    stats = {
        'total_contigs': 0,
        'hq_removed': 0,
        'too_short_removed': 0,
        'kept': 0,
        'total_bp': 0,
        'removed_bp': 0,
        'kept_bp': 0
    }
    
    output_sequences = []
    
    print_section("Processing assembly")
    
    for contig_name, sequence in read_fasta_streaming(assembly):
        stats['total_contigs'] += 1
        contig_len = len(sequence)
        stats['total_bp'] += contig_len
        
        # Extract ID and check if it's HQ
        contig_id = extract_contig_id(contig_name)
        
        if contig_id in hq_ids:
            # This is an HQ contig - remove it
            stats['hq_removed'] += 1
            stats['removed_bp'] += contig_len
            continue
        
        if contig_len < min_length:
            # Too short - remove it
            stats['too_short_removed'] += 1
            stats['removed_bp'] += contig_len
            continue
        
        # Keep this contig
        stats['kept'] += 1
        stats['kept_bp'] += contig_len
        
        if not dry_run:
            output_sequences.append((contig_name, sequence))
    
    # Write output
    if not dry_run and output_sequences:
        log_info(f"\nWriting output to: {output}")
        write_fasta(output_sequences, output)
    
    # Print statistics
    print_section("RESULTS")
    print(f"  Total contigs in input:     {stats['total_contigs']:>10,}")
    print(f"  HQ contigs removed:         {stats['hq_removed']:>10,}")
    print(f"  Too short (removed):        {stats['too_short_removed']:>10,}")
    print(f"  Contigs kept:               {stats['kept']:>10,}")
    print()
    print(f"  Total bases in input:       {stats['total_bp']:>10,} bp")
    print(f"  Bases removed:              {stats['removed_bp']:>10,} bp")
    print(f"  Bases kept:                 {stats['kept_bp']:>10,} bp")
    
    if stats['total_bp'] > 0:
        pct_kept = 100 * stats['kept_bp'] / stats['total_bp']
        pct_removed = 100 - pct_kept
        print(f"  Retention rate:             {pct_kept:>9.2f}%")
        print(f"  Removal rate:               {pct_removed:>9.2f}%")
    
    if dry_run:
        log_warning("\nDRY-RUN mode: No output file generated")
    
    return stats


def register_parser(subparsers):
    """Register this command's parser"""
    parser = subparsers.add_parser(
        'remove-hq-contigs',
        help='Remove high-quality genome contigs from assembly',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Remove HQ contigs using directory
  mamisa remove-hq-contigs \\
    --assembly assembly.fa \\
    --hq-dir hq_genomes/ \\
    --output filtered.fa
  
  # Remove HQ contigs using list file
  mamisa remove-hq-contigs \\
    --assembly assembly.fa \\
    --hq-list hq_ids.txt \\
    --output filtered.fa \\
    --min-length 1000
  
  # Dry-run to see statistics
  mamisa remove-hq-contigs \\
    --assembly assembly.fa \\
    --hq-dir hq_genomes/ \\
    --dry-run
        """
    )
    
    parser.add_argument('-a', '--assembly', type=Path, required=True,
                        help='Input assembly FASTA file')
    
    hq_group = parser.add_mutually_exclusive_group(required=True)
    hq_group.add_argument('--hq-dir', type=Path,
                          help='Directory containing HQ genome files')
    hq_group.add_argument('--hq-list', type=Path,
                          help='Text file with list of HQ contig IDs')
    
    parser.add_argument('-o', '--output', type=Path,
                        help='Output filtered assembly file')
    parser.add_argument('-l', '--min-length', type=int, default=0,
                        help='Minimum contig length to keep (default: 0)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show statistics without generating output')
    parser.add_argument('--stats', type=Path,
                        help='Save statistics to TSV file')
    
    parser.set_defaults(func=run)
    return parser


def run(args):
    """Execute the remove-hq-contigs command"""
    import argparse
    
    # Validation
    if not args.dry_run and not args.output:
        raise ValueError("--output is required unless using --dry-run")
    
    validate_file_exists(args.assembly, "Assembly file")
    
    # Load HQ IDs
    if args.hq_dir:
        validate_dir_exists(args.hq_dir, "HQ directory")
        hq_ids = load_hq_ids_from_dir(args.hq_dir)
    else:
        validate_file_exists(args.hq_list, "HQ list file")
        hq_ids = load_hq_ids_from_list(args.hq_list)
    
    if not hq_ids:
        log_error("No HQ contig IDs found!")
        sys.exit(1)
    
    # Process
    stats = process_assembly(
        assembly=args.assembly,
        hq_ids=hq_ids,
        output=args.output,
        min_length=args.min_length,
        dry_run=args.dry_run
    )
    
    # Save stats
    if args.stats and not args.dry_run:
        with open(args.stats, 'w') as f:
            f.write("metric\tvalue\n")
            for key, value in stats.items():
                f.write(f"{key}\t{value}\n")
        log_info(f"\nStatistics saved to: {args.stats}")
    
    log_info("\n✓ Done!")


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    register_parser(parser.add_subparsers(dest='command'))
    args = parser.parse_args()
    run(args)
