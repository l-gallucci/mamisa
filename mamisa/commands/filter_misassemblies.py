#!/usr/bin/env python3
"""
MaMISA - filter-misassemblies command
Integrated assembly filtering with HQ genome awareness and misassembly handling
"""

import sys
import argparse
from pathlib import Path
from collections import defaultdict
from typing import Set, Dict, List, Tuple
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
    - scaffold_789 → 789
    Returns the first number found or the original text if no number
    """
    match = re.search(r'(\d+)', text)
    if match:
        return match.group(1)
    return text.strip()


def parse_clipping_files(clipping_files: List[Path]) -> Dict[str, List[int]]:
    """
    Parse *-clipping.txt files from anvi'o
    Returns: dict[contig_name] = [list of clipping positions]
    """
    clipping_positions = defaultdict(list)
    
    for file_path in clipping_files:
        log_info(f"  Reading: {file_path.name}")
        with open(file_path) as f:
            header = f.readline()
            
            for line in f:
                if line.strip():
                    fields = line.strip().split('\t')
                    contig_name = fields[0]
                    position = int(fields[2])
                    clipping_positions[contig_name].append(position)
    
    for contig in clipping_positions:
        clipping_positions[contig].sort()
    
    return dict(clipping_positions)


def get_hq_contig_ids(hq_dir: Path) -> Set[str]:
    """
    Extract contig IDs from HQ genome filenames
    Handles: genome_123.fa, contig_456.fasta, etc.
    """
    hq_ids = set()
    
    if not hq_dir or not hq_dir.exists():
        log_warning(f"HQ directory not found: {hq_dir}")
        return hq_ids
    
    for filename in hq_dir.iterdir():
        if filename.name.startswith("."):
            continue
        
        file_id = extract_contig_id(filename.stem)
        if file_id:
            hq_ids.add(file_id)
    
    return hq_ids


def split_contig(contig_name: str, sequence: str, split_positions: List[int], 
                 min_length: int) -> List[Tuple[str, str]]:
    """
    Split a contig at given positions and return fragments >= min_length
    """
    fragments = []
    positions = [0] + split_positions + [len(sequence)]
    
    for i in range(len(positions) - 1):
        start = positions[i]
        end = positions[i + 1]
        fragment_seq = sequence[start:end]
        fragment_len = len(fragment_seq)
        
        if fragment_len >= min_length:
            if len(split_positions) > 0:
                fragment_name = f"{contig_name}_split{i+1}_len{fragment_len}"
            else:
                fragment_name = contig_name
            fragments.append((fragment_name, fragment_seq))
    
    return fragments


def classify_contigs(assembly: Path, hq_ids: Set[str], 
                     clipping_positions: Dict[str, List[int]]) -> Dict[str, Set[str]]:
    """
    Classify contigs into categories:
    - hq_clean: HQ genomes without misassemblies
    - hq_problematic: HQ genomes WITH misassemblies
    - nonhq_misassembly: Non-HQ with misassemblies
    - nonhq_clean: Non-HQ without misassemblies
    """
    categories = {
        'hq_clean': set(),
        'hq_problematic': set(),
        'nonhq_misassembly': set(),
        'nonhq_clean': set()
    }
    
    log_info("Classifying contigs...")
    
    for contig_name, sequence in read_fasta_streaming(assembly):
        contig_id = extract_contig_id(contig_name)
        is_hq = contig_id in hq_ids
        has_misassembly = contig_name in clipping_positions
        
        if is_hq and not has_misassembly:
            categories['hq_clean'].add(contig_name)
        elif is_hq and has_misassembly:
            categories['hq_problematic'].add(contig_name)
        elif not is_hq and has_misassembly:
            categories['nonhq_misassembly'].add(contig_name)
        else:
            categories['nonhq_clean'].add(contig_name)
    
    return categories


def process_assembly(assembly: Path, hq_dir: Path, misassembly_dir: Path,
                     output: Path, mode: str, min_length: int,
                     preserve_hq_with_issues: bool, dry_run: bool = False) -> dict:
    """
    Main processing function with integrated HQ handling
    """
    
    print_header("MaMISA - Integrated Assembly Filtering")
    
    # 1. Load HQ contig IDs
    print_section("Loading HQ genome IDs")
    log_info(f"HQ directory: {hq_dir}")
    hq_ids = get_hq_contig_ids(hq_dir) if hq_dir else set()
    log_info(f"Found {len(hq_ids):,} HQ genome IDs")
    
    # 2. Parse misassembly files
    print_section("Loading misassembly data")
    log_info(f"Misassembly directory: {misassembly_dir}")
    clipping_files = list(misassembly_dir.glob('*-clipping.txt'))
    if not clipping_files:
        log_error(f"No *-clipping.txt files found in {misassembly_dir}")
        sys.exit(1)
    
    clipping_positions = parse_clipping_files(clipping_files)
    log_info(f"Contigs with misassemblies: {len(clipping_positions):,}")
    
    # 3. Classify all contigs
    categories = classify_contigs(assembly, hq_ids, clipping_positions)
    
    # 4. Print classification summary
    print_section("CONTIG CLASSIFICATION")
    print(f"  HQ genomes (clean):              {len(categories['hq_clean']):>8,}")
    print(f"  HQ genomes (with misassemblies): {len(categories['hq_problematic']):>8,}")
    print(f"  Non-HQ (with misassemblies):     {len(categories['nonhq_misassembly']):>8,}")
    print(f"  Non-HQ (clean):                  {len(categories['nonhq_clean']):>8,}")
    
    # 5. Decision logic
    print_section("PROCESSING STRATEGY")
    
    to_remove = categories['hq_clean'].copy()
    to_split = categories['nonhq_misassembly'].copy()
    to_keep_intact = categories['nonhq_clean'].copy()
    
    if preserve_hq_with_issues:
        log_info("Preserving HQ genomes with misassemblies")
        to_split.update(categories['hq_problematic'])
        log_info(f"  → {len(categories['hq_problematic']):,} HQ contigs will be split")
    else:
        log_warning("Removing ALL HQ genomes (including problematic ones)")
        to_remove.update(categories['hq_problematic'])
    
    print(f"\n  Actions:")
    print(f"     Remove (HQ clean):        {len(to_remove):>8,}")
    print(f"     Split (with misasm):      {len(to_split):>8,}")
    print(f"     Keep intact (clean):      {len(to_keep_intact):>8,}")
    
    if dry_run:
        log_info("\nDRY-RUN MODE: No output file generated")
        
        if categories['hq_problematic']:
            print("\nExamples of HQ genomes with misassemblies (would be preserved):")
            for i, contig in enumerate(list(categories['hq_problematic'])[:5], 1):
                n_positions = len(clipping_positions.get(contig, []))
                print(f"  {i}. {contig} (clipping at {n_positions} position(s))")
            if len(categories['hq_problematic']) > 5:
                print(f"  ... and {len(categories['hq_problematic']) - 5} more")
        
        return {
            'hq_clean': len(categories['hq_clean']),
            'hq_problematic': len(categories['hq_problematic']),
            'nonhq_misassembly': len(categories['nonhq_misassembly']),
            'nonhq_clean': len(categories['nonhq_clean']),
            'to_remove': len(to_remove),
            'to_split': len(to_split),
            'to_keep': len(to_keep_intact)
        }
    
    # 6. Process assembly
    print_section(f"Processing assembly (mode: {mode.upper()})")
    
    output_sequences = []
    stats = {
        'removed': 0,
        'split': 0,
        'kept_intact': 0,
        'fragments_generated': 0,
        'fragments_kept': 0,
        'original_bp': 0,
        'final_bp': 0
    }
    
    for contig_name, sequence in read_fasta_streaming(assembly):
        contig_len = len(sequence)
        stats['original_bp'] += contig_len
        
        if contig_name in to_remove:
            stats['removed'] += 1
            continue
        
        if contig_name in to_split and mode == 'split':
            split_pos = clipping_positions[contig_name]
            fragments = split_contig(contig_name, sequence, split_pos, min_length)
            
            stats['split'] += 1
            stats['fragments_generated'] += len(split_pos) + 1
            stats['fragments_kept'] += len(fragments)
            
            output_sequences.extend(fragments)
            stats['final_bp'] += sum(len(seq) for _, seq in fragments)
        
        elif contig_name in to_split and mode == 'remove':
            stats['removed'] += 1
            continue
        
        else:
            if contig_len >= min_length:
                output_sequences.append((contig_name, sequence))
                stats['kept_intact'] += 1
                stats['final_bp'] += contig_len
            else:
                stats['removed'] += 1
    
    # 7. Write output
    log_info(f"\nWriting output to: {output}")
    write_fasta(output_sequences, output)
    
    # 8. Final statistics
    print_section("FINAL STATISTICS")
    print(f"  Contigs removed:         {stats['removed']:>10,}")
    print(f"  Contigs split:           {stats['split']:>10,}")
    print(f"  Contigs kept intact:     {stats['kept_intact']:>10,}")
    print(f"  Total output contigs:    {len(output_sequences):>10,}")
    print(f"\n  Original bases:          {stats['original_bp']:>10,} bp")
    print(f"  Final bases:             {stats['final_bp']:>10,} bp")
    pct = 100 * stats['final_bp'] / stats['original_bp'] if stats['original_bp'] > 0 else 0
    print(f"  Retained:                {pct:>9.2f}%")
    
    if categories['hq_problematic'] and preserve_hq_with_issues:
        print(f"\n  ⚠  {len(categories['hq_problematic']):,} HQ genomes had misassemblies")
        print(f"     and were PRESERVED (split into fragments)")
    
    return stats


def register_parser(subparsers):
    """Register this command's parser"""
    import argparse
    
    parser = subparsers.add_parser(
        'filter-misassemblies',
        help='Filter assembly based on misassembly detection',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Dry-run to see what would happen
  mamisa filter-misassemblies \\
    --assembly assembly.fa \\
    --misassemblies misasm_dir/ \\
    --hq-genomes hq_dir/ \\
    --dry-run
  
  # Process with HQ preservation
  mamisa filter-misassemblies \\
    --assembly assembly.fa \\
    --misassemblies misasm_dir/ \\
    --hq-genomes hq_dir/ \\
    --output clean.fa \\
    --mode split \\
    --min-length 2500 \\
    --preserve-hq-with-issues
        """
    )
    
    parser.add_argument('-a', '--assembly', type=Path, required=True,
                        help='Input assembly FASTA file')
    parser.add_argument('-m', '--misassemblies', type=Path, required=True,
                        help='Directory with *-clipping.txt files')
    parser.add_argument('-g', '--hq-genomes', type=Path,
                        help='Directory with HQ genome files (optional)')
    parser.add_argument('-o', '--output', type=Path,
                        help='Output filtered assembly file')
    parser.add_argument('-l', '--min-length', type=int, default=2500,
                        help='Minimum contig/fragment length (default: 2500)')
    parser.add_argument('--mode', choices=['remove', 'split'], default='split',
                        help='Misassembly handling mode (default: split)')
    parser.add_argument('--preserve-hq-with-issues', action='store_true',
                        help='Keep HQ genomes with misassemblies')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show statistics without generating output')
    parser.add_argument('--stats', type=Path,
                        help='Save statistics to TSV file')
    
    parser.set_defaults(func=run)
    return parser


def run(args):
    """Execute the filter-misassemblies command"""
    
    # Validation
    if not args.dry_run and not args.output:
        raise ValueError("--output is required unless using --dry-run")
    
    validate_file_exists(args.assembly, "Assembly file")
    validate_dir_exists(args.misassemblies, "Misassemblies directory")
    
    # Process
    stats = process_assembly(
        assembly=args.assembly,
        hq_dir=args.hq_genomes,
        misassembly_dir=args.misassemblies,
        output=args.output,
        mode=args.mode,
        min_length=args.min_length,
        preserve_hq_with_issues=args.preserve_hq_with_issues,
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
