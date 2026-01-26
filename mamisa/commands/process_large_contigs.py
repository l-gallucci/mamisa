#!/usr/bin/env python3
"""
MaMISA - process-large-contigs command
Extract and QC large contigs, then intelligently filter based on quality and misassemblies
"""

import sys
import argparse
import subprocess
import shutil
from pathlib import Path
from typing import Set, Dict, List, Tuple
import csv

from ..utils.fasta import read_fasta_streaming, write_fasta
from ..utils.validation import validate_file_exists, validate_dir_exists, check_dependencies
from ..utils.logging import log_info, log_error, log_warning, print_header, print_section


def extract_large_contigs(assembly: Path, max_length: int, output_dir: Path) -> Tuple[int, int]:
    """
    Separate contigs into large (>max_length) and regular
    
    Returns:
        Tuple of (n_large, n_regular)
    """
    large_contigs = []
    regular_contigs = []
    
    for name, seq in read_fasta_streaming(assembly):
        if len(seq) > max_length:
            large_contigs.append((name, seq))
        else:
            regular_contigs.append((name, seq))
    
    # Write files
    large_file = output_dir / "large_contigs.fa"
    regular_file = output_dir / "assembly_regular.fa"
    
    if large_contigs:
        write_fasta(large_contigs, large_file)
        log_info(f"Extracted {len(large_contigs):,} large contigs to: {large_file}")
    else:
        log_warning("No large contigs found")
    
    write_fasta(regular_contigs, regular_file)
    log_info(f"Regular contigs: {len(regular_contigs):,} → {regular_file}")
    
    return len(large_contigs), len(regular_contigs)


def split_large_contigs_individually(large_contigs_file: Path, output_dir: Path) -> int:
    """
    Split large contigs into individual files for CheckM2
    
    Returns:
        Number of files created
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    count = 0
    for name, seq in read_fasta_streaming(large_contigs_file):
        # Sanitize filename
        safe_name = name.replace('/', '_').replace('\\', '_')
        out_file = output_dir / f"{safe_name}.fa"
        write_fasta([(name, seq)], out_file)
        count += 1
    
    log_info(f"Split into {count:,} individual files")
    return count


def run_checkm2(genome_dir: Path, output_dir: Path, threads: int = 1,
                checkm2_env: str = None) -> bool:
    """
    Run CheckM2 on genomes directory
    
    Returns:
        True if successful
    """
    # Check if checkm2 is available
    deps = check_dependencies(['checkm2'])
    
    if deps['checkm2'] is None:
        log_error("CheckM2 not found in PATH")
        if checkm2_env:
            log_error(f"Please activate conda environment: conda activate {checkm2_env}")
        return False
    
    log_info(f"Found CheckM2: {deps['checkm2']}")
    
    # Run CheckM2
    cmd = [
        'checkm2', 'predict',
        '--threads', str(threads),
        '--input', str(genome_dir),
        '--output-directory', str(output_dir),
        '-x', 'fa',
        '--force'
    ]
    
    log_info(f"Running: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        log_info("CheckM2 completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        log_error(f"CheckM2 failed: {e}")
        log_error(f"stderr: {e.stderr}")
        return False


def parse_checkm2_results(checkm2_output_dir: Path) -> Dict[str, Dict]:
    """
    Parse CheckM2 quality_report.tsv
    
    Returns:
        Dict[contig_name] = {'completeness': float, 'contamination': float, 'passes_qc': bool}
    """
    quality_file = checkm2_output_dir / "quality_report.tsv"
    
    if not quality_file.exists():
        log_error(f"CheckM2 quality report not found: {quality_file}")
        return {}
    
    results = {}
    
    with open(quality_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            name = row['Name']
            comp = float(row['Completeness'].replace('%', ''))
            cont = float(row['Contamination'].replace('%', ''))
            
            results[name] = {
                'completeness': comp,
                'contamination': cont
            }
    
    log_info(f"Parsed quality data for {len(results):,} genomes")
    return results


def parse_misassembly_data(misassembly_dir: Path) -> Set[str]:
    """
    Parse anvi'o clipping files to get contigs with misassemblies
    
    Returns:
        Set of contig names with misassemblies
    """
    contigs_with_misasm = set()
    
    clipping_files = list(misassembly_dir.glob('*-clipping.txt'))
    
    if not clipping_files:
        log_warning(f"No *-clipping.txt files found in {misassembly_dir}")
        return contigs_with_misasm
    
    for clip_file in clipping_files:
        with open(clip_file) as f:
            next(f)  # Skip header
            for line in f:
                if line.strip():
                    contig_name = line.split('\t')[0]
                    contigs_with_misasm.add(contig_name)
    
    log_info(f"Found {len(contigs_with_misasm):,} contigs with misassemblies")
    return contigs_with_misasm


def make_filtering_decisions(quality_data: Dict[str, Dict], 
                            misasm_contigs: Set[str],
                            min_completeness: float,
                            max_contamination: float) -> Dict[str, List[str]]:
    """
    Decide what to do with each large contig
    
    Returns:
        Dict with categories: 'extract_hq', 'keep_hq_misasm', 'keep_low_quality'
    """
    decisions = {
        'extract_hq': [],           # HQ without misassemblies → extract
        'keep_hq_misasm': [],       # HQ with misassemblies → keep for splitting
        'keep_low_quality': []      # Low quality → keep
    }
    
    print_section("LARGE CONTIGS ANALYSIS")
    
    for contig_name, qc_data in quality_data.items():
        comp = qc_data['completeness']
        cont = qc_data['contamination']
        passes_qc = (comp >= min_completeness and cont <= max_contamination)
        has_misasm = contig_name in misasm_contigs
        
        print(f"\n{contig_name}:")
        print(f"  Completeness:  {comp:>6.1f}%")
        print(f"  Contamination: {cont:>6.1f}%")
        print(f"  Misassemblies: {'YES' if has_misasm else 'NO'}")
        
        if passes_qc and not has_misasm:
            decisions['extract_hq'].append(contig_name)
            print(f"  → EXTRACT (HQ clean genome)")
        elif passes_qc and has_misasm:
            decisions['keep_hq_misasm'].append(contig_name)
            print(f"  → KEEP (HQ with misassemblies - will be split)")
        else:
            decisions['keep_low_quality'].append(contig_name)
            print(f"  → KEEP (low quality)")
    
    # Summary
    print_section("DECISION SUMMARY")
    print(f"  HQ to extract:               {len(decisions['extract_hq']):>6,}")
    print(f"  HQ with misasm (keep):       {len(decisions['keep_hq_misasm']):>6,}")
    print(f"  Low quality (keep):          {len(decisions['keep_low_quality']):>6,}")
    
    return decisions


def save_decisions(decisions: Dict[str, List[str]], quality_data: Dict[str, Dict],
                  misasm_contigs: Set[str], output_file: Path):
    """Save filtering decisions to TSV file"""
    
    with open(output_file, 'w') as f:
        f.write("contig\tdecision\tcompleteness\tcontamination\thas_misassembly\n")
        
        for contig in decisions['extract_hq']:
            qc = quality_data[contig]
            has_misasm = contig in misasm_contigs
            f.write(f"{contig}\tEXTRACT_HQ\t{qc['completeness']:.2f}\t{qc['contamination']:.2f}\t{has_misasm}\n")
        
        for contig in decisions['keep_hq_misasm']:
            qc = quality_data[contig]
            has_misasm = contig in misasm_contigs
            f.write(f"{contig}\tKEEP_HQ_MISASM\t{qc['completeness']:.2f}\t{qc['contamination']:.2f}\t{has_misasm}\n")
        
        for contig in decisions['keep_low_quality']:
            qc = quality_data[contig]
            has_misasm = contig in misasm_contigs
            f.write(f"{contig}\tKEEP_LOW_QUAL\t{qc['completeness']:.2f}\t{qc['contamination']:.2f}\t{has_misasm}\n")
    
    log_info(f"Decisions saved to: {output_file}")


def extract_hq_genomes(large_contigs_individual_dir: Path, hq_list: List[str],
                      output_dir: Path) -> int:
    """Extract HQ genomes to separate directory"""
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    extracted = 0
    for contig_name in hq_list:
        # Sanitize filename
        safe_name = contig_name.replace('/', '_').replace('\\', '_')
        source = large_contigs_individual_dir / f"{safe_name}.fa"
        dest = output_dir / f"{safe_name}.fa"
        
        if source.exists():
            shutil.copy2(source, dest)
            extracted += 1
        else:
            log_warning(f"Source file not found for {contig_name}: {source}")
    
    log_info(f"Extracted {extracted:,} HQ genomes to: {output_dir}")
    return extracted


def create_updated_assembly(regular_assembly: Path, large_contigs_individual_dir: Path,
                           contigs_to_keep: List[str], output_file: Path):
    """Create updated assembly with regular + selected large contigs"""
    
    sequences = []
    
    # Add all regular contigs
    for name, seq in read_fasta_streaming(regular_assembly):
        sequences.append((name, seq))
    
    log_info(f"Added {len(sequences):,} regular contigs")
    
    # Add selected large contigs
    added = 0
    for contig_name in contigs_to_keep:
        safe_name = contig_name.replace('/', '_').replace('\\', '_')
        contig_file = large_contigs_individual_dir / f"{safe_name}.fa"
        
        if contig_file.exists():
            for name, seq in read_fasta_streaming(contig_file):
                sequences.append((name, seq))
                added += 1
        else:
            log_warning(f"Contig file not found: {contig_file}")
    
    log_info(f"Added {added:,} large contigs to keep")
    
    # Write updated assembly
    write_fasta(sequences, output_file)
    log_info(f"Updated assembly: {len(sequences):,} contigs → {output_file}")


def register_parser(subparsers):
    """Register this command's parser"""
    
    parser = subparsers.add_parser(
        'process-large-contigs',
        help='Extract, QC, and filter large contigs intelligently',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Complete workflow
  mamisa process-large-contigs \\
    --assembly assembly.fa \\
    --misassemblies misasm_dir/ \\
    --output-dir filtered/ \\
    --max-length 300000 \\
    --min-completeness 50 \\
    --max-contamination 10 \\
    --threads 40
  
  # Skip CheckM2 if already run
  mamisa process-large-contigs \\
    --assembly assembly.fa \\
    --misassemblies misasm_dir/ \\
    --output-dir filtered/ \\
    --skip-checkm2
        """
    )
    
    # Input/Output
    parser.add_argument('-a', '--assembly', type=Path, required=True,
                        help='Input assembly FASTA file')
    parser.add_argument('-m', '--misassemblies', type=Path, required=True,
                        help='Directory with anvi\'o *-clipping.txt files')
    parser.add_argument('-o', '--output-dir', type=Path, required=True,
                        help='Output directory')
    
    # Thresholds
    parser.add_argument('--max-length', type=int, default=300000,
                        help='Maximum contig length for extraction (default: 300000 bp)')
    parser.add_argument('--min-completeness', type=float, default=50.0,
                        help='Minimum completeness for HQ (default: 50%%)')
    parser.add_argument('--max-contamination', type=float, default=10.0,
                        help='Maximum contamination for HQ (default: 10%%)')
    
    # CheckM2 options
    parser.add_argument('--skip-checkm2', action='store_true',
                        help='Skip CheckM2 step (use existing results)')
    parser.add_argument('--checkm2-results', type=Path,
                        help='Existing CheckM2 output directory (if --skip-checkm2)')
    parser.add_argument('--threads', type=int, default=1,
                        help='Number of threads for CheckM2 (default: 1)')
    parser.add_argument('--checkm2-env', default='checkm2',
                        help='CheckM2 conda environment name (default: checkm2)')
    
    parser.set_defaults(func=run)
    return parser


def run(args):
    """Execute the process-large-contigs command"""
    
    # Validation
    validate_file_exists(args.assembly, "Assembly file")
    validate_dir_exists(args.misassemblies, "Misassemblies directory")
    
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    print_header("MaMISA - Process Large Contigs")
    
    # Step 1: Extract large contigs
    print_section("STEP 1: Extracting Large Contigs")
    
    extract_dir = args.output_dir / "01_extracted"
    extract_dir.mkdir(exist_ok=True)
    
    n_large, n_regular = extract_large_contigs(
        args.assembly,
        args.max_length,
        extract_dir
    )
    
    if n_large == 0:
        log_warning("No large contigs found - nothing to process")
        log_info("Assembly can be used as-is for binning")
        return
    
    # Step 2: Split large contigs
    print_section("STEP 2: Splitting Large Contigs")
    
    individual_dir = args.output_dir / "02_individual"
    large_contigs_file = extract_dir / "large_contigs.fa"
    
    n_files = split_large_contigs_individually(large_contigs_file, individual_dir)
    
    # Step 3: Run CheckM2
    print_section("STEP 3: Quality Assessment (CheckM2)")
    
    checkm2_dir = args.output_dir / "03_checkm2"
    
    if args.skip_checkm2:
        if args.checkm2_results:
            checkm2_dir = args.checkm2_results
            log_info(f"Using existing CheckM2 results: {checkm2_dir}")
        else:
            log_error("--skip-checkm2 requires --checkm2-results")
            sys.exit(1)
    else:
        log_info(f"Running CheckM2 with {args.threads} threads...")
        success = run_checkm2(individual_dir, checkm2_dir, args.threads, args.checkm2_env)
        
        if not success:
            log_error("CheckM2 failed")
            sys.exit(1)
    
    # Step 4: Parse results
    print_section("STEP 4: Analyzing Results")
    
    quality_data = parse_checkm2_results(checkm2_dir)
    
    if not quality_data:
        log_error("No quality data found")
        sys.exit(1)
    
    misasm_contigs = parse_misassembly_data(args.misassemblies)
    
    # Step 5: Make decisions
    decisions = make_filtering_decisions(
        quality_data,
        misasm_contigs,
        args.min_completeness,
        args.max_contamination
    )
    
    # Step 6: Save decisions
    print_section("STEP 5: Saving Results")
    
    decisions_file = args.output_dir / "filtering_decisions.tsv"
    save_decisions(decisions, quality_data, misasm_contigs, decisions_file)
    
    # Step 7: Extract HQ genomes
    hq_dir = args.output_dir / "HQ_extracted"
    n_extracted = extract_hq_genomes(individual_dir, decisions['extract_hq'], hq_dir)
    
    # Step 8: Create updated assembly
    contigs_to_keep = decisions['keep_hq_misasm'] + decisions['keep_low_quality']
    
    updated_assembly = args.output_dir / "assembly_for_filtering.fa"
    create_updated_assembly(
        extract_dir / "assembly_regular.fa",
        individual_dir,
        contigs_to_keep,
        updated_assembly
    )
    
    # Final summary
    print_section("RESULTS")
    print(f"  Large contigs processed:     {n_large:>10,}")
    print(f"  HQ genomes extracted:        {n_extracted:>10,}")
    print(f"  Contigs in updated assembly: {n_regular + len(contigs_to_keep):>10,}")
    print(f"\n  Decisions:         {decisions_file}")
    print(f"  HQ genomes:        {hq_dir}/")
    print(f"  Updated assembly:  {updated_assembly}")
    
    print("\nNext step:")
    print("  Run mamisa filter-misassemblies on the updated assembly:")
    print(f"  mamisa filter-misassemblies \\")
    print(f"    --assembly {updated_assembly} \\")
    print(f"    --misassemblies {args.misassemblies} \\")
    print(f"    --hq-genomes {hq_dir} \\")
    print(f"    --output final_clean_assembly.fa \\")
    print(f"    --mode split \\")
    print(f"    --preserve-hq-with-issues")
    
    log_info("\n✓ Done!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    register_parser(parser.add_subparsers(dest='command'))
    args = parser.parse_args()
    run(args)
