#!/usr/bin/env python3
"""
MaMISA - Manage Misassemblies
Main entry point for package execution
"""

import sys
import argparse
from pathlib import Path

from mamisa.commands import filter_misassemblies
from mamisa.commands import remove_hq_contigs
from mamisa.commands import filter_checkm2
from mamisa.commands import run_gtdbtk
from mamisa.commands import process_large_contigs
from mamisa.commands import check_chimeras
from mamisa.commands import check_read_chimeras
from mamisa.commands import classify_clipping
from mamisa.commands import check_zero_coverage
from mamisa import __version__


def main():
    parser = argparse.ArgumentParser(
        prog='mamisa',
        description='MaMISA - Manage Misassemblies: Toolkit for assembly quality control',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Available commands:
  filter-misassemblies    Filter assembly based on misassembly detection
  remove-hq-contigs       Remove high-quality genome contigs from assembly
  filter-checkm2          Filter genomes based on CheckM2 quality reports
  run-gtdbtk              Run GTDB-Tk taxonomy classification
  process-large-contigs   Extract, QC, and filter large contigs intelligently
  check-chimeras          Detect chimeric MAGs and circular contigs (GC-based)
  check-read-chimeras     Detect chimeric contigs via read-level taxonomy (Kraken2+BAM)
  classify-clipping       Classify each clipping position with BAM evidence
  check-zero-coverage     Validate assembly regions with no read coverage via BLAST

Examples:
  mamisa filter-misassemblies --help
  mamisa filter-checkm2 --help
  
For more information, visit: https://github.com/yourusername/mamisa
        """
    )
    
    parser.add_argument('-v', '--version', action='version', 
                        version=f'MaMISA version {__version__}')
    
    subparsers = parser.add_subparsers(dest='command', help='Command to run')
    
    # Register all subcommands
    filter_misassemblies.register_parser(subparsers)
    remove_hq_contigs.register_parser(subparsers)
    filter_checkm2.register_parser(subparsers)
    run_gtdbtk.register_parser(subparsers)
    process_large_contigs.register_parser(subparsers)
    check_chimeras.register_parser(subparsers)
    check_read_chimeras.register_parser(subparsers)
    classify_clipping.register_parser(subparsers)
    check_zero_coverage.register_parser(subparsers)
    
    # Parse arguments
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        sys.exit(0)
    
    # Execute command
    try:
        args.func(args)
    except KeyboardInterrupt:
        print("\n[!] Interrupted by user", file=sys.stderr)
        sys.exit(130)
    except Exception as e:
        print(f"\n[ERROR] {e}", file=sys.stderr)
        if '--debug' in sys.argv:
            raise
        sys.exit(1)


if __name__ == '__main__':
    main()
