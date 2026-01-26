"""
FASTA file handling utilities
"""

from pathlib import Path
from typing import Tuple, List, Iterator
import gzip


def is_gzipped(filepath: Path) -> bool:
    """Check if file is gzip compressed"""
    return filepath.suffix.lower() in {'.gz', '.gzip'}


def open_file(filepath: Path, mode: str = 'r'):
    """Open file, handling gzipped files automatically"""
    if is_gzipped(filepath):
        if 'b' not in mode and 't' not in mode:
            mode = mode + 't'
        return gzip.open(filepath, mode, encoding='utf-8')
    return open(filepath, mode)


def read_fasta_streaming(fasta_path: Path) -> Iterator[Tuple[str, str]]:
    """
    Generator that yields (name, sequence) tuples without loading all into memory
    Handles both plain and gzipped FASTA files
    
    Args:
        fasta_path: Path to FASTA file
        
    Yields:
        Tuple of (contig_name, sequence)
    """
    current_name = None
    current_seq = []
    
    with open_file(fasta_path, 'rt') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    yield current_name, ''.join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_name:
            yield current_name, ''.join(current_seq)


def read_fasta(fasta_path: Path) -> dict:
    """
    Read entire FASTA file into memory as dictionary
    
    Args:
        fasta_path: Path to FASTA file
        
    Returns:
        Dictionary mapping contig names to sequences
    """
    sequences = {}
    for name, seq in read_fasta_streaming(fasta_path):
        sequences[name] = seq
    return sequences


def write_fasta(sequences: List[Tuple[str, str]], output_path: Path, 
                line_width: int = 80):
    """
    Write sequences to FASTA file
    
    Args:
        sequences: List of (name, sequence) tuples
        output_path: Output file path
        line_width: Number of characters per line (default: 80)
    """
    with open(output_path, 'w') as f:
        for name, seq in sequences:
            f.write(f">{name}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i+line_width] + '\n')


def count_sequences(fasta_path: Path) -> Tuple[int, int]:
    """
    Count sequences and total bases in FASTA file
    
    Args:
        fasta_path: Path to FASTA file
        
    Returns:
        Tuple of (number_of_sequences, total_bases)
    """
    n_seqs = 0
    n_bases = 0
    
    for name, seq in read_fasta_streaming(fasta_path):
        n_seqs += 1
        n_bases += len(seq)
    
    return n_seqs, n_bases


def get_sequence_lengths(fasta_path: Path) -> dict:
    """
    Get lengths of all sequences in FASTA file
    
    Args:
        fasta_path: Path to FASTA file
        
    Returns:
        Dictionary mapping sequence names to lengths
    """
    lengths = {}
    for name, seq in read_fasta_streaming(fasta_path):
        lengths[name] = len(seq)
    return lengths


def filter_by_length(input_fasta: Path, output_fasta: Path, 
                     min_length: int = 0, max_length: int = float('inf')) -> int:
    """
    Filter FASTA file by sequence length
    
    Args:
        input_fasta: Input FASTA file
        output_fasta: Output FASTA file
        min_length: Minimum sequence length (inclusive)
        max_length: Maximum sequence length (inclusive)
        
    Returns:
        Number of sequences kept
    """
    kept = []
    
    for name, seq in read_fasta_streaming(input_fasta):
        seq_len = len(seq)
        if min_length <= seq_len <= max_length:
            kept.append((name, seq))
    
    write_fasta(kept, output_fasta)
    return len(kept)
