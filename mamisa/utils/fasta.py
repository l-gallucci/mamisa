"""
FASTA file handling utilities
"""

import re
import gzip
from pathlib import Path
from typing import Tuple, List, Iterator


def is_gzipped(filepath: Path) -> bool:
    """Check if file is gzip compressed."""
    return Path(filepath).suffix.lower() in {'.gz', '.gzip'}


def open_file(filepath: Path, mode: str = 'r'):
    """Open file, handling gzipped files automatically."""
    if is_gzipped(filepath):
        if 'b' not in mode and 't' not in mode:
            mode = mode + 't'
        return gzip.open(filepath, mode, encoding='utf-8')
    return open(filepath, mode)


def extract_contig_id(text: str) -> str:
    """
    Extract a numeric contig identifier from various naming conventions:
      - contig_123              → '123'
      - stdin.part_contig_7420  → '7420'
      - NODE_456_length_...     → '456'
      - scaffold_789            → '789'

    Returns the first integer found in the string, or the stripped original
    text if no integer is present.
    """
    match = re.search(r'(\d+)', text)
    return match.group(1) if match else text.strip()


def read_fasta_streaming(fasta_path: Path) -> Iterator[Tuple[str, str]]:
    """
    Generator yielding (name, sequence) tuples one record at a time.
    Handles both plain and gzipped FASTA files.
    Only the first word of the header line is used as the name.
    """
    current_name = None
    current_seq: List[str] = []

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
    """Read entire FASTA file into memory as {name: sequence} dict."""
    return {name: seq for name, seq in read_fasta_streaming(fasta_path)}


def write_fasta(sequences: List[Tuple[str, str]], output_path: Path,
                line_width: int = 80):
    """
    Write sequences to a FASTA file.
    Automatically writes gzip-compressed output when output_path ends in .gz or .gzip.
    """
    if is_gzipped(output_path):
        ctx = gzip.open(output_path, 'wt', encoding='utf-8')
    else:
        ctx = open(output_path, 'w')

    with ctx as f:
        for name, seq in sequences:
            f.write(f">{name}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i + line_width] + '\n')


def count_sequences(fasta_path: Path) -> Tuple[int, int]:
    """
    Count sequences and total bases in a FASTA file.
    Returns (n_sequences, total_bases).
    """
    n_seqs = 0
    n_bases = 0
    for _, seq in read_fasta_streaming(fasta_path):
        n_seqs += 1
        n_bases += len(seq)
    return n_seqs, n_bases


def get_sequence_lengths(fasta_path: Path) -> dict:
    """Return {name: length} dict for all sequences in a FASTA file."""
    return {name: len(seq) for name, seq in read_fasta_streaming(fasta_path)}


def filter_by_length(input_fasta: Path, output_fasta: Path,
                     min_length: int = 0, max_length: int = float('inf')) -> int:
    """
    Filter a FASTA file by sequence length and write to output_fasta.
    Returns the number of sequences kept.
    """
    kept = [
        (name, seq)
        for name, seq in read_fasta_streaming(input_fasta)
        if min_length <= len(seq) <= max_length
    ]
    write_fasta(kept, output_fasta)
    return len(kept)
