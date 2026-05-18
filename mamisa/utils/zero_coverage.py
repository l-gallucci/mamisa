"""
Zero-coverage region extraction and read-support validation.

Workflow
--------
  1. Find zero-coverage intervals in the assembly
       - primary: parse anvi'o *-zero_cov.txt (already computed)
       - fallback: stream samtools depth -a and detect zero runs
  2. Extract sequences for those intervals
  3. BLAST them against the original read set
       Regions with zero hits = assembler-invented sequence with no read backing
  4. (optional) Meryl 21-mer support fraction
       Fraction of k-mers in the region that are present in the source reads.
       Useful cross-check when BLAST finds no hits for compositionally simple seqs.

Reference: Meren Lab long-read assembly benchmarking workflow
  https://merenlab.org/data/benchmarking-long-read-assemblers/
"""

import gzip
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .logging import log_info, log_warning, log_error
from .validation import check_dependencies
from .fasta import read_fasta_streaming


# (contig, start_0based, end_exclusive)
ZeroCovRegion = Tuple[str, int, int]


# ──────────────────────────────────────────────────────────────────────────────
# Dependency checks
# ──────────────────────────────────────────────────────────────────────────────

def check_blast() -> bool:
    deps = check_dependencies(['blastn', 'makeblastdb'])
    missing = [k for k, v in deps.items() if v is None]
    if missing:
        log_error(f"BLAST+ not found ({', '.join(missing)}). "
                  "Install: conda install -c bioconda blast")
        return False
    return True


def check_meryl() -> bool:
    if check_dependencies(['meryl'])['meryl'] is None:
        log_warning("meryl not found — k-mer support check skipped. "
                    "Install: conda install -c bioconda meryl")
        return False
    return True


# ──────────────────────────────────────────────────────────────────────────────
# Source 1: parse anvi'o *-zero_cov.txt files
# ──────────────────────────────────────────────────────────────────────────────

def parse_zero_cov_files(zero_cov_dir: Path) -> List[ZeroCovRegion]:
    """
    Parse *-zero_cov.txt files produced by anvi-script-find-misassemblies.

    Expected format (tab-separated, one header line):
        contig_name  sample_name  start  end
    Coordinates are 1-based; converted to 0-based half-open here.
    """
    regions: List[ZeroCovRegion] = []
    files = sorted(zero_cov_dir.glob('*-zero_cov.txt'))

    if not files:
        log_warning(f"No *-zero_cov.txt files in {zero_cov_dir}")
        return regions

    for fpath in files:
        log_info(f"  Reading: {fpath.name}")
        with open(fpath) as f:
            f.readline()  # header
            for lineno, line in enumerate(f, start=2):
                if not line.strip():
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 4:
                    log_warning(
                        f"Skipping line {lineno} in {fpath.name}: "
                        f"expected ≥4 fields, got {len(fields)}"
                    )
                    continue
                contig = fields[0]
                try:
                    start = int(fields[2]) - 1   # 1-based → 0-based
                    end   = int(fields[3])        # end is already exclusive
                except ValueError:
                    log_warning(f"Bad coordinates on line {lineno} in {fpath.name}")
                    continue
                if end > start:
                    regions.append((contig, start, end))

    log_info(f"Parsed {len(regions):,} zero-coverage regions from {len(files)} file(s)")
    return regions


# ──────────────────────────────────────────────────────────────────────────────
# Source 2: compute zero-coverage intervals from BAM
# ──────────────────────────────────────────────────────────────────────────────

def compute_zero_cov_from_bam(
    bam_file: Path,
    min_mapq: int = 20,
    min_length: int = 500,
) -> List[ZeroCovRegion]:
    """
    Stream ``samtools depth -a`` and collect runs of zero depth.
    Returns (contig, start_0based, end_exclusive) for every run >= min_length.
    """
    cmd = ['samtools', 'depth', '-a', '-Q', str(min_mapq), str(bam_file)]
    regions: List[ZeroCovRegion] = []

    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.DEVNULL, text=True)
    except Exception as e:
        log_error(f"samtools depth failed: {e}")
        return regions

    cur_contig: Optional[str] = None
    run_start: Optional[int] = None
    prev_pos: Optional[int] = None

    def _flush(contig, start, end):
        if end - start + 1 >= min_length:
            regions.append((contig, start - 1, end))   # convert to 0-based

    for line in proc.stdout:
        parts = line.split('\t')
        if len(parts) < 3:
            continue
        contig = parts[0]
        try:
            pos   = int(parts[1])
            depth = int(parts[2])
        except ValueError:
            continue

        if contig != cur_contig:
            if run_start is not None:
                _flush(cur_contig, run_start, prev_pos)
            cur_contig = contig
            run_start  = None
            prev_pos   = None

        if depth == 0:
            if run_start is None:
                run_start = pos
            prev_pos = pos
        else:
            if run_start is not None:
                _flush(contig, run_start, prev_pos)
            run_start = None
            prev_pos  = None

    if run_start is not None:
        _flush(cur_contig, run_start, prev_pos)

    proc.stdout.close()
    proc.wait()

    log_info(f"BAM depth scan: {len(regions):,} zero-coverage regions ≥{min_length} bp")
    return regions


# ──────────────────────────────────────────────────────────────────────────────
# Extract sequences
# ──────────────────────────────────────────────────────────────────────────────

def extract_region_sequences(
    assembly: Path,
    regions: List[ZeroCovRegion],
    min_length: int = 500,
) -> Dict[str, Tuple[str, int, int, str]]:
    """
    Pull subsequences for each region from the assembly.

    Returns {region_id: (contig, start, end, sequence)}.
    region_id = "contig:start-end"  (0-based, half-open).
    Regions shorter than min_length after slicing are skipped.
    """
    by_contig: Dict[str, List[Tuple[int, int]]] = {}
    for contig, start, end in regions:
        if end - start >= min_length:
            by_contig.setdefault(contig, []).append((start, end))

    result: Dict[str, Tuple[str, int, int, str]] = {}
    for contig_name, seq in read_fasta_streaming(assembly):
        if contig_name not in by_contig:
            continue
        for start, end in by_contig[contig_name]:
            subseq = seq[start:end]
            if len(subseq) >= min_length:
                rid = f"{contig_name}:{start}-{end}"
                result[rid] = (contig_name, start, end, subseq)

    log_info(f"Extracted {len(result):,} region sequences for validation")
    return result


def _write_region_fasta(
    region_seqs: Dict[str, Tuple[str, int, int, str]],
    path: Path,
) -> None:
    with open(path, 'w') as f:
        for rid, (_, _, _, seq) in region_seqs.items():
            f.write(f">{rid}\n{seq}\n")


# ──────────────────────────────────────────────────────────────────────────────
# Reads file handling (FASTQ → FASTA conversion)
# ──────────────────────────────────────────────────────────────────────────────

def _is_fastq(reads_file: Path) -> bool:
    name = reads_file.name.lower()
    return any(name.endswith(sfx) for sfx in
               ('.fq', '.fastq', '.fq.gz', '.fastq.gz'))


def _fastq_to_fasta(reads_file: Path, fasta_path: Path) -> None:
    opener = gzip.open if str(reads_file).endswith('.gz') else open
    with opener(reads_file, 'rt') as fin, open(fasta_path, 'w') as fout:
        while True:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline()
            fin.readline()   # +
            fin.readline()   # quality
            if header.startswith('@'):
                fout.write('>' + header[1:])
                fout.write(seq)


def prepare_reads_fasta(reads_file: Path, work_dir: Path) -> Path:
    """Return a FASTA version of the reads (converts FASTQ/FASTQ.gz if needed)."""
    if _is_fastq(reads_file):
        fasta_path = work_dir / 'reads_converted.fasta'
        log_info("Converting reads to FASTA…")
        _fastq_to_fasta(reads_file, fasta_path)
        return fasta_path
    return reads_file


# ──────────────────────────────────────────────────────────────────────────────
# BLAST
# ──────────────────────────────────────────────────────────────────────────────

def make_blast_db(reads_fasta: Path, db_path: Path) -> bool:
    """Build a nucleotide BLAST database from the reads FASTA."""
    cmd = [
        'makeblastdb',
        '-in', str(reads_fasta),
        '-dbtype', 'nucl',
        '-out', str(db_path),
        '-parse_seqids',
    ]
    try:
        subprocess.run(cmd, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        log_info(f"BLAST database built: {db_path}")
        return True
    except subprocess.CalledProcessError as e:
        log_error(f"makeblastdb failed: {e.stderr.decode()[:300]}")
        return False


def blast_regions_vs_reads(
    query_fasta: Path,
    blast_db: Path,
    threads: int = 4,
    evalue: float = 1e-5,
    max_target_seqs: int = 1,
) -> Dict[str, int]:
    """
    Run blastn (dust=no, as in the anvi'o workflow) and return
    {region_id: n_hits}.  Regions absent from the dict have zero hits.
    """
    cmd = [
        'blastn',
        '-query',           str(query_fasta),
        '-db',              str(blast_db),
        '-dust',            'no',
        '-evalue',          str(evalue),
        '-max_target_seqs', str(max_target_seqs),
        '-num_threads',     str(threads),
        '-outfmt',          '6 qseqid',
        '-out',             '-',
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        log_error(f"blastn failed: {e.stderr[:300]}")
        return {}

    hits: Dict[str, int] = {}
    for line in result.stdout.splitlines():
        qid = line.strip()
        if qid:
            hits[qid] = hits.get(qid, 0) + 1

    n_with_hits = len(hits)
    log_info(f"BLAST: {n_with_hits:,} regions with at least one hit")
    return hits


# ──────────────────────────────────────────────────────────────────────────────
# Meryl k-mer support (optional)
# ──────────────────────────────────────────────────────────────────────────────

def _meryl_count(src: Path, db: Path, k: int, threads: int) -> bool:
    try:
        subprocess.run(
            ['meryl', 'count', f'k={k}', f'threads={threads}',
             str(src), 'output', str(db)],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        return True
    except subprocess.CalledProcessError:
        return False


def _meryl_intersect(db_a: Path, db_b: Path, out: Path) -> bool:
    try:
        subprocess.run(
            ['meryl', 'intersect', str(db_a), str(db_b), 'output', str(out)],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        return True
    except subprocess.CalledProcessError:
        return False


def _meryl_count_kmers(db: Path) -> int:
    """Count distinct k-mers in a Meryl database via `meryl print`."""
    try:
        r = subprocess.run(
            ['meryl', 'print', str(db)],
            capture_output=True, text=True, check=True,
        )
        return sum(1 for line in r.stdout.splitlines() if line.strip())
    except subprocess.CalledProcessError:
        return 0


def run_meryl_kmer_support(
    region_seqs: Dict[str, Tuple[str, int, int, str]],
    reads_fasta: Path,
    k: int = 21,
    threads: int = 4,
    work_dir: Optional[Path] = None,
) -> Dict[str, float]:
    """
    Compute the fraction of k-mers in each zero-coverage region that are
    present in the source reads.

    Returns {region_id: support_fraction}  (nan = meryl run failed for this region)
    0.0 = no read k-mers; 1.0 = fully supported.
    """
    import math

    if work_dir is None:
        work_dir = reads_fasta.parent / '_meryl_tmp'
    work_dir.mkdir(parents=True, exist_ok=True)

    # Build reads k-mer database once
    reads_db = work_dir / 'reads.meryl'
    log_info(f"Meryl: building read k-mer database (k={k})…")
    if not _meryl_count(reads_fasta, reads_db, k, threads):
        log_error("meryl count on reads failed")
        return {}

    results: Dict[str, float] = {}

    for rid, (_, _, _, seq) in region_seqs.items():
        region_fa  = work_dir / '_region.fasta'
        region_db  = work_dir / '_region.meryl'
        support_db = work_dir / '_supported.meryl'

        with open(region_fa, 'w') as f:
            f.write(f">{rid}\n{seq}\n")

        if not _meryl_count(region_fa, region_db, k, 1):
            results[rid] = float('nan')
            continue

        total = _meryl_count_kmers(region_db)
        if total == 0:
            results[rid] = float('nan')
            continue

        if not _meryl_intersect(region_db, reads_db, support_db):
            results[rid] = float('nan')
            continue

        supported = _meryl_count_kmers(support_db)
        results[rid] = supported / total

    log_info(f"Meryl: k-mer support computed for {len(results):,} regions")
    return results


# ──────────────────────────────────────────────────────────────────────────────
# Verdict
# ──────────────────────────────────────────────────────────────────────────────

def assign_verdict(
    n_blast_hits: int,
    kmer_support: float,
    use_meryl: bool,
) -> str:
    """
    unsupported  no BLAST hits and (if Meryl) < 10 % k-mer support
    partial      BLAST hit but low k-mer support, or k-mer support 10–50 %
    supported    BLAST hit found; region has read backing
    """
    import math
    nan = float('nan')
    has_kmer = use_meryl and not math.isnan(kmer_support)

    if n_blast_hits == 0:
        if has_kmer and kmer_support >= 0.10:
            return 'partial'
        return 'unsupported'
    else:
        if has_kmer and kmer_support < 0.10:
            return 'partial'
        return 'supported'
