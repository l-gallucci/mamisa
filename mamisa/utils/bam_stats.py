"""
BAM-level statistics for clipping position classification.

Alignment categories used
-------------------------
Primary only (FLAG & 0x900 == 0)
    Used for: discordant-pair fraction, insert-size distribution, strand ratio.
    Pair-level fields (TLEN, RNEXT, proper-pair flag 0x2) are only reliable here.

Primary + secondary (FLAG & 0x800 == 0, supplementary excluded)
    Used for: local read depth around clipping positions.
    Including secondary (multi-mappers) inflates coverage at repeat-collapsed
    loci — which is precisely the signal we want to detect.

Supplementary (FLAG 0x800) are excluded throughout.
    They require SA-tag parsing, are absent from many metagenomic BAMs, and
    add complexity that is better handled by dedicated SV callers.

All parsing is done via ``samtools view`` subprocess — no pysam required.
"""

import math
import subprocess
import bisect
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from .logging import log_info, log_warning, log_error
from .validation import check_dependencies


# ──────────────────────────────────────────────────────────────────────────────
# Dependency check
# ──────────────────────────────────────────────────────────────────────────────

def check_samtools() -> bool:
    if check_dependencies(['samtools'])['samtools'] is None:
        log_error("samtools not found in PATH — required for BAM stats")
        log_error("Install: conda install -c bioconda samtools")
        return False
    return True


# ──────────────────────────────────────────────────────────────────────────────
# Contig-level mean depth  (samtools coverage)
# ──────────────────────────────────────────────────────────────────────────────

def get_contig_mean_depths(bam_file: Path, min_mapq: int = 20) -> Dict[str, float]:
    """
    Run ``samtools coverage`` and return {contig: mean_depth}.

    samtools coverage already excludes supplementary alignments internally
    (it counts only primary alignments by default).
    We pass --min-MQ to apply the same MAPQ filter used in the rest of the
    pipeline so all coverage estimates are on a consistent base.
    """
    cmd = ['samtools', 'coverage', '--min-MQ', str(min_mapq), str(bam_file)]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        log_error(f"samtools coverage failed: {e.stderr[:300]}")
        return {}

    depths: Dict[str, float] = {}
    for line in result.stdout.splitlines():
        if line.startswith('#') or not line.strip():
            continue
        fields = line.split('\t')
        if len(fields) < 7:
            continue
        contig = fields[0]
        try:
            mean_depth = float(fields[6])   # column 7 = meandepth
        except (ValueError, IndexError):
            continue
        depths[contig] = mean_depth

    log_info(f"samtools coverage: mean depths for {len(depths):,} contigs")
    return depths


# ──────────────────────────────────────────────────────────────────────────────
# Library insert-size estimation
# ──────────────────────────────────────────────────────────────────────────────

def estimate_insert_size(bam_file: Path, n_reads: int = 20_000,
                         min_mapq: int = 20) -> Tuple[float, float]:
    """
    Estimate library mean and std of insert sizes from the first ``n_reads``
    primary, properly-paired, concordant reads.

    Returns (mean, std). Returns (0, 0) for single-end or if estimation fails.
    """
    # -f 2 = properly paired; -F 2308 = not unmapped (4), not mate-unmapped (8),
    #        not secondary (256), not supplementary (2048)
    cmd = ['samtools', 'view', '-f', '2', '-F', '2308',
           '-q', str(min_mapq), str(bam_file)]
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.DEVNULL, text=True)
    except Exception as e:
        log_warning(f"Could not estimate insert size: {e}")
        return 0.0, 0.0

    sizes: List[int] = []
    for line in proc.stdout:
        fields = line.split('\t')
        if len(fields) < 9:
            continue
        try:
            tlen = int(fields[8])
        except ValueError:
            continue
        if tlen > 0:
            sizes.append(tlen)
        if len(sizes) >= n_reads:
            break

    proc.stdout.close()
    proc.wait()

    if not sizes:
        return 0.0, 0.0

    mean = sum(sizes) / len(sizes)
    variance = sum((x - mean) ** 2 for x in sizes) / len(sizes)
    std = math.sqrt(variance)
    log_info(f"Library insert size: mean={mean:.0f} bp, std={std:.0f} bp "
             f"(estimated from {len(sizes):,} pairs)")
    return mean, std


# ──────────────────────────────────────────────────────────────────────────────
# Per-position stats via single BAM stream
# ──────────────────────────────────────────────────────────────────────────────

class _PosAccumulator:
    """Collects read-level evidence for one clipping position."""
    __slots__ = (
        'clip_pos',
        # primary + secondary (for depth)
        'n_depth',
        # primary only
        'n_primary',
        'n_proper_pair',
        'n_discordant',
        'n_large_insert',
        'n_forward',
        'n_reverse',
        # clipped-base entropy (from reads soft-clipping at exactly this position)
        '_clipped_bases',
    )

    def __init__(self, clip_pos: int):
        self.clip_pos = clip_pos
        self.n_depth = 0
        self.n_primary = 0
        self.n_proper_pair = 0
        self.n_discordant = 0
        self.n_large_insert = 0
        self.n_forward = 0
        self.n_reverse = 0
        self._clipped_bases: List[str] = []

    def add_read(self, flag: int, tlen: int, seq: str, cigar: str,
                 read_pos: int, insert_mean: float, insert_std: float):
        is_secondary = bool(flag & 0x100)
        is_primary = not is_secondary          # (supplementary already excluded upstream)

        # All non-supplementary reads count toward depth
        self.n_depth += 1

        if not is_primary:
            return

        self.n_primary += 1

        # Proper pair
        if flag & 0x2:
            self.n_proper_pair += 1
        else:
            self.n_discordant += 1

        # Large insert (only for proper pairs to avoid noise)
        if (flag & 0x2) and insert_mean > 0 and abs(tlen) > insert_mean + 3 * insert_std:
            self.n_large_insert += 1

        # Strand
        if flag & 0x10:
            self.n_reverse += 1
        else:
            self.n_forward += 1

        # Collect clipped bases from reads whose soft-clip edge is AT this position.
        # CIGAR 'S' at the start → clip ends at read_pos-1 (position before alignment).
        # CIGAR 'S' at the end  → clip starts at read_pos + aligned_length.
        # We only collect up to 200 clipped bases per position to bound memory.
        if len(self._clipped_bases) < 200 and seq and cigar and cigar != '*':
            clipped = _extract_clipped_bases(cigar, seq, read_pos, self.clip_pos)
            if clipped:
                self._clipped_bases.append(clipped)

    def clipped_base_entropy(self) -> float:
        """Shannon entropy (bits) of the base composition of clipped sequences."""
        if not self._clipped_bases:
            return float('nan')
        combined = ''.join(self._clipped_bases).upper()
        n = len(combined)
        if n == 0:
            return float('nan')
        counts = {b: combined.count(b) for b in 'ACGTN'}
        entropy = 0.0
        for b, c in counts.items():
            if c > 0:
                p = c / n
                entropy -= p * math.log2(p)
        return entropy

    def to_dict(self, contig: str, contig_length: int,
                contig_mean_depth: float) -> Dict:
        depth_ratio = (self.n_depth / contig_mean_depth
                       if contig_mean_depth > 0 else float('nan'))
        clip_frac = float('nan')
        if self.n_depth > 0:
            clip_frac = self.n_depth / self.n_depth   # placeholder — real clip frac below
        disc_frac = (self.n_discordant / self.n_primary
                     if self.n_primary > 0 else float('nan'))
        large_ins_frac = (self.n_large_insert / self.n_primary
                          if self.n_primary > 0 else float('nan'))
        strand_bias = float('nan')
        if self.n_forward + self.n_reverse > 0:
            strand_bias = self.n_forward / (self.n_forward + self.n_reverse)

        near_end = (
            self.clip_pos < 500
            or (contig_length > 0 and self.clip_pos > contig_length - 500)
        )

        return {
            'contig': contig,
            'clip_pos': self.clip_pos,
            'contig_length': contig_length,
            'local_depth': self.n_depth,
            'primary_reads_in_window': self.n_primary,
            'contig_mean_depth': round(contig_mean_depth, 2),
            'depth_ratio': round(depth_ratio, 3),
            'discordant_fraction': round(disc_frac, 3) if not math.isnan(disc_frac) else 'N/A',
            'large_insert_fraction': round(large_ins_frac, 3) if not math.isnan(large_ins_frac) else 'N/A',
            'strand_fwd_fraction': round(strand_bias, 3) if not math.isnan(strand_bias) else 'N/A',
            'clipped_base_entropy': round(self.clipped_base_entropy(), 3)
                                    if not math.isnan(self.clipped_base_entropy()) else 'N/A',
            'near_contig_end': near_end,
        }


def _extract_clipped_bases(cigar: str, seq: str,
                            read_pos: int, clip_pos: int,
                            tol: int = 10) -> str:
    """
    If the read has a soft-clip whose boundary is within `tol` bp of `clip_pos`,
    return the clipped bases; otherwise return ''.
    Uses a simple CIGAR parser — handles S, M, D, I, N, H operations.
    """
    import re
    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    if not ops:
        return ''

    seq_pos = 0       # position in the read sequence
    ref_pos = read_pos

    for length_s, op in ops:
        length = int(length_s)
        if op == 'S':
            # Soft-clip: read_pos before or after alignment
            # The boundary on the reference side is ref_pos
            if abs(ref_pos - clip_pos) <= tol:
                return seq[seq_pos:seq_pos + length]
            seq_pos += length
        elif op in ('M', '=', 'X'):
            seq_pos += length
            ref_pos += length
        elif op == 'D' or op == 'N':
            ref_pos += length
        elif op == 'I':
            seq_pos += length
        elif op == 'H':
            pass   # hard-clip: bases not in seq

    return ''


def collect_position_stats(
    bam_file: Path,
    clipping_data: Dict[str, List],   # contig → [(pos, coverage), ...]
    window: int,
    min_mapq: int,
    contig_mean_depths: Dict[str, float],
    contig_lengths: Dict[str, int],
    insert_mean: float,
    insert_std: float,
) -> Dict[Tuple[str, int], Dict]:
    """
    Single streaming pass through the BAM.
    Returns {(contig, pos): stats_dict} for every clipping position.

    Uses primary + secondary alignments for depth (supplementary excluded via -F 2048).
    Pair statistics (discordant, insert size) use primary only (checked per-read via FLAG).
    """
    if not clipping_data:
        return {}

    # Build fast lookup: contig → sorted list of clip positions
    pos_index: Dict[str, List[int]] = {
        contig: sorted(r[0] for r in records)
        for contig, records in clipping_data.items()
    }
    # Accumulators: (contig, clip_pos) → _PosAccumulator
    accumulators: Dict[Tuple[str, int], _PosAccumulator] = {}
    for contig, records in clipping_data.items():
        for pos, _ in records:
            accumulators[(contig, pos)] = _PosAccumulator(pos)

    # Stream BAM: exclude supplementary (-F 2048), apply MAPQ filter
    cmd = ['samtools', 'view', '-F', '2048', '-q', str(min_mapq), str(bam_file)]
    log_info("Streaming BAM for per-position stats…")
    n_reads = 0

    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True, bufsize=1)
        for line in proc.stdout:
            if line.startswith('@'):
                continue
            fields = line.split('\t')
            if len(fields) < 10:
                continue

            flag = int(fields[1])
            contig = fields[2]
            if contig == '*' or contig not in pos_index:
                continue

            try:
                read_pos = int(fields[3])
            except ValueError:
                continue

            seq = fields[9]
            cigar = fields[5]
            try:
                tlen = int(fields[8])
            except ValueError:
                tlen = 0

            # Find clipping positions within window of this read's start
            positions = pos_index[contig]
            lo = bisect.bisect_left(positions, read_pos - window)
            hi = bisect.bisect_right(positions, read_pos + window)

            for clip_pos in positions[lo:hi]:
                key = (contig, clip_pos)
                if key in accumulators:
                    accumulators[key].add_read(
                        flag, tlen, seq, cigar,
                        read_pos, insert_mean, insert_std,
                    )

            n_reads += 1
            if n_reads % 1_000_000 == 0:
                log_info(f"  Processed {n_reads:,} reads…")

        proc.stdout.close()
        proc.wait()
    except Exception as e:
        log_error(f"BAM streaming failed: {e}")
        return {}

    log_info(f"BAM stream complete: {n_reads:,} reads processed")

    # Build output dicts
    results: Dict[Tuple[str, int], Dict] = {}
    for (contig, clip_pos), acc in accumulators.items():
        mean_depth = contig_mean_depths.get(contig, 0.0)
        contig_len = contig_lengths.get(contig, 0)
        results[(contig, clip_pos)] = acc.to_dict(contig, contig_len, mean_depth)

    return results
