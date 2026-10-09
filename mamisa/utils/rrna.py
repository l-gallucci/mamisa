"""
rRNA / tRNA parsing utilities for the MIMAG high-quality criterion.

MIMAG (Bowers et al. 2017) requires, for a High-Quality MAG, in addition to
completeness >90% / contamination <5%:
  - presence of the 5S, 16S and 23S rRNA genes, and
  - tRNAs for at least 18 of the 20 standard amino acids.

These functions parse the output of:
  - barrnap        (rRNA genes, GFF3)
  - tRNAscan-SE    (tRNA genes, tabular -o output)
and combine them into a per-genome MIMAG rRNA/tRNA verdict. They are pure
parsers (no external calls), so they can be unit-tested without the tools.
"""

import csv
import re
from pathlib import Path
from typing import Dict, Set, Tuple

from .logging import log_warning


# The 20 standard amino acids tRNAscan-SE reports in its "Type" column.
STANDARD_AA = {
    'Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile',
    'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val',
}

# tRNA "Type" values that are NOT one of the 20 standard amino acids.
NON_STANDARD_AA = {'Undet', 'Pseudo', 'Sup', 'SeC', 'fMet', 'Ile2'}

MIMAG_MIN_AA = 18


def parse_barrnap_gff(gff_path: Path) -> Dict[str, int]:
    """
    Parse a barrnap GFF3 file and count 5S / 16S / 23S rRNA gene hits.

    Returns a dict with keys '5S', '16S', '23S' (counts). Partial hits (barrnap
    annotates these with a 'note=aligned only ...' / '(partial)' product) are
    NOT counted, matching MIMAG's requirement for the full-length gene.
    """
    counts = {'5S': 0, '16S': 0, '23S': 0}
    if not Path(gff_path).exists():
        return counts

    with open(gff_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'rRNA':
                continue
            attrs = fields[8]
            # skip partial / fragmentary rRNA genes
            if 'partial' in attrs.lower():
                continue
            m = re.search(r'Name=(\d+S)_rRNA', attrs)
            if not m:
                continue
            name = m.group(1)
            if name in counts:
                counts[name] += 1
    return counts


def parse_trnascan(trnascan_path: Path) -> Set[str]:
    """
    Parse a tRNAscan-SE tabular output file (the `-o` output) and return the
    SET of distinct standard amino acids with at least one non-pseudogene tRNA.

    tRNAscan-SE output has three header lines, then columns:
      Name  tRNA#  Begin  End  Type  Codon  ...  Score  Note
    The amino acid is the 5th column (index 4). Pseudogenes (Note contains
    'pseudo') and non-standard types are excluded.
    """
    found: Set[str] = set()
    if not Path(trnascan_path).exists():
        return found

    with open(trnascan_path) as f:
        for line in f:
            fields = line.split('\t') if '\t' in line else line.split()
            if len(fields) < 5:
                continue
            aa = fields[4].strip()
            if aa not in STANDARD_AA:
                continue
            note = fields[-1].strip().lower() if fields else ''
            if 'pseudo' in note:
                continue
            found.add(aa)
    return found


def assess_mimag_rna(rrna_counts: Dict[str, int], trna_aas: Set[str]) -> Dict:
    """
    Combine rRNA counts and tRNA amino-acid set into a MIMAG rRNA/tRNA verdict.

    Returns a dict with:
      n_5S, n_16S, n_23S, n_trna_aa, rrna_complete, trna_complete, mimag_rna_ok
    """
    n_5s = rrna_counts.get('5S', 0)
    n_16s = rrna_counts.get('16S', 0)
    n_23s = rrna_counts.get('23S', 0)
    n_aa = len(trna_aas & STANDARD_AA)

    rrna_complete = n_5s >= 1 and n_16s >= 1 and n_23s >= 1
    trna_complete = n_aa >= MIMAG_MIN_AA

    return {
        'n_5S': n_5s,
        'n_16S': n_16s,
        'n_23S': n_23s,
        'n_trna_aa': n_aa,
        'rrna_complete': rrna_complete,
        'trna_complete': trna_complete,
        'mimag_rna_ok': rrna_complete and trna_complete,
    }


SUMMARY_FIELDS = [
    'genome', 'n_5S', 'n_16S', 'n_23S', 'n_trna_aa',
    'rrna_complete', 'trna_complete', 'mimag_rna_ok',
]


def write_mimag_rna_summary(records: list, out_path: Path) -> None:
    """Write per-genome MIMAG rRNA/tRNA records to a TSV."""
    with open(out_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=SUMMARY_FIELDS, delimiter='\t')
        writer.writeheader()
        for rec in records:
            writer.writerow({k: rec.get(k, '') for k in SUMMARY_FIELDS})


def _as_bool(value: str) -> bool:
    return str(value).strip().lower() in ('1', 'true', 'yes', 't')


def parse_mimag_rna_summary(path: Path) -> Dict[str, Dict]:
    """
    Read a mimag_rna_summary.tsv (written by run-mimag-rna) back into
    {genome_name: record}. Used by organize-mags to gate the HQ tier.
    """
    result: Dict[str, Dict] = {}
    if not Path(path).exists():
        log_warning(f"MIMAG rRNA/tRNA summary not found: {path}")
        return result

    with open(path, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not reader.fieldnames or 'genome' not in reader.fieldnames:
            log_warning(f"Unexpected MIMAG rRNA/tRNA summary format: {path}")
            return result
        for row in reader:
            name = (row.get('genome') or '').strip()
            if not name:
                continue
            result[name] = {
                'n_5S': int(row.get('n_5S') or 0),
                'n_16S': int(row.get('n_16S') or 0),
                'n_23S': int(row.get('n_23S') or 0),
                'n_trna_aa': int(row.get('n_trna_aa') or 0),
                'rrna_complete': _as_bool(row.get('rrna_complete', '')),
                'trna_complete': _as_bool(row.get('trna_complete', '')),
                'mimag_rna_ok': _as_bool(row.get('mimag_rna_ok', '')),
            }
    return result
