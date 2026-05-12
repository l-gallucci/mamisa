"""
CheckM2 result parsing utilities
"""

import csv
import sys
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple

from .logging import log_info, log_warning, log_error


def parse_quality_value(value: str) -> float:
    """Safely parse a completeness or contamination value to float."""
    if value is None:
        return float('nan')
    try:
        return float(value.strip().replace('%', ''))
    except ValueError:
        return float('nan')


def assign_tier(completeness: float, contamination: float,
                hq_comp: float = 90.0, hq_cont: float = 5.0,
                mq_comp: float = 70.0, mq_cont: float = 10.0,
                lq_comp: float = 50.0, lq_cont: float = 10.0) -> str:
    """
    Assign MiMAG quality tier based on completeness and contamination.
    Returns: 'HQ', 'MQ', 'LQ', or 'Fail'
    """
    if completeness >= hq_comp and contamination <= hq_cont:
        return "HQ"
    if completeness >= mq_comp and contamination <= mq_cont:
        return "MQ"
    if completeness >= lq_comp and contamination <= lq_cont:
        return "LQ"
    return "Fail"


def parse_quality_report(report_path: Path) -> List[Dict]:
    """
    Parse a single CheckM2 quality_report.tsv.

    Returns:
        List of dicts with keys: name, completeness, contamination
    """
    records = []

    with open(report_path, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        required = {'Name', 'Completeness', 'Contamination'}
        if not required.issubset(set(reader.fieldnames or [])):
            log_warning(f"Missing required columns in {report_path}, skipping")
            return records

        for row in reader:
            name = (row.get('Name') or '').strip()
            if not name:
                continue
            records.append({
                'name': name,
                'completeness': parse_quality_value(row.get('Completeness')),
                'contamination': parse_quality_value(row.get('Contamination')),
            })

    return records


def find_quality_reports(root_dir: Path) -> List[Path]:
    """Recursively find all quality_report.tsv files under root_dir."""
    return sorted(root_dir.rglob("quality_report.tsv"))


def parse_all_reports(root_dir: Path, thresholds: Dict) -> Tuple[List[Dict], Dict]:
    """
    Parse all CheckM2 reports under root_dir, assign tiers, and count by tier.

    Args:
        root_dir:   Root directory to search for quality_report.tsv files
        thresholds: Dict with keys hq_comp, hq_cont, mq_comp, mq_cont, lq_comp, lq_cont

    Returns:
        (all_records, tier_counts)
        Each record has keys: report_path, name, completeness, contamination, tier
    """
    reports = find_quality_reports(root_dir)

    if not reports:
        log_error(f"No quality_report.tsv files found in {root_dir}")
        sys.exit(1)

    log_info(f"Found {len(reports):,} CheckM2 report(s)")

    all_records = []
    tier_counts: Dict[str, int] = defaultdict(int)

    for report_path in reports:
        log_info(f"  Processing: {report_path}")
        for rec in parse_quality_report(report_path):
            rec['tier'] = assign_tier(
                rec['completeness'], rec['contamination'],
                thresholds.get('hq_comp', 90.0), thresholds.get('hq_cont', 5.0),
                thresholds.get('mq_comp', 70.0), thresholds.get('mq_cont', 10.0),
                thresholds.get('lq_comp', 50.0), thresholds.get('lq_cont', 10.0),
            )
            rec['report_path'] = str(report_path)
            all_records.append(rec)
            tier_counts[rec['tier']] += 1

    return all_records, dict(tier_counts)
