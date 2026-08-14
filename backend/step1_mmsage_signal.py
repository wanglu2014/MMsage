"""
Step 1: MMSage Pseudotime Signal Processing
=============================================
Load coordinates CSV files, parse bacteria/metabolite pairs,
compute normalized MMSage signal, output top-N candidates.
"""

import re
import json
import glob
from pathlib import Path
from typing import Optional

import pandas as pd
import numpy as np


PROJECT_DIR = Path(__file__).resolve().parent.parent


def portable_source_path(filepath: str | Path) -> str:
    """Use a project-relative source reference when the input is packaged locally."""
    source = Path(filepath)
    try:
        return source.resolve().relative_to(PROJECT_DIR).as_posix()
    except ValueError:
        return source.as_posix()


def extract_bacteria_from_filename(filename: str) -> Optional[str]:
    """
    Extract bacteria name from coordinates filename.

    Expected pattern:
      pluscombno1V0317_min0_Akkermansia_muciniphila_seed_1_dim_...
    Bacteria name sits between the third underscore-delimited token
    and '_seed_'.
    """
    base = Path(filename).stem  # remove .csv extension
    # Match everything between min0_ (or similar prefix) and _seed_
    m = re.search(r'_min\d+_(.+?)_seed_', base)
    if m:
        return m.group(1)
    # Fallback: try to find genus_species pattern before _seed_
    m = re.search(r'_([A-Z][a-z]+_[a-z]+(?:_[a-z]+)*)_seed_', base)
    if m:
        return m.group(1)
    return None


def parse_metabolite_from_rowname(rowname: str, bacteria: str) -> str:
    """
    Extract metabolite name from Row.names column.

    Row.names format: "Isobutyric acid-Akkermansia_muciniphila"
    """
    # Remove quotes
    rowname = rowname.strip().strip('"')
    # Split by the bacteria suffix
    suffix = f"-{bacteria}"
    if suffix in rowname:
        return rowname.split(suffix)[0]
    # Fallback: split by last hyphen followed by uppercase
    parts = rowname.rsplit('-', 1)
    if len(parts) == 2:
        return parts[0]
    return rowname


def process_single_file(
    filepath: str,
    top_n: int = 200,
) -> list:
    """
    Process a single coordinates CSV file and return top-N candidates.

    Args:
        filepath: Path to a single *_coordinates_tunek.csv file
        top_n: Number of top candidates to return

    Returns:
        List of candidate dicts sorted by mmsage_norm descending
    """
    fpath = Path(filepath)
    if not fpath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    fname = fpath.name
    source_file = portable_source_path(fpath)
    bacteria = extract_bacteria_from_filename(fname)
    if not bacteria:
        raise ValueError(f"Cannot extract bacteria name from filename: {fname}")

    print(f"Processing single file: {fname}")
    print(f"Bacteria: {bacteria}")

    df = pd.read_csv(filepath)

    if 'Pseudotime' not in df.columns:
        raise ValueError(f"'Pseudotime' column not found in {fname}")

    # Find the Row.names column
    rowname_col = None
    if 'Row.names' in df.columns:
        rowname_col = 'Row.names'
    elif df.columns[0] == '' or df.columns[0].startswith('Unnamed'):
        if len(df.columns) > 1 and 'Row.names' in df.columns.tolist():
            rowname_col = 'Row.names'
        else:
            rowname_col = df.columns[1] if len(df.columns) > 1 else df.columns[0]
    else:
        rowname_col = df.columns[0]

    max_pt = df['Pseudotime'].max()
    if max_pt <= 0 or pd.isna(max_pt):
        raise ValueError(f"Invalid max pseudotime: {max_pt}")

    all_candidates = []
    for _, row in df.iterrows():
        rowname = str(row.get(rowname_col, ''))
        metabolite = parse_metabolite_from_rowname(rowname, bacteria)
        pseudotime = row['Pseudotime']

        if pd.isna(pseudotime):
            continue

        all_candidates.append({
            'bacteria': bacteria,
            'metabolite': metabolite,
            'pseudotime': float(pseudotime),
            'max_pseudotime': float(max_pt),
            'source_file': source_file,
        })

    print(f"Total raw pairs: {len(all_candidates)}")

    if not all_candidates:
        return []

    df_all = pd.DataFrame(all_candidates)
    # Rank-based dual-end signal: distance to nearest endpoint (root or tip),
    # then rank so that closest-to-endpoint = highest score.
    # This avoids the linear normalization problem where most values cluster near 1.0.
    def _rank_dual_end(group):
        pt = group['pseudotime']
        pt_max = pt.max()
        dist_to_end = pt.apply(lambda x: min(x, pt_max - x))
        n = len(group)
        rank = dist_to_end.rank(ascending=True, method='first')
        group['mmsage_norm'] = 1.0 - (rank - 1) / n
        group['rank_in_microbe'] = pt.rank(ascending=False, method='first').astype(int)
        group['total_in_microbe'] = n
        return group

    # Apply directly for one file to avoid Pandas 2.0+ groupby index loss.
    df_all = _rank_dual_end(df_all)
    if 'bacteria' not in df_all.columns:
        df_all = df_all.reset_index()

    df_all = df_all.sort_values('mmsage_norm', ascending=False)
    df_dedup = df_all.drop_duplicates(subset=['bacteria', 'metabolite'], keep='first')
    df_top = df_dedup.head(top_n)

    results = []
    for _, row in df_top.iterrows():
        results.append({
            'bacteria': row['bacteria'],
            'metabolite': row['metabolite'],
            'pseudotime': round(row['pseudotime'], 4),
            'mmsage_norm': round(row['mmsage_norm'], 4),
            'rank_in_microbe': int(row['rank_in_microbe']),
            'total_in_microbe': int(row['total_in_microbe']),
            'source_file': row['source_file'],
        })

    print(f"Returning top {len(results)} candidates")
    return results


def process_coordinates_dir(
    coordinates_dir: str,
    top_n: int = 200,
    bacteria_filter: Optional[str] = None,
) -> list:
    """
    Process all coordinates CSV files and return top-N candidates.

    Args:
        coordinates_dir: Path to directory containing *_coordinates_tunek.csv files
        top_n: Number of top candidates to return
        bacteria_filter: If set, only process files containing this string

    Returns:
        List of candidate dicts sorted by mmsage_norm descending
    """
    coord_dir = Path(coordinates_dir)
    pattern = str(coord_dir / "*_coordinates_tunek.csv")
    files = glob.glob(pattern)

    if not files:
        raise FileNotFoundError(
            f"No *_coordinates_tunek.csv files found in {coordinates_dir}"
        )

    print(f"Found {len(files)} coordinates files")

    if bacteria_filter:
        files = [f for f in files if bacteria_filter in f]
        print(f"Filtered to {len(files)} files containing '{bacteria_filter}'")

    all_candidates = []

    for i, filepath in enumerate(files):
        fname = Path(filepath).name
        source_file = portable_source_path(filepath)
        bacteria = extract_bacteria_from_filename(fname)
        if not bacteria:
            continue

        try:
            df = pd.read_csv(filepath)
        except Exception:
            continue

        # Ensure required columns exist
        if 'Pseudotime' not in df.columns:
            continue

        # Find the Row.names column (could be unnamed first col or 'Row.names')
        rowname_col = None
        if 'Row.names' in df.columns:
            rowname_col = 'Row.names'
        elif df.columns[0] == '' or df.columns[0].startswith('Unnamed'):
            # Check second column
            if len(df.columns) > 1 and 'Row.names' in df.columns.tolist():
                rowname_col = 'Row.names'
            else:
                rowname_col = df.columns[1] if len(df.columns) > 1 else df.columns[0]
        else:
            rowname_col = df.columns[0]

        # Parse each row
        max_pt = df['Pseudotime'].max()
        if max_pt <= 0 or pd.isna(max_pt):
            continue

        for _, row in df.iterrows():
            rowname = str(row.get(rowname_col, ''))
            metabolite = parse_metabolite_from_rowname(rowname, bacteria)
            pseudotime = row['Pseudotime']

            if pd.isna(pseudotime):
                continue

            all_candidates.append({
                'bacteria': bacteria,
                'metabolite': metabolite,
                'pseudotime': float(pseudotime),
                'max_pseudotime': float(max_pt),
                'source_file': source_file,
            })

        if (i + 1) % 500 == 0:
            print(f"  Processed {i + 1}/{len(files)} files, "
                  f"{len(all_candidates)} pairs so far")

    print(f"Total raw pairs: {len(all_candidates)}")

    if not all_candidates:
        return []

    # Build DataFrame for efficient processing
    df_all = pd.DataFrame(all_candidates)

    # Rank-based dual-end signal (same logic as process_single_file)
    def _rank_dual_end(group):
        pt = group['pseudotime']
        pt_max = pt.max()
        dist_to_end = pt.apply(lambda x: min(x, pt_max - x))
        n = len(group)
        rank = dist_to_end.rank(ascending=True, method='first')
        group['mmsage_norm'] = 1.0 - (rank - 1) / n
        group['rank_in_microbe'] = pt.rank(ascending=False, method='first').astype(int)
        group['total_in_microbe'] = n
        return group

    df_all = df_all.groupby('bacteria', group_keys=False).apply(_rank_dual_end)
    # Restore bacteria from the index after grouped processing if needed.
    if 'bacteria' not in df_all.columns:
        df_all = df_all.reset_index()

    # Sort globally by mmsage_norm descending, take top_n
    df_all = df_all.sort_values('mmsage_norm', ascending=False)

    # Deduplicate: keep best score per (bacteria, metabolite) pair
    df_dedup = df_all.drop_duplicates(subset=['bacteria', 'metabolite'], keep='first')
    df_top = df_dedup.head(top_n)

    # Build output
    results = []
    for _, row in df_top.iterrows():
        results.append({
            'bacteria': row['bacteria'],
            'metabolite': row['metabolite'],
            'pseudotime': round(row['pseudotime'], 4),
            'mmsage_norm': round(row['mmsage_norm'], 4),
            'rank_in_microbe': int(row['rank_in_microbe']),
            'total_in_microbe': int(row['total_in_microbe']),
            'source_file': row['source_file'],
        })

    print(f"Returning top {len(results)} candidates")
    return results


def run_step1(
    coordinates_dir: str,
    output_path: str,
    top_n: int = 200,
    bacteria_filter: Optional[str] = None,
    coordinates_file: Optional[str] = None,
) -> list:
    """
    Run step 1 and save results.

    Args:
        coordinates_dir: Path to coordinates CSV directory
        output_path: Path to save output JSON
        top_n: Number of top candidates
        bacteria_filter: Optional filter for bacteria name in filename
        coordinates_file: If set, process this single file instead of directory

    Returns:
        List of candidate dicts
    """
    print("=" * 60)
    print("  STEP 1: MMSage Signal Processing")
    print("=" * 60)

    if coordinates_file:
        candidates = process_single_file(coordinates_file, top_n)
    else:
        candidates = process_coordinates_dir(coordinates_dir, top_n, bacteria_filter)

    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(candidates, f, indent=2, ensure_ascii=False)

    print(f"Saved {len(candidates)} candidates to {output_path}")
    return candidates


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Step 1: MMSage Signal Processing')
    _project_dir = str(Path(__file__).parent.parent)
    parser.add_argument('--coordinates-dir',
                        default=_project_dir + '/data/sample_coordinates',
                        help='Path to coordinates CSV directory')
    parser.add_argument('--output', default='outputs/step1_candidates.json')
    parser.add_argument('--top-n', type=int, default=200)
    parser.add_argument('--bacteria-filter', default=None,
                        help='Only process files containing this string')
    parser.add_argument('--coordinates-file',
                        default=_project_dir + '/data/sample_coordinates/pluscombno1V0317_min0_Akkermansia_muciniphila_seed_1_dim_10_neighbor_2_dist_0.4_metric_euclidean_rank_rootknow_cor0303_1_top_50_Pthre_0.1_pair.csv_clu_1_coordinates_tunek.csv',
                        help='Process a single coordinates CSV file instead of directory')
    args = parser.parse_args()

    run_step1(args.coordinates_dir, args.output, args.top_n, args.bacteria_filter,
              coordinates_file=args.coordinates_file)
