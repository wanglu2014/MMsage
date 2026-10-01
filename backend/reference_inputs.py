"""Stage uploaded tables in the layout consumed by the M1–M4 scripts."""
import csv
import math
from pathlib import Path
from scipy.stats import spearmanr


def read_table(path):
    path = Path(path)
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.reader(stream, delimiter="," if path.suffix.lower() == ".csv" else "\t"))
    if len(rows) < 2 or len(rows[0]) < 2:
        raise ValueError(f"{path.name}: requires a header and data rows")
    if any(len(row) != len(rows[0]) for row in rows):
        raise ValueError(f"{path.name}: inconsistent column counts")
    return rows


def matrix(path, allow_na=False):
    rows = read_table(path)
    for ids in (rows[0][1:], [row[0] for row in rows[1:]]):
        if any(not value.strip() for value in ids) or len(set(ids)) != len(ids):
            raise ValueError(f"{Path(path).name}: missing or duplicate IDs")
    for row in rows[1:]:
        for value in row[1:]:
            if allow_na and value in ("NA", ""):
                continue
            if not math.isfinite(float(value)):
                raise ValueError(f"{Path(path).name}: all values must be finite")
    return rows


def write_table(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        csv.writer(stream, delimiter="\t", lineterminator="\n").writerows(rows)


def align_common_samples(mic, met, minimum=3):
    """Keep shared samples using the metabolite-file order used by TreeMM."""
    mic_positions = {sample: i for i, sample in enumerate(mic[0][1:], start=1)}
    shared = [sample for sample in met[0][1:] if sample in mic_positions]
    if len(shared) < minimum:
        raise ValueError(f"At least {minimum} shared samples are required")
    mic_cols = [mic_positions[sample] for sample in shared]
    met_positions = {sample: i for i, sample in enumerate(met[0][1:], start=1)}
    met_cols = [met_positions[sample] for sample in shared]
    aligned_mic = [[mic[0][0]] + shared]
    aligned_met = [[met[0][0]] + shared]
    aligned_mic.extend([[row[0]] + [row[i] for i in mic_cols] for row in mic[1:]])
    aligned_met.extend([[row[0]] + [row[i] for i in met_cols] for row in met[1:]])
    return aligned_mic, aligned_met


def derive_reference_tables(microbe_path, metabolite_path, job_dir, taxonomy_source=None):
    """Create the M1/M2 auxiliary tables from the two uploaded matrices.

    Taxonomy is inferred from feature IDs. This is useful when IDs already contain
    species names; CAG-only IDs remain valid inputs but will not match a named M1
    target unless the ID itself contains that target name.
    """
    mic, met = align_common_samples(matrix(microbe_path), matrix(metabolite_path))
    samples = mic[0][1:]
    mic_ids = [row[0] for row in mic[1:]]
    met_ids = [row[0] for row in met[1:]]
    taxonomy = [["CAG", "species"]]
    source_map = {}
    if taxonomy_source and Path(taxonomy_source).exists():
        source_rows = read_table(taxonomy_source)
        headers = {name.strip().lower(): i for i, name in enumerate(source_rows[0])}
        id_col = headers.get("id", 0)
        species_col = headers.get("speciescluster")
        if species_col is not None:
            for row in source_rows[1:]:
                if len(row) > max(id_col, species_col) and row[species_col] not in ("", "NA"):
                    source_map[row[id_col]] = row[species_col]
    taxonomy.extend([[feature, source_map.get(feature, feature)] for feature in mic_ids])
    corr = [["CAG"] + met_ids]
    pvals = [["CAG"] + met_ids]
    for mrow in mic[1:]:
        x = [float(v) for v in mrow[1:]]
        crow, prow = [mrow[0]], [mrow[0]]
        for met_row in met[1:]:
            met_id = met_row[0]
            y = [float(v) for v in met_row[1:]]
            rho, p = spearmanr(x, y, nan_policy="omit")
            crow.append("NA" if not math.isfinite(float(rho)) else str(float(rho)))
            prow.append("NA" if not math.isfinite(float(p)) else str(float(p)))
        corr.append(crow)
        pvals.append(prow)
    job_dir = Path(job_dir)
    write_table(job_dir / "data/taxon_names.tsv", taxonomy)
    write_table(job_dir / "results/method02_spearman_matrix.tsv", corr)
    write_table(job_dir / "results/method02_spearman_pvalues.tsv", pvals)
    return taxonomy, corr, pvals


def stage_inputs(job_dir, microbe_path, metabolite_path, taxonomy_path=None, correlation_path=None, pvalues_path=None, taxonomy_source=None):
    job_dir = Path(job_dir)
    mic, met = align_common_samples(matrix(microbe_path), matrix(metabolite_path))
    if any(float(v) < 0 for row in mic[1:] for v in row[1:]) or not any(float(v) > 0 for row in mic[1:] for v in row[1:]):
        raise ValueError("Microbial abundance must be non-negative with at least one positive value")
    if taxonomy_path is None or correlation_path is None or pvalues_path is None:
        derive_reference_tables(microbe_path, metabolite_path, job_dir, taxonomy_source=taxonomy_source)
        taxonomy_path = job_dir / "data/taxon_names.tsv"
        correlation_path = job_dir / "results/method02_spearman_matrix.tsv"
        pvalues_path = job_dir / "results/method02_spearman_pvalues.tsv"
    tax = read_table(taxonomy_path)
    if len(tax[0]) != 2:
        raise ValueError("Taxonomy must contain exactly two columns: CAG and species")
    # Preserve identifiers used in output filenames; reject path components.
    for row in mic[1:]:
        if any(c in row[0] for c in '/\\:*?"<>|') or row[0] in ('.', '..'):
            raise ValueError("Microbial IDs must be valid filename components")
    corr, pvals = matrix(correlation_path, allow_na=True), matrix(pvalues_path, allow_na=True)
    if corr[0][1:] != pvals[0][1:] or [r[0] for r in corr[1:]] != [r[0] for r in pvals[1:]]:
        raise ValueError("Correlation and P-value matrices must have identical row and column IDs in identical order")
    if not set(corr[0][1:]).issubset({r[0] for r in met[1:]}):
        raise ValueError("Correlation metabolite IDs must exist in the abundance matrix")
    if any(not -1 <= float(v) <= 1 for r in corr[1:] for v in r[1:] if v not in ("NA", "")):
        raise ValueError("Spearman coefficients must lie between -1 and 1")
    if any(not 0 <= float(v) <= 1 for r in pvals[1:] for v in r[1:] if v not in ("NA", "")):
        raise ValueError("P-values must lie between 0 and 1")
    for name, rows in (("data/microbes_wide.tsv", mic), ("data/metabolites_wide.tsv", met),
                       ("data/taxon_names.tsv", tax), ("results/method02_spearman_matrix.tsv", corr),
                       ("results/method02_spearman_pvalues.tsv", pvals)):
        write_table(job_dir / name, rows)
