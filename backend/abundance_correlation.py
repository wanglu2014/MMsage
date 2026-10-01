"""Independent microbe--metabolite correlation ranking."""
from pathlib import Path
import pandas as pd
from scipy.stats import spearmanr


def _read_matrix(path: str | Path, label: str) -> pd.DataFrame:
    path = Path(path)
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    frame = pd.read_csv(path, sep=sep)
    if frame.shape[1] < 3:
        raise ValueError(f"{label} must contain an ID column and at least two samples")
    frame = frame.rename(columns={frame.columns[0]: "feature_id"})
    if frame["feature_id"].isna().any() or frame["feature_id"].astype(str).duplicated().any():
        raise ValueError(f"{label} contains missing or duplicate feature IDs")
    frame["feature_id"] = frame["feature_id"].astype(str)
    frame = frame.set_index("feature_id")
    try:
        frame = frame.apply(pd.to_numeric, errors="raise")
    except Exception as exc:
        raise ValueError(f"{label} contains non-numeric abundance values") from exc
    return frame


def rank_correlations(microbe_path: str | Path, metabolite_path: str | Path) -> pd.DataFrame:
    microbes = _read_matrix(microbe_path, "Microbe matrix")
    metabolites = _read_matrix(metabolite_path, "Metabolite matrix")
    samples = [s for s in microbes.columns if s in metabolites.columns]
    if len(samples) < 3:
        raise ValueError("At least three shared sample columns are required")
    microbes, metabolites = microbes[samples], metabolites[samples]
    rows = []
    for microbe, x in microbes.iterrows():
        if x.nunique(dropna=True) < 2:
            continue
        for metabolite, y in metabolites.iterrows():
            valid = x.notna() & y.notna()
            n = int(valid.sum())
            if n < 3 or x[valid].nunique() < 2 or y[valid].nunique() < 2:
                continue
            rho, p_value = spearmanr(x[valid].to_numpy(), y[valid].to_numpy())
            if pd.notna(rho) and pd.notna(p_value):
                rows.append({"microbe_id": microbe, "metabolite_id": metabolite,
                             "rho": float(rho), "p_value": float(p_value), "n_samples": n})
    result = pd.DataFrame(rows)
    if result.empty:
        return pd.DataFrame(columns=["rank", "microbe_id", "metabolite_id", "rho", "p_value", "fdr", "direction", "n_samples"])
    result = result.sort_values(["p_value", "rho"], ascending=[True, False]).reset_index(drop=True)
    # Benjamini-Hochberg correction without an additional dependency.
    m = len(result)
    ranks = pd.Series(range(1, m + 1), index=result.index)
    result["fdr"] = (result["p_value"] * m / ranks).iloc[::-1].cummin().iloc[::-1].clip(upper=1.0)
    result = result.sort_values("rho", key=lambda s: s.abs(), ascending=False).reset_index(drop=True)
    result["rank"] = range(1, len(result) + 1)
    result["direction"] = result["rho"].map(lambda value: "Positive" if value > 0 else "Negative")
    return result[["rank", "microbe_id", "metabolite_id", "rho", "p_value", "fdr", "direction", "n_samples"]]
