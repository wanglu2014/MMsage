from pathlib import Path
import hashlib
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from scipy import stats

ROOT_DIR = Path(__file__).resolve().parents[1]
PAPER_ROOT = ROOT_DIR.parent
DATA_ROOT = ROOT_DIR.parent / "data"
SRC = DATA_ROOT / "soil_0512" / "soilori-ROC_0412result_replaced.csv"
SIM_RELEASE_SRC = DATA_ROOT / "sim_0512" / "Sim release ROC_0206resmodi.csv"
SIM_CONSUMPTION_SRC = DATA_ROOT / "sim_0512" / "Sim cons_ROC_0303_resneg.csv"
TARGET_DIR = ROOT_DIR / "output" / "input_cache"
TARGET_FILE = TARGET_DIR / "soil_ROC_0412result_replaced.csv"
OUTPUT_DIR = ROOT_DIR / "output"
CSV_DIR = OUTPUT_DIR / "csv"
FIG_DIR = OUTPUT_DIR / "figures"


def sha256sum(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def trace_path(path: Path) -> str:
    """Return a package-relative path suitable for archived provenance records."""
    return path.resolve().relative_to(PAPER_ROOT.resolve()).as_posix()


def trapz_metric(x: np.ndarray, y: np.ndarray) -> float:
    d = pd.DataFrame({"x": x, "y": y}).dropna()
    d = d.sort_values("x")
    d = d.groupby("x", as_index=False)["y"].mean()
    if len(d) < 2:
        return float("nan")
    return float(np.trapezoid(d["y"].to_numpy(), d["x"].to_numpy()))


def holm_adjust(p_values: np.ndarray) -> np.ndarray:
    p_values = np.asarray(p_values, dtype=float)
    m = len(p_values)
    order = np.argsort(p_values)
    adjusted = np.empty(m, dtype=float)
    for rank, idx in enumerate(order):
        adjusted[idx] = (m - rank) * p_values[idx]
    for i in range(1, m):
        adjusted[order[i]] = max(adjusted[order[i]], adjusted[order[i - 1]])
    return np.clip(adjusted, 0.0, 1.0)


def rank_biserial_from_paired(diff: np.ndarray) -> float:
    nonzero = diff[diff != 0]
    n = len(nonzero)
    if n == 0:
        return 0.0
    ranks = stats.rankdata(np.abs(nonzero), method="average")
    w_pos = float(np.sum(ranks[nonzero > 0]))
    w_neg = float(np.sum(ranks[nonzero < 0]))
    denom = w_pos + w_neg
    if denom == 0:
        return 0.0
    return (w_pos - w_neg) / denom


def save_plot(fig, path: Path):
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def analyze_metric_vs_best(df: pd.DataFrame, metric_col: str, best_method: str, higher_is_better: bool = True) -> pd.DataFrame:
    # Explicitly pair by X and then run paired Wilcoxon tests.
    rows = []

    anchor_df = (
        df.loc[df["Files"] == best_method, ["X", metric_col]]
        .rename(columns={metric_col: "anchor_value"})
        .copy()
    )

    compared_methods = sorted([m for m in df["Files"].unique().tolist() if m != best_method])
    for m in compared_methods:
        m_df = (
            df.loc[df["Files"] == m, ["X", metric_col]]
            .rename(columns={metric_col: "method_value"})
            .copy()
        )

        paired = pd.merge(anchor_df, m_df, on="X", how="inner").sort_values("X")
        anchor_vals = paired["anchor_value"].to_numpy()
        method_vals = paired["method_value"].to_numpy()

        if higher_is_better:
            oriented_diff = anchor_vals - method_vals
            x_for_test = anchor_vals
            y_for_test = method_vals
        else:
            oriented_diff = method_vals - anchor_vals
            x_for_test = method_vals
            y_for_test = anchor_vals

        nonzero = oriented_diff[oriented_diff != 0]
        if len(nonzero) == 0:
            stat_val = 0.0
            p_val = 1.0
        else:
            res = stats.wilcoxon(x_for_test, y_for_test, alternative="greater", zero_method="wilcox", correction=False, mode="auto")
            stat_val = float(res.statistic)
            p_val = float(res.pvalue)

        rows.append(
            {
                "metric": metric_col,
                "direction": "higher_better" if higher_is_better else "lower_better",
                "best_method": best_method,
                "compared_method": m,
                "paired_by": "X",
                "paired_n": int(len(paired)),
                "paired_x_min": float(paired["X"].min()) if len(paired) > 0 else np.nan,
                "paired_x_max": float(paired["X"].max()) if len(paired) > 0 else np.nan,
                "test_name": "wilcoxon_signed_rank",
                "alternative": "greater",
                "zero_method": "wilcox",
                "wilcoxon_statistic": stat_val,
                "p_value_raw": p_val,
                "hl_median_diff_oriented": float(np.median(oriented_diff)) if len(oriented_diff) > 0 else np.nan,
                "rank_biserial": float(rank_biserial_from_paired(oriented_diff)) if len(oriented_diff) > 0 else np.nan,
                "best_better_rate": float(np.mean(oriented_diff > 0)) if len(oriented_diff) > 0 else np.nan,
                "best_equal_rate": float(np.mean(oriented_diff == 0)) if len(oriented_diff) > 0 else np.nan,
                "best_worse_rate": float(np.mean(oriented_diff < 0)) if len(oriented_diff) > 0 else np.nan,
                "n_points": int(len(oriented_diff)),
                "n_nonzero_diff": int(len(nonzero)),
            }
        )

    out = pd.DataFrame(rows)
    if len(out) > 0:
        out["p_value_holm"] = holm_adjust(out["p_value_raw"].to_numpy())
        out["significant_best_better_alpha_0_05"] = out["p_value_holm"] < 0.05
        out = out.sort_values(["p_value_holm", "hl_median_diff_oriented"], ascending=[True, False]).reset_index(drop=True)
    return out


def format_p_value(p: float) -> str:
    if pd.isna(p):
        return "p=NA"
    if p < 1e-4:
        return f"p={p:.1e}"
    return f"p={p:.4f}"


def build_metric_panel_specs_8():
    return [
        {"metric_col": "TP", "metric_label": "TP", "higher_is_better": True, "panel": "A"},
        {"metric_col": "recall", "metric_label": "Recall", "higher_is_better": True, "panel": "B"},
        {"metric_col": "precision", "metric_label": "Precision", "higher_is_better": True, "panel": "C"},
        {"metric_col": "f1_score", "metric_label": "F1-Score", "higher_is_better": True, "panel": "D"},
        {"metric_col": "FP", "metric_label": "FP", "higher_is_better": False, "panel": "E"},
        {"metric_col": "FN", "metric_label": "FN", "higher_is_better": False, "panel": "F"},
        {"metric_col": "TPR", "metric_label": "TPR", "higher_is_better": True, "panel": "G"},
        {"metric_col": "TNR", "metric_label": "TNR", "higher_is_better": True, "panel": "H"},
    ]


def normalize_method_name(method_name: str) -> str:
    name = str(method_name).strip()
    low = name.lower()
    if low == "back":
        return "MMSage"
    if low == "noback":
        return "MMSage-noback"
    return name


def classify_method_group(method_name: str) -> str:
    name = str(method_name).strip().lower()
    if name in {"mmsage", "back"}:
        return "mmsage"
    if name in {"mmsage-noback", "noback"}:
        return "mmsage-noback"

    if any(k in name for k in ["mmvec", "spiec-easi", "spieceasi", "sparcc", "mminp", "mimosa"]):
        return "network"
    if any(k in name for k in ["pearson", "spearman", "phi", "rho"]):
        return "traditional"
    return "other"


def build_method_color_map(metric_values: pd.DataFrame, ordered_methods: list, higher_is_better: bool, reference_method: str) -> dict:
    if not reference_method:
        return {m: "#D0D0D0" for m in ordered_methods}

    ref_vals = metric_values.loc[metric_values["Files"] == reference_method, "value"].dropna()
    if len(ref_vals) == 0:
        return {m: "#D0D0D0" for m in ordered_methods}

    ref_med = float(ref_vals.median())
    color_map = {}
    for m in ordered_methods:
        vals = metric_values.loc[metric_values["Files"] == m, "value"].dropna()
        med = float(vals.median()) if len(vals) > 0 else np.nan
        if pd.isna(med):
            color_map[m] = "#D0D0D0"
        else:
            better = med > ref_med if higher_is_better else med < ref_med
            color_map[m] = "#4C78A8" if better else "#D0D0D0"
    return color_map


def sort_methods_for_plot(methods: list) -> list:
    methods = [normalize_method_name(m) for m in methods]
    unique = []
    seen = set()
    for m in methods:
        if m not in seen:
            unique.append(m)
            seen.add(m)
    front = [m for m in unique if m in {"MMSage", "MMSage-noback"}]
    rest = sorted([m for m in unique if m not in {"MMSage", "MMSage-noback"}], key=lambda x: x.lower())
    ordered = []
    for m in ["MMSage", "MMSage-noback"]:
        if m in front:
            ordered.append(m)
    ordered.extend(rest)
    return ordered


def p_to_star(p: float) -> str:
    if pd.isna(p):
        return ""
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


def add_significance_bracket(ax, x1: int, x2: int, y: float, height: float, text: str):
    ax.plot([x1, x1, x2, x2], [y, y + height, y + height, y], lw=0.9, c="black", clip_on=False)
    ax.text((x1 + x2) / 2, y + height, text, ha="center", va="bottom", fontsize=8, color="black", fontweight="bold")


def resolve_anchor_method(ordered_methods: list, target: str):
    target = str(target).strip().lower()
    alias = {"mmsage": "back", "mmsage-noback": "noback"}.get(target, target)
    for m in ordered_methods:
        if str(m).strip().lower() == target:
            return m
    for m in ordered_methods:
        if str(m).strip().lower() == alias:
            return m
    return None


def build_star_map(anchor_df: pd.DataFrame, ordered_methods: list) -> dict:
    star_map = {m: "" for m in ordered_methods}
    if anchor_df is None or len(anchor_df) == 0:
        return star_map
    for _, r in anchor_df.iterrows():
        method = normalize_method_name(r.get("compared_method"))
        if method in star_map:
            star_map[method] = p_to_star(float(r.get("p_value_holm", np.nan)))
    return star_map


def draw_single_metric_panel(
    ax,
    metric_values: pd.DataFrame,
    ordered_methods: list,
    metric_label: str,
    panel_letter: str,
    color_map: dict,
    noback_star_map: dict,
    back_star_map: dict,
    noback_name: str,
    back_name: str,
    mmsage_pair_p: float,
):
    palette = [mcolors.to_rgba(color_map.get(m, "#BDBDBD"), alpha=0.45) for m in ordered_methods]

    sns.boxplot(
        data=metric_values,
        x="Files",
        y="value",
        hue="Files",
        hue_order=ordered_methods,
        order=ordered_methods,
        dodge=False,
        legend=False,
        ax=ax,
        showfliers=False,
        palette=palette,
        width=0.38,
        linewidth=0.8,
        boxprops={"edgecolor": "black", "linewidth": 0.8},
        whiskerprops={"color": "black", "linewidth": 0.8},
        capprops={"color": "black", "linewidth": 0.8},
        medianprops={"color": "black", "linewidth": 0.9},
    )

    sns.stripplot(
        data=metric_values,
        x="Files",
        y="value",
        order=ordered_methods,
        ax=ax,
        color="gray",
        size=1.8,
        jitter=0.18,
        alpha=0.28,
        linewidth=0,
    )

    values_all = metric_values["value"].dropna()
    if len(values_all) > 0:
        ymin = float(values_all.min())
        ymax = float(values_all.max())
        yr = ymax - ymin if ymax > ymin else max(abs(ymax), 1.0)
    else:
        ymin, ymax, yr = 0.0, 1.0, 1.0

    if noback_name in ordered_methods:
        noback_vals = metric_values.loc[metric_values["Files"] == noback_name, "value"].dropna()
        if len(noback_vals) > 0:
            ax.axhline(float(noback_vals.median()), ls="--", lw=1, alpha=0.7, color="black")

    for idx, m in enumerate(ordered_methods):
        if m == noback_name:
            continue
        m_vals = metric_values.loc[metric_values["Files"] == m, "value"].dropna()
        if len(m_vals) == 0:
            continue
        base_y = float(m_vals.max()) + 0.04 * yr

        star_noback = noback_star_map.get(m, "")
        if star_noback:
            ax.text(idx, base_y, star_noback, ha="center", va="bottom", fontsize=9, color="black", fontweight="bold")

        star_back = back_star_map.get(m, "")
        if star_back:
            plus = "+" * len(star_back)
            ax.text(idx, base_y + 0.085 * yr, plus, ha="center", va="bottom", fontsize=8, color="#444444")

    top_pad = 0.22 * yr
    if (
        back_name
        and noback_name
        and (back_name in ordered_methods)
        and (noback_name in ordered_methods)
        and (not pd.isna(mmsage_pair_p))
    ):
        p_text = p_to_star(float(mmsage_pair_p))
        if p_text == "":
            p_text = format_p_value(float(mmsage_pair_p))
        x1 = ordered_methods.index(back_name)
        x2 = ordered_methods.index(noback_name)
        bracket_y = ymax + 0.06 * yr
        add_significance_bracket(ax, min(x1, x2), max(x1, x2), bracket_y, 0.035 * yr, p_text)
        top_pad = 0.32 * yr

    ax.set_ylim(ymin - 0.02 * yr, ymax + top_pad)
    ax.set_title(metric_label, fontsize=7)
    ax.set_xlabel("")
    ax.set_ylabel(metric_label, fontsize=7)
    ax.tick_params(axis="x", rotation=60, labelsize=5.5)
    ax.tick_params(axis="y", labelsize=7)
    for label in ax.get_xticklabels():
        label.set_ha("right")
        label.set_rotation_mode("anchor")
    ax.text(0.01, 0.98, panel_letter, transform=ax.transAxes, ha="left", va="top", fontsize=7, fontweight="bold")
    ax.set_facecolor("white")
    sns.despine(ax=ax, top=True, right=True)


def plot_nejm_3x4_metrics_boxplot(
    df: pd.DataFrame,
    ordered_methods: list,
    existing_anchors: list,
    csv_dir: Path,
    fig_dir: Path,
    out_name: str,
):
    panel_specs = build_metric_panel_specs_8()

    all_metric_test_rows = []
    metric_plot_summary_rows = []
    metric_boxplot_data_rows = []
    pairwise_excel_rows = []

    noback_name = resolve_anchor_method(ordered_methods, "MMSage-noback")
    back_name = resolve_anchor_method(ordered_methods, "MMSage")

    fig, axes = plt.subplots(4, 2, figsize=(10.8, 8.27), facecolor="white")
    axes = axes.flatten()

    for ax, spec in zip(axes, panel_specs):
        metric_col = spec["metric_col"]
        metric_label = spec["metric_label"]
        higher_is_better = spec["higher_is_better"]

        metric_values = df[["Files", "X", metric_col]].copy()
        metric_values = metric_values.rename(columns={metric_col: "value"})
        metric_values["metric"] = metric_col
        metric_values["metric_label"] = metric_label
        metric_values.to_csv(csv_dir / f"boxplot_data_{metric_col}.csv", index=False)
        metric_boxplot_data_rows.append(metric_values)

        mat = df.pivot(index="X", columns="Files", values=metric_col).sort_index()
        var_per_method = mat.var(axis=0, ddof=1)
        if np.all(np.nan_to_num(var_per_method.to_numpy(), nan=0.0) < 1e-15):
            fried_stat, fried_p = np.nan, np.nan
        else:
            fried_stat, fried_p = stats.friedmanchisquare(*[mat[c].to_numpy() for c in mat.columns])

        anchor_results = {}
        for anchor in existing_anchors:
            anchor_df = analyze_metric_vs_best(df, metric_col, anchor, higher_is_better=higher_is_better)
            anchor_df.to_csv(csv_dir / f"wilcoxon_vs_{anchor}_{metric_col}.csv", index=False)
            anchor_results[anchor] = anchor_df
            if len(anchor_df) > 0:
                tmp = anchor_df.copy()
                tmp["metric_label"] = metric_label
                tmp["scenario_metric"] = metric_col
                tmp["scenario"] = out_name
                all_metric_test_rows.append(tmp)
                pairwise_excel_rows.append(tmp)

        per_anchor_summary = {}
        for anchor in existing_anchors:
            adf = anchor_results.get(anchor, pd.DataFrame())
            if len(adf) == 0:
                per_anchor_summary[anchor] = {"sig": 0, "n": 0, "minp": np.nan}
            else:
                per_anchor_summary[anchor] = {
                    "sig": int(adf["significant_best_better_alpha_0_05"].sum()),
                    "n": int(len(adf)),
                    "minp": float(adf["p_value_holm"].min()),
                }

        mmsage_pair_p = np.nan
        if back_name and noback_name and back_name in mat.columns and noback_name in mat.columns:
            paired2 = mat[[back_name, noback_name]].dropna()
            d2 = paired2[back_name].to_numpy() - paired2[noback_name].to_numpy()
            nz2 = d2[d2 != 0]
            if len(nz2) == 0:
                mmsage_pair_p = 1.0
            else:
                r2 = stats.wilcoxon(
                    paired2[back_name].to_numpy(),
                    paired2[noback_name].to_numpy(),
                    alternative="two-sided",
                    zero_method="wilcox",
                    correction=False,
                    mode="auto",
                )
                mmsage_pair_p = float(r2.pvalue)

        metric_plot_summary_rows.append(
            {
                "scenario": out_name,
                "metric": metric_col,
                "metric_label": metric_label,
                "direction": "higher_better" if higher_is_better else "lower_better",
                "friedman_statistic": float(fried_stat) if not pd.isna(fried_stat) else np.nan,
                "friedman_p_value": float(fried_p) if not pd.isna(fried_p) else np.nan,
                "mmsage_noback_sig": per_anchor_summary.get("MMSage-noback", {}).get("sig", np.nan),
                "mmsage_noback_n": per_anchor_summary.get("MMSage-noback", {}).get("n", np.nan),
                "mmsage_noback_min_holm_p": per_anchor_summary.get("MMSage-noback", {}).get("minp", np.nan),
                "mmsage_sig": per_anchor_summary.get("MMSage", {}).get("sig", np.nan),
                "mmsage_n": per_anchor_summary.get("MMSage", {}).get("n", np.nan),
                "mmsage_min_holm_p": per_anchor_summary.get("MMSage", {}).get("minp", np.nan),
                "mmsage_vs_noback_p_raw": float(mmsage_pair_p) if not pd.isna(mmsage_pair_p) else np.nan,
                "boxplot_data_csv": f"boxplot_data_{metric_col}.csv",
                "mmsage_noback_csv": f"wilcoxon_vs_MMSage-noback_{metric_col}.csv" if "MMSage-noback" in existing_anchors else "",
                "mmsage_csv": f"wilcoxon_vs_MMSage_{metric_col}.csv" if "MMSage" in existing_anchors else "",
            }
        )

        noback_star_map = build_star_map(anchor_results.get(noback_name, pd.DataFrame()), ordered_methods) if noback_name else {m: "" for m in ordered_methods}
        back_star_map = build_star_map(anchor_results.get(back_name, pd.DataFrame()), ordered_methods) if back_name else {m: "" for m in ordered_methods}

        color_map = build_method_color_map(
            metric_values=metric_values,
            ordered_methods=ordered_methods,
            higher_is_better=higher_is_better,
            reference_method=noback_name,
        )

        draw_single_metric_panel(
            ax=ax,
            metric_values=metric_values,
            ordered_methods=ordered_methods,
            metric_label=metric_label,
            panel_letter=spec["panel"],
            color_map=color_map,
            noback_star_map=noback_star_map,
            back_star_map=back_star_map,
            noback_name=noback_name if noback_name else "",
            back_name=back_name if back_name else "",
            mmsage_pair_p=mmsage_pair_p,
        )

    for ax in axes[len(panel_specs):]:
        ax.axis("off")

    handles = [
        plt.Line2D([0], [0], color="#4C78A8", lw=8, label="Box > MMSage-noback"),
        plt.Line2D([0], [0], color="#D0D0D0", lw=8, label="Box <= MMSage-noback"),
        plt.Line2D([0], [0], color="#444444", marker="$+$", linestyle="None", markersize=10, label="+ = significant vs MMSage"),
        plt.Line2D([0], [0], color="#000000", marker="$*$", linestyle="None", markersize=10, label="* = significant vs MMSage-noback"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=8)
    fig.subplots_adjust(left=0.06, right=0.99, top=0.965, bottom=0.16, wspace=0.28, hspace=0.42)
    fig.savefig(fig_dir / out_name, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return all_metric_test_rows, metric_plot_summary_rows, metric_boxplot_data_rows, pairwise_excel_rows


def main():
    TARGET_DIR.mkdir(parents=True, exist_ok=True)
    CSV_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    if SRC.resolve() != TARGET_FILE.resolve():
        shutil.copy2(SRC, TARGET_FILE)

    src_hash = sha256sum(SRC)
    dst_hash = sha256sum(TARGET_FILE)
    copy_df = pd.DataFrame(
        [
            {
                "source": trace_path(SRC),
                "target": trace_path(TARGET_FILE),
                "source_size": SRC.stat().st_size,
                "target_size": TARGET_FILE.stat().st_size,
                "source_sha256": src_hash,
                "target_sha256": dst_hash,
                "size_match": SRC.stat().st_size == TARGET_FILE.stat().st_size,
                "hash_match": src_hash == dst_hash,
            }
        ]
    )
    copy_df.to_csv(CSV_DIR / "copy_verification.csv", index=False)

    df = pd.read_csv(TARGET_FILE)
    required_cols = {
        "Files",
        "X",
        "TP",
        "TN",
        "FP",
        "FN",
        "TPR",
        "TNR",
        "precision",
        "recall",
        "f1_score",
    }
    missing = sorted(list(required_cols - set(df.columns)))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = df.copy()
    df["Files"] = df["Files"].astype(str).map(normalize_method_name)
    df = df.loc[~df["Files"].isin(["phi Clr", "rho Clr"])].copy()
    df["X"] = pd.to_numeric(df["X"], errors="coerce")
    num_cols = ["TP", "TN", "FP", "FN", "TPR", "TNR", "precision", "recall", "f1_score"]
    for c in num_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    overview_df = pd.DataFrame(
        [
            {
                "n_rows": int(df.shape[0]),
                "n_cols": int(df.shape[1]),
                "n_methods": int(df["Files"].nunique()),
                "missing_values_total": int(df.isna().sum().sum()),
                "x_min": float(df["X"].min()),
                "x_max": float(df["X"].max()),
            }
        ]
    )
    overview_df.to_csv(CSV_DIR / "data_overview.csv", index=False)

    coverage = (
        df.groupby("Files")
        .agg(
            n_rows=("X", "size"),
            unique_x=("X", "nunique"),
            x_min=("X", "min"),
            x_max=("X", "max"),
        )
        .reset_index()
    )
    coverage["x_complete_1_100"] = (coverage["unique_x"] == 100) & (coverage["x_min"] == 1) & (coverage["x_max"] == 100)
    coverage.to_csv(CSV_DIR / "per_method_x_coverage.csv", index=False)

    eps = 1e-12
    calc_tpr = df["TP"] / (df["TP"] + df["FN"] + eps)
    calc_tnr = df["TN"] / (df["TN"] + df["FP"] + eps)
    calc_precision = df["TP"] / (df["TP"] + df["FP"] + eps)
    calc_recall = calc_tpr
    calc_f1 = 2 * calc_precision * calc_recall / (calc_precision + calc_recall + eps)

    consistency = pd.DataFrame(
        {
            "Files": df["Files"],
            "X": df["X"],
            "diff_TPR": (df["TPR"] - calc_tpr).abs(),
            "diff_TNR": (df["TNR"] - calc_tnr).abs(),
            "diff_precision": (df["precision"] - calc_precision).abs(),
            "diff_recall": (df["recall"] - calc_recall).abs(),
            "diff_f1_score": (df["f1_score"] - calc_f1).abs(),
        }
    )
    consistency.to_csv(CSV_DIR / "metric_recompute_diff_detail.csv", index=False)
    consistency_summary = pd.DataFrame(
        [
            {
                "max_diff_TPR": consistency["diff_TPR"].max(),
                "max_diff_TNR": consistency["diff_TNR"].max(),
                "max_diff_precision": consistency["diff_precision"].max(),
                "max_diff_recall": consistency["diff_recall"].max(),
                "max_diff_f1_score": consistency["diff_f1_score"].max(),
            }
        ]
    )
    consistency_summary.to_csv(CSV_DIR / "metric_recompute_diff_summary.csv", index=False)

    df["fpr"] = 1 - df["TNR"]
    df["balanced_accuracy"] = (df["TPR"] + df["TNR"]) / 2
    df["NPV"] = df["TN"] / (df["TN"] + df["FN"] + eps)

    mcc_num = df["TP"] * df["TN"] - df["FP"] * df["FN"]
    mcc_den = np.sqrt((df["TP"] + df["FP"]) * (df["TP"] + df["FN"]) * (df["TN"] + df["FP"]) * (df["TN"] + df["FN"]))
    df["mcc"] = np.where(mcc_den > 0, mcc_num / mcc_den, 0.0)

    summary_rows = []
    for method, g in df.groupby("Files"):
        g = g.sort_values("X")
        auroc = trapz_metric(g["fpr"].to_numpy(), g["TPR"].to_numpy())
        auprc = trapz_metric(g["recall"].to_numpy(), g["precision"].to_numpy())
        summary_rows.append(
            {
                "Files": method,
                "Median_F1": g["f1_score"].median(),
                "AUPRC": auprc,
                "AUROC": auroc,
                "Max_F1": g["f1_score"].max(),
                "Median_Balanced_Accuracy": g["balanced_accuracy"].median(),
                "Median_MCC": g["mcc"].median(),
                "Q10_F1": g["f1_score"].quantile(0.1),
                "IQR_F1": g["f1_score"].quantile(0.75) - g["f1_score"].quantile(0.25),
                "Segment1_Median_F1": g.loc[g["X"].between(1, 33), "f1_score"].median(),
                "Segment2_Median_F1": g.loc[g["X"].between(34, 66), "f1_score"].median(),
                "Segment3_Median_F1": g.loc[g["X"].between(67, 100), "f1_score"].median(),
            }
        )

    summary = pd.DataFrame(summary_rows)
    summary["Files"] = summary["Files"].map(normalize_method_name)
    summary = summary.groupby("Files", as_index=False).agg(
        Median_F1=("Median_F1", "max"),
        AUPRC=("AUPRC", "max"),
        AUROC=("AUROC", "max"),
        Max_F1=("Max_F1", "max"),
        Median_Balanced_Accuracy=("Median_Balanced_Accuracy", "max"),
        Median_MCC=("Median_MCC", "max"),
        Q10_F1=("Q10_F1", "max"),
        IQR_F1=("IQR_F1", "min"),
        Segment1_Median_F1=("Segment1_Median_F1", "max"),
        Segment2_Median_F1=("Segment2_Median_F1", "max"),
        Segment3_Median_F1=("Segment3_Median_F1", "max"),
    )
    summary = summary.sort_values(
        ["Median_F1", "AUPRC", "Q10_F1", "IQR_F1", "Files"],
        ascending=[False, False, False, True, True],
    ).reset_index(drop=True)
    summary["rank"] = np.arange(1, len(summary) + 1)
    summary.to_csv(CSV_DIR / "method_summary_metrics.csv", index=False)

    ordered_for_all = sort_methods_for_plot(summary["Files"].tolist())
    summary["Files"] = pd.Categorical(summary["Files"], categories=ordered_for_all, ordered=True)
    summary = summary.sort_values("Files").reset_index(drop=True)

    best_method = "MMSage-noback" if "MMSage-noback" in summary["Files"].astype(str).tolist() else summary.iloc[0]["Files"]

    f1_matrix = df.pivot(index="X", columns="Files", values="f1_score").sort_index()
    ba_matrix = df.pivot(index="X", columns="Files", values="balanced_accuracy").sort_index()
    f1_matrix = f1_matrix.reindex(columns=[m for m in ordered_for_all if m in f1_matrix.columns])
    ba_matrix = ba_matrix.reindex(columns=[m for m in ordered_for_all if m in ba_matrix.columns])

    f1_matrix.to_csv(CSV_DIR / "plot_f1_vs_x_data.csv")
    ba_matrix.to_csv(CSV_DIR / "plot_balanced_accuracy_vs_x_data.csv")

    roc_plot_data = df[["Files", "X", "fpr", "TPR"]].sort_values(["Files", "fpr", "TPR"])
    pr_plot_data = df[["Files", "X", "recall", "precision"]].sort_values(["Files", "recall", "precision"])
    roc_plot_data.to_csv(CSV_DIR / "plot_roc_data.csv", index=False)
    pr_plot_data.to_csv(CSV_DIR / "plot_pr_data.csv", index=False)
    df[["Files", "X", "f1_score"]].to_csv(CSV_DIR / "plot_f1_box_data.csv", index=False)

    top3 = summary["Files"].head(3).tolist()
    methods = summary["Files"].tolist()

    fig = plt.figure(figsize=(13, 7))
    for m in methods:
        alpha = 1.0 if m in top3 else 0.35
        lw = 2.5 if m in top3 else 1.1
        plt.plot(f1_matrix.index, f1_matrix[m], linewidth=lw, alpha=alpha, label=m)
    plt.xlabel("X")
    plt.ylabel("F1")
    plt.title("F1 vs X")
    plt.legend(ncol=2, fontsize=8)
    save_plot(fig, FIG_DIR / "fig_f1_vs_x.png")

    fig = plt.figure(figsize=(13, 7))
    for m in methods:
        alpha = 1.0 if m in top3 else 0.35
        lw = 2.5 if m in top3 else 1.1
        plt.plot(ba_matrix.index, ba_matrix[m], linewidth=lw, alpha=alpha, label=m)
    plt.xlabel("X")
    plt.ylabel("Balanced Accuracy")
    plt.title("Balanced Accuracy vs X")
    plt.legend(ncol=2, fontsize=8)
    save_plot(fig, FIG_DIR / "fig_balanced_accuracy_vs_x.png")

    fig = plt.figure(figsize=(12, 7))
    for m in methods:
        g = roc_plot_data[roc_plot_data["Files"] == m]
        alpha = 1.0 if m in top3 else 0.35
        lw = 2.5 if m in top3 else 1.1
        plt.plot(g["fpr"], g["TPR"], linewidth=lw, alpha=alpha, label=m)
    plt.xlabel("FPR")
    plt.ylabel("TPR")
    plt.title("ROC Curves")
    plt.legend(ncol=2, fontsize=8)
    save_plot(fig, FIG_DIR / "fig_roc_curves.png")

    fig = plt.figure(figsize=(12, 7))
    for m in methods:
        g = pr_plot_data[pr_plot_data["Files"] == m]
        alpha = 1.0 if m in top3 else 0.35
        lw = 2.5 if m in top3 else 1.1
        plt.plot(g["recall"], g["precision"], linewidth=lw, alpha=alpha, label=m)
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title("PR Curves")
    plt.legend(ncol=2, fontsize=8)
    save_plot(fig, FIG_DIR / "fig_pr_curves.png")

    box_df = df[["Files", "f1_score"]].copy()
    ordered_methods = summary["Files"].tolist()
    data_for_box = [box_df.loc[box_df["Files"] == m, "f1_score"].to_numpy() for m in ordered_methods]
    fig = plt.figure(figsize=(14, 7))
    plt.boxplot(data_for_box, tick_labels=ordered_methods, showfliers=False)
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("F1")
    plt.title("F1 Distribution by Method")
    save_plot(fig, FIG_DIR / "fig_f1_boxplot.png")

    friedman_stat, friedman_p = stats.friedmanchisquare(*[f1_matrix[c].to_numpy() for c in f1_matrix.columns])
    n_blocks = f1_matrix.shape[0]
    k_methods = f1_matrix.shape[1]
    kendall_w = float(friedman_stat / (n_blocks * (k_methods - 1)))

    ranks = f1_matrix.rank(axis=1, ascending=False, method="average")
    mean_ranks = ranks.mean(axis=0).sort_values()

    friedman_df = pd.DataFrame(
        [
            {
                "statistic": friedman_stat,
                "p_value": friedman_p,
                "n_blocks": n_blocks,
                "k_methods": k_methods,
                "kendall_w": kendall_w,
            }
        ]
    )
    friedman_df.to_csv(CSV_DIR / "friedman_result.csv", index=False)

    mean_rank_df = mean_ranks.reset_index()
    mean_rank_df.columns = ["Files", "mean_rank"]
    mean_rank_df.to_csv(CSV_DIR / "friedman_mean_ranks.csv", index=False)

    fig = plt.figure(figsize=(10, 6))
    plt.barh(mean_rank_df["Files"], mean_rank_df["mean_rank"])
    plt.xlabel("Mean Rank (lower is better)")
    plt.title(f"Friedman Mean Ranks (p={friedman_p:.3e}, W={kendall_w:.4f})")
    save_plot(fig, FIG_DIR / "fig_friedman_mean_ranks.png")

    best_series = f1_matrix[best_method].to_numpy()
    wilcoxon_rows = []
    pvals = []
    others = [m for m in f1_matrix.columns if m != best_method]

    for m in others:
        cur = f1_matrix[m].to_numpy()
        diff = best_series - cur
        nonzero = diff[diff != 0]
        if len(nonzero) == 0:
            stat_val = 0.0
            p_val = 1.0
        else:
            w_res = stats.wilcoxon(best_series, cur, alternative="greater", zero_method="wilcox", correction=False, mode="auto")
            stat_val = float(w_res.statistic)
            p_val = float(w_res.pvalue)
        pvals.append(p_val)
        wilcoxon_rows.append(
            {
                "best_method": best_method,
                "compared_method": m,
                "wilcoxon_statistic": stat_val,
                "p_value_raw": p_val,
                "hl_median_diff_f1": float(np.median(diff)),
                "rank_biserial": float(rank_biserial_from_paired(diff)),
                "best_higher_rate": float(np.mean(diff > 0)),
                "best_equal_rate": float(np.mean(diff == 0)),
                "best_lower_rate": float(np.mean(diff < 0)),
                "n_points": int(len(diff)),
                "n_nonzero_diff": int(len(nonzero)),
            }
        )

    wilcoxon_df = pd.DataFrame(wilcoxon_rows)
    if len(wilcoxon_df) > 0:
        wilcoxon_df["p_value_holm"] = holm_adjust(wilcoxon_df["p_value_raw"].to_numpy())
        wilcoxon_df["significant_best_higher_alpha_0_05"] = wilcoxon_df["p_value_holm"] < 0.05
    else:
        wilcoxon_df["p_value_holm"] = []
        wilcoxon_df["significant_best_higher_alpha_0_05"] = []

    wilcoxon_df = wilcoxon_df.sort_values(["p_value_holm", "hl_median_diff_f1"], ascending=[True, False]).reset_index(drop=True)
    wilcoxon_df.to_csv(CSV_DIR / "wilcoxon_vs_best.csv", index=False)

    if len(wilcoxon_df) > 0:
        plot_df = wilcoxon_df.sort_values("hl_median_diff_f1", ascending=True)
        colors = ["#1b9e77" if sig else "#7570b3" for sig in plot_df["significant_best_higher_alpha_0_05"]]
        fig = plt.figure(figsize=(10, 7))
        plt.barh(plot_df["compared_method"], plot_df["hl_median_diff_f1"], color=colors)
        plt.axvline(0, color="black", linewidth=1)
        plt.xlabel(f"HL Median Difference in F1 ({best_method} - method)")
        plt.title("Wilcoxon vs Best Method (Holm-corrected)")
        save_plot(fig, FIG_DIR / "fig_wilcoxon_vs_best_hl_diff.png")

    pair_rows = []
    pair_p = []
    pair_keys = []
    cols = list(f1_matrix.columns)
    for i in range(len(cols)):
        for j in range(i + 1, len(cols)):
            a = cols[i]
            b = cols[j]
            s1 = f1_matrix[a].to_numpy()
            s2 = f1_matrix[b].to_numpy()
            d = s1 - s2
            nz = d[d != 0]
            if len(nz) == 0:
                p_val = 1.0
                stat_val = 0.0
            else:
                res = stats.wilcoxon(s1, s2, alternative="two-sided", zero_method="wilcox", correction=False, mode="auto")
                p_val = float(res.pvalue)
                stat_val = float(res.statistic)
            pair_p.append(p_val)
            pair_keys.append((a, b))
            pair_rows.append(
                {
                    "method_a": a,
                    "method_b": b,
                    "wilcoxon_statistic": stat_val,
                    "p_value_raw": p_val,
                    "hl_median_diff_f1_a_minus_b": float(np.median(d)),
                    "rank_biserial_a_minus_b": float(rank_biserial_from_paired(d)),
                }
            )

    if len(pair_rows) > 0:
        pair_df = pd.DataFrame(pair_rows)
        pair_df["p_value_holm"] = holm_adjust(pair_df["p_value_raw"].to_numpy())
        pair_df["significant_alpha_0_05"] = pair_df["p_value_holm"] < 0.05
        pair_df.to_csv(CSV_DIR / "pairwise_wilcoxon_all_methods.csv", index=False)

        padj_matrix = pd.DataFrame(np.ones((len(cols), len(cols))), index=cols, columns=cols)
        for _, r in pair_df.iterrows():
            a = r["method_a"]
            b = r["method_b"]
            p = r["p_value_holm"]
            padj_matrix.loc[a, b] = p
            padj_matrix.loc[b, a] = p
        np.fill_diagonal(padj_matrix.values, 0.0)
        padj_matrix.to_csv(CSV_DIR / "pairwise_wilcoxon_holm_p_matrix.csv")

        fig = plt.figure(figsize=(11, 9))
        im = plt.imshow(padj_matrix.loc[summary["Files"], summary["Files"]], vmin=0, vmax=1, cmap="viridis_r")
        plt.colorbar(im, fraction=0.046, pad=0.04, label="Holm-adjusted p-value")
        plt.xticks(range(len(summary)), summary["Files"], rotation=45, ha="right")
        plt.yticks(range(len(summary)), summary["Files"])
        plt.title("Pairwise Wilcoxon Holm-adjusted p-value Matrix")
        save_plot(fig, FIG_DIR / "fig_pairwise_wilcoxon_holm_heatmap.png")

    co_best = [best_method]
    if len(wilcoxon_df) > 0:
        tie_methods = wilcoxon_df.loc[wilcoxon_df["significant_best_higher_alpha_0_05"] == False, "compared_method"].tolist()
        co_best = [best_method] + tie_methods

    decision_df = pd.DataFrame(
        [
            {
                "best_method": best_method,
                "co_best_methods": ";".join(sorted(co_best)),
                "co_best_count": len(co_best),
                "decision_rule": "Median_F1 > AUPRC > Q10_F1 > IQR_F1(asc) > Files(asc)",
            }
        ]
    )
    decision_df.to_csv(CSV_DIR / "final_decision.csv", index=False)

    ranking_plot_df = summary[["Files", "Median_F1", "AUPRC", "rank"]].copy().sort_values("Median_F1", ascending=True)
    ranking_plot_df.to_csv(CSV_DIR / "ranking_plot_data.csv", index=False)
    fig = plt.figure(figsize=(10, 7))
    colors = ["#d95f02" if m == best_method else "#7570b3" for m in ranking_plot_df["Files"]]
    plt.barh(ranking_plot_df["Files"], ranking_plot_df["Median_F1"], color=colors)
    plt.xlabel("Median F1 over X")
    plt.title("Method Ranking by Median F1")
    save_plot(fig, FIG_DIR / "fig_method_ranking_median_f1.png")

    scenario_specs = [
        {
            "scenario": "soil_release",
            "src": TARGET_FILE,
            "out_name": "fig_soil_release_8metrics_nejm_0504.pdf",
            "ordered_methods": summary["Files"].tolist(),
            "anchors": ["MMSage-noback", "MMSage"],
        },
        {
            "scenario": "simulation_release",
            "src": SIM_RELEASE_SRC,
            "out_name": "fig_simulation_release_8metrics_nejm_0504.pdf",
            "ordered_methods": summary["Files"].tolist(),
            "anchors": ["MMSage-noback", "MMSage"],
        },
        {
            "scenario": "simulation_consumption",
            "src": SIM_CONSUMPTION_SRC,
            "out_name": "fig_simulation_consumption_8metrics_nejm_0504.pdf",
            "ordered_methods": summary["Files"].tolist(),
            "anchors": ["MMSage-noback", "MMSage"],
        },
    ]

    all_metric_test_rows = []
    metric_plot_summary_rows = []
    metric_boxplot_data_rows = []
    pairwise_excel_rows = []

    for spec in scenario_specs:
        scenario_df = pd.read_csv(spec["src"])
        for c in ["X", "TP", "TN", "FP", "FN", "TPR", "TNR", "precision", "recall", "f1_score"]:
            scenario_df[c] = pd.to_numeric(scenario_df[c], errors="coerce")
        scenario_df = scenario_df.copy()
        scenario_df["Files"] = scenario_df["Files"].astype(str)
        scenario_df["Files"] = scenario_df["Files"].map(normalize_method_name)
        scenario_df = scenario_df.loc[~scenario_df["Files"].isin(["phi Clr", "rho Clr"])].copy()

        scenario_methods = scenario_df["Files"].unique().tolist()
        scenario_methods_norm = [normalize_method_name(m) for m in scenario_methods]
        scenario_ordered = [normalize_method_name(m) for m in spec["ordered_methods"] if normalize_method_name(m) in scenario_methods_norm]
        for anchor_name in ["MMSage", "MMSage-noback"]:
            resolved_anchor = resolve_anchor_method(scenario_methods, anchor_name)
            if resolved_anchor is not None:
                resolved_anchor = normalize_method_name(resolved_anchor)
                if resolved_anchor not in scenario_ordered:
                    scenario_ordered = [resolved_anchor] + scenario_ordered
        if len(scenario_ordered) == 0:
            scenario_ordered = sort_methods_for_plot(scenario_methods)
        else:
            scenario_ordered = sort_methods_for_plot(scenario_ordered)

        scenario_anchors = []
        for a in spec["anchors"]:
            resolved = resolve_anchor_method(scenario_methods, a)
            if resolved is not None:
                scenario_anchors.append(resolved)

        rows, summaries, box_rows, excel_rows = plot_nejm_3x4_metrics_boxplot(
            df=scenario_df,
            ordered_methods=scenario_ordered,
            existing_anchors=scenario_anchors,
            csv_dir=CSV_DIR,
            fig_dir=FIG_DIR,
            out_name=spec["out_name"],
        )
        all_metric_test_rows.extend(rows)
        metric_plot_summary_rows.extend(summaries)
        metric_boxplot_data_rows.extend(box_rows)
        pairwise_excel_rows.extend(excel_rows)

    if len(all_metric_test_rows) > 0:
        all_metric_df = pd.concat(all_metric_test_rows, ignore_index=True)
        all_metric_df.to_csv(CSV_DIR / "wilcoxon_vs_anchors_all_metrics_long.csv", index=False)
    else:
        pd.DataFrame(columns=[
            "metric", "direction", "best_method", "compared_method", "wilcoxon_statistic",
            "p_value_raw", "hl_median_diff_oriented", "rank_biserial", "best_better_rate", "best_equal_rate",
            "best_worse_rate", "n_points", "n_nonzero_diff", "p_value_holm", "significant_best_better_alpha_0_05",
            "metric_label", "scenario", "scenario_metric"
        ]).to_csv(CSV_DIR / "wilcoxon_vs_anchors_all_metrics_long.csv", index=False)

    if len(metric_boxplot_data_rows) > 0:
        all_box_df = pd.concat(metric_boxplot_data_rows, ignore_index=True)
        all_box_df.to_csv(CSV_DIR / "boxplot_data_all_metrics_long.csv", index=False)

        iqr_diag = (
            all_box_df.groupby(["metric", "metric_label", "Files"], as_index=False)
            .agg(
                n=("value", "size"),
                nunique=("value", "nunique"),
                q1=("value", lambda s: s.quantile(0.25)),
                median=("value", "median"),
                q3=("value", lambda s: s.quantile(0.75)),
                min_value=("value", "min"),
                max_value=("value", "max"),
            )
        )
        iqr_diag["iqr"] = iqr_diag["q3"] - iqr_diag["q1"]
        iqr_diag["is_degenerate_box"] = iqr_diag["iqr"] == 0
        iqr_diag.to_csv(CSV_DIR / "boxplot_iqr_diagnosis.csv", index=False)

    if len(all_metric_test_rows) > 0:
        ann_pairs = pd.concat(all_metric_test_rows, ignore_index=True).copy()
        ann_pairs["p_label"] = ann_pairs["p_value_holm"].apply(format_p_value)
        ann_pairs.to_csv(CSV_DIR / "annotation_pairs_all_metrics.csv", index=False)

    if len(pairwise_excel_rows) > 0:
        pairwise_df = pd.concat(pairwise_excel_rows, ignore_index=True)
        pairwise_df.to_excel(CSV_DIR / "pairwise_pvalues_all_scenarios.xlsx", index=False, engine="openpyxl")
    else:
        pd.DataFrame().to_excel(CSV_DIR / "pairwise_pvalues_all_scenarios.xlsx", index=False, engine="openpyxl")

    metric_plot_summary_df = pd.DataFrame(metric_plot_summary_rows)
    metric_plot_summary_df.to_csv(CSV_DIR / "multi_metrics_pdf_index.csv", index=False)
    if len(metric_plot_summary_df) > 0:
        metric_plot_summary_df[["scenario", "metric", "metric_label", "mmsage_vs_noback_p_raw"]].to_csv(
            CSV_DIR / "mmsage_vs_noback_pairwise_pvalues.csv",
            index=False,
        )

    print(f"Done. Best method: {best_method}")
    print(f"Target file: {TARGET_FILE}")
    print(f"CSV outputs: {CSV_DIR}")
    print(f"Figure outputs: {FIG_DIR}")


if __name__ == "__main__":
    main()
