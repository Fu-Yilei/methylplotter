from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
import numpy as np
import pandas as pd

from .model import Region

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_modkit_bed(path: str, *, drop_first_col: bool = False, percent_col: int = 10) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    if drop_first_col:
        df = df.iloc[:, 1:]

    required = max(2, percent_col)
    if df.shape[1] <= required:
        raise ValueError(
            f"File '{path}' does not have requested percent column index {percent_col}. "
            f"Detected {df.shape[1]} columns."
        )

    return pd.DataFrame(
        {
            "chr": df.iloc[:, 0].astype(str),
            "start": df.iloc[:, 1].astype(int),
            "end": df.iloc[:, 2].astype(int),
            "percent_modified": df.iloc[:, percent_col].astype(float),
        }
    )


def filter_region(df: pd.DataFrame, region: Region) -> pd.DataFrame:
    mask = (
        (df["chr"] == region.chrom)
        & (df["start"] >= region.start)
        & (df["end"] <= region.end)
    )
    return df.loc[mask].sort_values("start").reset_index(drop=True)


def sliding_mean(
    xs: np.ndarray,
    ys: np.ndarray,
    window_size: int,
    min_points_for_smooth: int = 5,
) -> Tuple[np.ndarray, np.ndarray]:
    if (
        len(xs) < window_size
        or window_size < 2
        or len(xs) < min_points_for_smooth
        or len(xs) != len(ys)
    ):
        return xs, ys

    order = np.argsort(xs)
    x_sorted = xs[order]
    y_sorted = ys[order]

    csum_y = np.cumsum(y_sorted, dtype=float)
    sums_y = csum_y[window_size - 1 :] - np.concatenate(([0.0], csum_y[:-window_size]))
    y_smooth = sums_y / window_size

    csum_x = np.cumsum(x_sorted, dtype=float)
    sums_x = csum_x[window_size - 1 :] - np.concatenate(([0.0], csum_x[:-window_size]))
    x_smooth = sums_x / window_size
    return x_smooth, y_smooth


def prepare_series(
    named_bed_paths: Dict[str, str],
    region: Region,
    window_size: int,
    min_points_for_smooth: int,
    *,
    drop_first_col: bool,
    percent_col: int,
) -> Tuple[List[Tuple[str, np.ndarray, np.ndarray]], Dict[str, pd.DataFrame]]:
    series: List[Tuple[str, np.ndarray, np.ndarray]] = []
    raw_region: Dict[str, pd.DataFrame] = {}

    for name, path in named_bed_paths.items():
        df = read_modkit_bed(path, drop_first_col=drop_first_col, percent_col=percent_col)
        subset = filter_region(df, region)
        raw_region[name] = subset

        if subset.empty:
            series.append((name, np.array([]), np.array([])))
            continue

        x = subset["start"].to_numpy()
        y = subset["percent_modified"].to_numpy()
        smoothed_x, smoothed_y = sliding_mean(x, y, window_size, min_points_for_smooth)
        series.append((name, smoothed_x, smoothed_y))

    return series, raw_region


def draw_series(
    series: List[Tuple[str, np.ndarray, np.ndarray]],
    region: Region,
    *,
    annotate_spans: Optional[List[Tuple[str, int, int]]] = None,
    annotate_vlines: Optional[List[Tuple[str, int]]] = None,
    figsize: Tuple[int, int] = (16, 6),
    linewidth: float = 2.0,
    grid: bool = True,
    font_label: int = 12,
    font_tick: int = 11,
    font_legend: int = 11,
    legend_loc: str = "best",
    out_path: Optional[str] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    fig, ax = plt.subplots(figsize=figsize)

    for name, xs, ys in series:
        xs = np.asarray(xs)
        ys = np.asarray(ys)
        if xs.size == 0:
            continue
        ax.plot(xs, ys, label=name, linewidth=linewidth)

    if annotate_spans:
        for label, start, end in annotate_spans:
            left, right = (start, end) if start <= end else (end, start)
            ax.axvspan(left, right, color="grey", alpha=0.12, label=label)

    if annotate_vlines:
        for label, position in annotate_vlines:
            ax.axvline(position, linestyle="--", linewidth=2, color="red", label=label)

    ax.set_xlabel("Genomic position", fontsize=font_label)
    ax.set_ylabel("% modified (windowed mean)", fontsize=font_label)
    ax.set_xlim(region.start, region.end)
    ax.set_ylim(0, 100)
    ax.tick_params(axis="both", labelsize=font_tick)

    handles, labels = ax.get_legend_handles_labels()
    seen = set()
    unique = [(h, l) for h, l in zip(handles, labels) if not (l in seen or seen.add(l))]
    if unique:
        ax.legend(*zip(*unique), loc=legend_loc, frameon=False, fontsize=font_legend)

    if grid:
        ax.grid(True, alpha=0.3)

    ax.get_xaxis().get_major_formatter().set_scientific(False)
    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    fig.tight_layout()

    if out_path:
        output = Path(out_path)
        if output.suffix:
            fig.savefig(output, bbox_inches="tight", pad_inches=0.1)
        else:
            fig.savefig(output.with_suffix(".png"), bbox_inches="tight", pad_inches=0.1)
            fig.savefig(output.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.1)

    return fig, ax
