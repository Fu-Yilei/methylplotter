from __future__ import annotations

import gzip
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError

from .model import Region

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Okabe-Ito colorblind-safe palette
_PALETTE = [
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73",  # green
    "#CC79A7",  # pink
    "#56B4E9",  # sky blue
    "#D55E00",  # vermillion
    "#F0E442",  # yellow
    "#000000",  # black
]


def _first_data_tokens(path: str) -> List[str]:
    opener = gzip.open if str(path).endswith((".gz", ".bgz")) else open
    with opener(path, "rt") as handle:
        for line in handle:
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            return text.split("\t")
    raise ValueError(f"File '{path}' contains no data rows.")


def _is_percent_token(token: str) -> bool:
    text = token.strip()
    if text in {"", ".", "m"}:
        return False
    try:
        value = float(text)
    except ValueError:
        return False
    return 0.0 <= value <= 100.0


def detect_percent_col(path: str, platform: str, *, drop_first_col: bool = False) -> int:
    tokens = _first_data_tokens(path)
    if drop_first_col:
        tokens = tokens[1:]

    preferred = {
        "ont": [10, 3, 8],
        "pb": [3, 10, 8],
    }.get(platform.lower(), [10, 3, 8])

    for idx in preferred:
        if idx < len(tokens) and _is_percent_token(tokens[idx]):
            return idx

    for idx in range(3, len(tokens)):
        if _is_percent_token(tokens[idx]):
            return idx

    raise ValueError(
        f"Could not infer methylation percent column for '{path}' (platform={platform}, columns={len(tokens)})."
    )


def read_modkit_bed(path: str, *, drop_first_col: bool = False, percent_col: int = 10) -> pd.DataFrame:
    if drop_first_col:
        try:
            df = pd.read_csv(path, sep="\t", header=None, comment="#", compression="infer")
        except EmptyDataError as exc:
            raise ValueError(f"File '{path}' contains no parseable rows.") from exc
        df = df.iloc[:, 1:]
        if df.shape[1] <= max(2, percent_col):
            raise ValueError(
                f"File '{path}' does not have requested percent column index {percent_col}. "
                f"Detected {df.shape[1]} columns after dropping first column."
            )
        get_col = lambda idx: df.iloc[:, idx]
    else:
        usecols = sorted({0, 1, 2, percent_col})
        try:
            df = pd.read_csv(
                path,
                sep="\t",
                header=None,
                comment="#",
                compression="infer",
                low_memory=False,
                usecols=usecols,
            )
        except EmptyDataError as exc:
            raise ValueError(f"File '{path}' contains no parseable rows.") from exc
        except ValueError as exc:
            raise ValueError(
                f"File '{path}' does not have requested percent column index {percent_col}."
            ) from exc
        missing = [idx for idx in [0, 1, 2, percent_col] if idx not in df.columns]
        if missing:
            raise ValueError(
                f"File '{path}' is missing required column index/indices: {', '.join(map(str, missing))}."
            )
        get_col = lambda idx: df[idx]

    try:
        out = pd.DataFrame(
            {
                "chr": get_col(0).astype(str),
                "start": get_col(1).astype(int),
                "end": get_col(2).astype(int),
                "percent_modified": get_col(percent_col).astype(float),
            }
        )
    except (IndexError, KeyError, ValueError) as exc:
        raise ValueError(
            f"Failed to parse BED file '{path}'. Verify column layout and numeric methylation values."
        ) from exc

    return out


def filter_region(df: pd.DataFrame, region: Region) -> pd.DataFrame:
    mask = (df["chr"] == region.chrom) & (df["start"] >= region.start) & (df["end"] <= region.end)
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


def prepare_series_detailed(
    named_bed_paths: Dict[str, str],
    region: Region,
    window_size: int,
    min_points_for_smooth: int,
    *,
    drop_first_col: bool,
    percent_col: Optional[int] = None,
    platform: Optional[str] = None,
) -> Tuple[List[Tuple[str, np.ndarray, np.ndarray]], Dict[str, pd.DataFrame], Dict[str, int]]:
    if percent_col is None and platform is None:
        raise ValueError("platform is required when percent_col is not provided.")

    series: List[Tuple[str, np.ndarray, np.ndarray]] = []
    raw_region: Dict[str, pd.DataFrame] = {}
    used_percent_col: Dict[str, int] = {}

    for name, path in named_bed_paths.items():
        resolved_percent_col = percent_col
        if resolved_percent_col is None:
            resolved_percent_col = detect_percent_col(
                path,
                platform or "unknown",
                drop_first_col=drop_first_col,
            )

        used_percent_col[name] = resolved_percent_col
        df = read_modkit_bed(path, drop_first_col=drop_first_col, percent_col=resolved_percent_col)
        subset = filter_region(df, region)
        raw_region[name] = subset

        if subset.empty:
            series.append((name, np.array([]), np.array([])))
            continue

        x = subset["start"].to_numpy()
        y = subset["percent_modified"].to_numpy()
        smoothed_x, smoothed_y = sliding_mean(x, y, window_size, min_points_for_smooth)
        series.append((name, smoothed_x, smoothed_y))

    return series, raw_region, used_percent_col


def prepare_series(
    named_bed_paths: Dict[str, str],
    region: Region,
    window_size: int,
    min_points_for_smooth: int,
    *,
    drop_first_col: bool,
    percent_col: int,
) -> Tuple[List[Tuple[str, np.ndarray, np.ndarray]], Dict[str, pd.DataFrame]]:
    series, raw_region, _ = prepare_series_detailed(
        named_bed_paths,
        region,
        window_size,
        min_points_for_smooth,
        drop_first_col=drop_first_col,
        percent_col=percent_col,
    )
    return series, raw_region


def draw_series(
    series: List[Tuple[str, np.ndarray, np.ndarray]],
    region: Region,
    *,
    annotate_spans: Optional[List[Tuple[str, int, int]]] = None,
    annotate_vlines: Optional[List[Tuple[str, int]]] = None,
    title: Optional[str] = None,
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

    for i, (name, xs, ys) in enumerate(series):
        xs = np.asarray(xs)
        ys = np.asarray(ys)
        if xs.size == 0:
            continue
        color = _PALETTE[i % len(_PALETTE)]
        ax.plot(xs, ys, label=name, linewidth=linewidth, color=color, zorder=3)
        ax.fill_between(xs, ys, alpha=0.08, color=color, zorder=2)

    if annotate_spans:
        for label, start, end in annotate_spans:
            left, right = (start, end) if start <= end else (end, start)
            ax.axvspan(left, right, color="#6BAED6", alpha=0.15, zorder=0)
            mid = (left + right) / 2
            ax.text(
                mid, 97, label,
                ha="center", va="top",
                fontsize=font_tick - 1, color="#2166AC",
                fontweight="bold", clip_on=True,
            )

    if annotate_vlines:
        for label, position in annotate_vlines:
            ax.axvline(position, linestyle="--", linewidth=1.5, color="#D62728", label=label, zorder=4)

    ax.set_xlabel("Genomic position", fontsize=font_label)
    ax.set_ylabel("% methylation (windowed mean)", fontsize=font_label)
    ax.set_xlim(region.start, region.end)
    ax.set_ylim(0, 100)
    ax.tick_params(axis="both", labelsize=font_tick)

    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    handles, labels = ax.get_legend_handles_labels()
    seen = set()
    unique = [(h, l) for h, l in zip(handles, labels) if not (l in seen or seen.add(l))]
    if unique:
        ax.legend(*zip(*unique), loc=legend_loc, frameon=False, fontsize=font_legend)

    if grid:
        ax.grid(True, alpha=0.2, linewidth=0.5)

    if title:
        ax.set_title(title, fontsize=font_label + 1, pad=8)

    fig.tight_layout()

    if out_path:
        output = Path(out_path)
        output.parent.mkdir(parents=True, exist_ok=True)

        if output.suffix:
            fig.savefig(output, bbox_inches="tight", pad_inches=0.1)
        else:
            fig.savefig(output.with_suffix(".png"), bbox_inches="tight", pad_inches=0.1)
            fig.savefig(output.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.1)

    return fig, ax
