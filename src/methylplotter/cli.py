from __future__ import annotations

import argparse
import sys

from . import __version__
from .core import draw_series, prepare_series_detailed
from .parsing import build_config


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="methylplotter",
        description="Plot methylation from haplotype-phased BED files.",
    )
    parser.add_argument("-v", "--version", action="version", version=f"methylplotter {__version__}")
    parser.add_argument("-p", "--platform", required=True, help="Sequencing platform: ont or pb")
    parser.add_argument("-b", "--bed", action="append", required=True, help="Input BED file(s)")
    parser.add_argument(
        "-s",
        "--sample",
        action="append",
        required=True,
        help="Sample name(s), in the same order as BED files",
    )
    parser.add_argument(
        "-g",
        "--gene",
        required=True,
        help="Gene info: chr:start-end:name, e.g. chr15:80143550-80197576:FAH",
    )
    parser.add_argument(
        "-r",
        "--region",
        default=None,
        help="Region to plot: chr:start-end. Default uses gene coordinates +/- 500 bp.",
    )
    parser.add_argument("-o", "--output", default="methylation_plot", help="Output file or prefix")
    parser.add_argument(
        "-l",
        "--line",
        default=None,
        help="A comma-separated (name,position) to draw vertical lines",
    )
    parser.add_argument(
        "-w",
        "--window_size",
        type=int,
        default=20,
        help="Window size for smoothing (number of points, not bp).",
    )
    parser.add_argument(
        "-m",
        "--min_points_for_smooth",
        type=int,
        default=3,
        help="Minimum points required for smoothing",
    )
    parser.add_argument(
        "--percent_col",
        type=int,
        default=None,
        help="Override methylation percent column index (0-based). Default uses platform-aware detection.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print per-sample parsing and smoothing details to stderr.",
    )
    return parser


def run(argv: list[str]) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        config = build_config(args)
    except ValueError as exc:
        parser.error(str(exc))

    try:
        series, raw_region, used_percent_cols = prepare_series_detailed(
            config.named_bed_paths,
            config.region,
            window_size=config.window_size,
            min_points_for_smooth=config.min_points_for_smooth,
            drop_first_col=False,
            percent_col=config.percent_col,
            platform=config.platform,
        )
    except (OSError, ValueError) as exc:
        parser.error(str(exc))

    if all(df.empty for df in raw_region.values()):
        parser.error(
            "No rows were found in the selected region for any sample. "
            "Check --region/--gene coordinates or input files."
        )

    vlines = [(config.vline.name, config.vline.position)] if config.vline else None

    fig, _ = draw_series(
        series,
        config.region,
        annotate_spans=[(config.gene.name, config.gene.start, config.gene.end)],
        annotate_vlines=vlines,
        title=f"{config.gene.name} — methylation profile",
        out_path=config.output,
    )
    fig.clear()

    if config.verbose:
        smoothed_counts = {name: len(xs) for name, xs, _ in series}
        for name, df in raw_region.items():
            points = len(df)
            smoothed = smoothed_counts.get(name, 0)
            percent_col = used_percent_cols[name]
            print(
                f"[methylplotter] sample={name} points={points} smoothed_points={smoothed} "
                f"percent_col={percent_col}",
                file=sys.stderr,
            )

    return 0


def main() -> None:
    raise SystemExit(run(sys.argv[1:]))


if __name__ == "__main__":
    main()
