from __future__ import annotations

import argparse
import sys

from . import __version__
from .core import draw_series, prepare_series
from .parsing import build_config, percent_col_for_platform


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
    parser.add_argument("-w", "--window_size", type=int, default=20, help="Window size for smoothing")
    parser.add_argument(
        "-m",
        "--min_points_for_smooth",
        type=int,
        default=3,
        help="Minimum points required for smoothing",
    )
    return parser


def run(argv: list[str]) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        config = build_config(args)
    except ValueError as exc:
        parser.error(str(exc))

    percent_col = percent_col_for_platform(config.platform)
    series, _ = prepare_series(
        config.named_bed_paths,
        config.region,
        window_size=config.window_size,
        min_points_for_smooth=config.min_points_for_smooth,
        drop_first_col=False,
        percent_col=percent_col,
    )

    vlines = [(config.vline.name, config.vline.position)] if config.vline else None

    draw_series(
        series,
        config.region,
        annotate_spans=[(config.gene.name, config.gene.start, config.gene.end)],
        annotate_vlines=vlines,
        out_path=config.output,
    )
    return 0


def main() -> None:
    raise SystemExit(run(sys.argv[1:]))


if __name__ == "__main__":
    main()
