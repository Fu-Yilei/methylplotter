#!/usr/bin/env python
import argparse
import sys, os
from src.__init__ import __version__ as version


def parse_args(argv):
    parser = argparse.ArgumentParser(
        prog="methylplotter",
        description="plot methylation from haplotype phased BED or BAM files",
        epilog=f"Version {version}",
    )
    parser.add_argument(
        "-v", "--version",
        action="version",
        version=f"methylplotter {version}"
    )
    parser.add_argument(
        "-p", "--platform", type=str, required=True,
        help="Sequencing platform: ont or pb"
    )
    parser.add_argument(
        "-b", "--bed", type=str, action="append", required=True,
        help="Input haplotype-phased BED file(s)"
    )
    parser.add_argument(
        "-s", "--sample", type=str, action="append", required=True,
        help="Sample name(s), in the same order as BED files"
    )
    parser.add_argument(
        "-g", "--gene", type=str, required=True,
        help="Gene info: format chr:start-end:name, e.g., chr15:80143550-80197576:FAH")
    parser.add_argument("-r", "--region", type=str, default=None,
                        help="Region to plot: format chr:start-end, e.g., chr15:80150000-80200000. "
                             "If not provided, use the gene coordinates +- 500 bps.")
    parser.add_argument(
        "-o", "--output", type=str, default="methylation_plot",
        help="Output prefix for plots (default: gene_methylation.pdf)"
    )
    parser.add_argument("-l", "--line", type=str, default=None,
                        help="A comma-separated (name,position) to draw vertical lines at (e.g., TR breakpoint)")
    parser.add_argument("-w", "--window_size", type=int, default=20,
                            help="Window size (bp) for smoothing (default: 20)")
    parser.add_argument("-m", "--min_points_for_smooth", type=int, default=3,
                            help="Minimum number of points required for smoothing (default: 3)")
    args = parser.parse_args(argv)
    return args
