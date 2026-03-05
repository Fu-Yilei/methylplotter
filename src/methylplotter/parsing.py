from __future__ import annotations

from argparse import Namespace
from typing import Dict, Optional

from .model import Gene, Region, RunConfig, VerticalLine


PLATFORMS = {"ont": 10, "pb": 8}


def parse_gene(gene_text: str) -> Gene:
    try:
        chrom, coordinates, name = gene_text.split(":", 2)
        start_s, end_s = coordinates.split("-", 1)
        start, end = int(start_s), int(end_s)
    except (AttributeError, ValueError) as exc:
        raise ValueError(
            "Gene must be formatted as chr:start-end:name (example: chr15:80143550-80197576:FAH)."
        ) from exc

    _validate_region_bounds(start, end, "Gene")
    return Gene(chrom=chrom, start=start, end=end, name=name)


def parse_region(region_text: str) -> Region:
    try:
        chrom, coordinates = region_text.split(":", 1)
        start_s, end_s = coordinates.split("-", 1)
        start, end = int(start_s), int(end_s)
    except (AttributeError, ValueError) as exc:
        raise ValueError(
            "Region must be formatted as chr:start-end (example: chr15:80150000-80200000)."
        ) from exc

    _validate_region_bounds(start, end, "Region")
    return Region(chrom=chrom, start=start, end=end)


def parse_vertical_line(line_text: str) -> VerticalLine:
    try:
        name, position_s = line_text.split(",", 1)
        position = int(position_s)
    except (AttributeError, ValueError) as exc:
        raise ValueError(
            "Line must be formatted as name,position (example: TR breakpoint,80170000)."
        ) from exc
    return VerticalLine(name=name, position=position)


def build_named_bed_paths(samples: list[str], beds: list[str]) -> Dict[str, str]:
    if len(samples) != len(beds):
        raise ValueError(
            f"The number of samples ({len(samples)}) must match the number of BED files ({len(beds)})."
        )
    return dict(zip(samples, beds))


def build_config(args: Namespace) -> RunConfig:
    platform = args.platform.lower()
    if platform not in PLATFORMS:
        allowed = ", ".join(sorted(PLATFORMS))
        raise ValueError(f"Unsupported platform '{args.platform}'. Allowed values: {allowed}.")

    if args.window_size < 1:
        raise ValueError("--window_size must be >= 1.")
    if args.min_points_for_smooth < 1:
        raise ValueError("--min_points_for_smooth must be >= 1.")

    gene = parse_gene(args.gene)
    region = (
        parse_region(args.region)
        if args.region
        else Region(gene.chrom, max(gene.start - 500, 0), gene.end + 500)
    )
    vline = parse_vertical_line(args.line) if args.line else None

    if region.chrom != gene.chrom:
        raise ValueError(
            f"Region chromosome ({region.chrom}) and gene chromosome ({gene.chrom}) must match."
        )

    named_bed_paths = build_named_bed_paths(args.sample, args.bed)

    return RunConfig(
        platform=platform,
        named_bed_paths=named_bed_paths,
        gene=gene,
        region=region,
        output=args.output,
        window_size=args.window_size,
        min_points_for_smooth=args.min_points_for_smooth,
        vline=vline,
    )


def percent_col_for_platform(platform: str) -> int:
    return PLATFORMS[platform]


def _validate_region_bounds(start: int, end: int, label: str) -> None:
    if start < 0 or end < 0:
        raise ValueError(f"{label} positions must be non-negative.")
    if end <= start:
        raise ValueError(f"{label} end must be greater than start.")
