from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple


@dataclass(frozen=True)
class Region:
    chrom: str
    start: int
    end: int

    def as_tuple(self) -> Tuple[str, int, int]:
        return (self.chrom, self.start, self.end)


@dataclass(frozen=True)
class Gene:
    chrom: str
    start: int
    end: int
    name: str


@dataclass(frozen=True)
class VerticalLine:
    name: str
    position: int


@dataclass(frozen=True)
class RunConfig:
    platform: str
    named_bed_paths: Dict[str, str]
    gene: Gene
    region: Region
    output: str
    window_size: int
    min_points_for_smooth: int
    vline: Optional[VerticalLine]
    percent_col: Optional[int]
    verbose: bool
