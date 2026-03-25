import argparse

import pytest

from methylplotter.parsing import build_config, parse_gene, parse_region


def test_parse_gene_and_region():
    gene = parse_gene("chr1:100-200:GENE1")
    region = parse_region("chr1:90-210")

    assert gene.chrom == "chr1"
    assert gene.start == 100
    assert gene.end == 200
    assert gene.name == "GENE1"
    assert region.start == 90
    assert region.end == 210


def test_build_config_default_region_and_pairs():
    args = argparse.Namespace(
        platform="ont",
        bed=["a.bed", "b.bed"],
        sample=["s1", "s2"],
        gene="chr1:100-200:GENE1",
        region=None,
        output="out",
        line="marker,150",
        window_size=20,
        min_points_for_smooth=3,
        percent_col=None,
        verbose=False,
    )

    config = build_config(args)
    assert config.region.start == 0
    assert config.region.end == 700
    assert config.named_bed_paths == {"s1": "a.bed", "s2": "b.bed"}
    assert config.vline is not None
    assert config.vline.position == 150


def test_build_config_rejects_mismatched_samples_and_beds():
    args = argparse.Namespace(
        platform="ont",
        bed=["a.bed"],
        sample=["s1", "s2"],
        gene="chr1:100-200:GENE1",
        region=None,
        output="out",
        line=None,
        window_size=20,
        min_points_for_smooth=3,
        percent_col=None,
        verbose=False,
    )

    with pytest.raises(ValueError, match="must match"):
        build_config(args)


def test_build_config_rejects_duplicate_sample_names():
    args = argparse.Namespace(
        platform="ont",
        bed=["a.bed", "b.bed"],
        sample=["dup", "dup"],
        gene="chr1:100-200:GENE1",
        region=None,
        output="out",
        line=None,
        window_size=20,
        min_points_for_smooth=3,
        percent_col=None,
        verbose=False,
    )

    with pytest.raises(ValueError, match="unique"):
        build_config(args)
