from pathlib import Path

import pytest

from methylplotter.cli import run


DATA_DIR = Path(__file__).parent / "data"


def test_cli_run_end_to_end(tmp_path):
    output = tmp_path / "figure"
    code = run(
        [
            "--platform",
            "ont",
            "--bed",
            str(DATA_DIR / "ont_sample_a.bed"),
            "--bed",
            str(DATA_DIR / "ont_sample_b.bed"),
            "--sample",
            "sampleA",
            "--sample",
            "sampleB",
            "--gene",
            "chr1:100-200:GENE1",
            "--region",
            "chr1:95-200",
            "--window_size",
            "3",
            "--min_points_for_smooth",
            "3",
            "--output",
            str(output),
        ]
    )

    assert code == 0
    assert (tmp_path / "figure.png").exists()
    assert (tmp_path / "figure.svg").exists()


def test_cli_pb_auto_detect_percent_col(tmp_path):
    output = tmp_path / "pb_figure"
    code = run(
        [
            "--platform",
            "pb",
            "--bed",
            str(DATA_DIR / "pb_sample_a.bed"),
            "--sample",
            "samplePB",
            "--gene",
            "chr5:500-560:GENE2",
            "--region",
            "chr5:495-565",
            "--output",
            str(output),
        ]
    )

    assert code == 0
    assert (tmp_path / "pb_figure.png").exists()
    assert (tmp_path / "pb_figure.svg").exists()


def test_cli_vline_annotation(tmp_path):
    output = tmp_path / "vline_figure"
    code = run(
        [
            "--platform",
            "ont",
            "--bed",
            str(DATA_DIR / "ont_sample_a.bed"),
            "--bed",
            str(DATA_DIR / "ont_sample_b.bed"),
            "--sample",
            "Sample A",
            "--sample",
            "Sample B",
            "--gene",
            "chr1:100-200:GENE1",
            "--region",
            "chr1:95-210",
            "--window_size",
            "3",
            "--min_points_for_smooth",
            "3",
            "--line",
            "Breakpoint,155",
            "--output",
            str(output),
        ]
    )

    assert code == 0
    assert (tmp_path / "vline_figure.png").exists()
    assert (tmp_path / "vline_figure.svg").exists()


def test_cli_errors_when_no_rows_in_region():
    with pytest.raises(SystemExit):
        run(
            [
                "--platform",
                "pb",
                "--bed",
                str(DATA_DIR / "pb_sample_a.bed"),
                "--sample",
                "samplePB",
                "--gene",
                "chr5:500-560:GENE2",
                "--region",
                "chr5:1-10",
            ]
        )
