from pathlib import Path

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
