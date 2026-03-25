from pathlib import Path

from methylplotter.core import detect_percent_col, draw_series, prepare_series, prepare_series_detailed
from methylplotter.model import Region


DATA_DIR = Path(__file__).parent / "data"


def test_prepare_series_ont_uses_percent_col_10():
    region = Region("chr1", 95, 200)
    named = {"A": str(DATA_DIR / "ont_sample_a.bed")}

    series, raw = prepare_series(
        named,
        region,
        window_size=3,
        min_points_for_smooth=3,
        drop_first_col=False,
        percent_col=10,
    )

    assert "A" in raw
    assert len(raw["A"]) == 5
    label, xs, ys = series[0]
    assert label == "A"
    assert len(xs) == 3
    assert round(float(ys[0]), 2) == 20.0


def test_prepare_series_pb_uses_percent_col_3():
    region = Region("chr5", 495, 565)
    named = {"PB": str(DATA_DIR / "pb_sample_a.bed")}

    series, raw = prepare_series(
        named,
        region,
        window_size=2,
        min_points_for_smooth=2,
        drop_first_col=False,
        percent_col=3,
    )

    assert len(raw["PB"]) == 4
    assert len(series[0][1]) == 3


def test_detect_percent_col_for_platforms():
    pb_col = detect_percent_col(str(DATA_DIR / "pb_sample_a.bed"), "pb")
    ont_col = detect_percent_col(str(DATA_DIR / "ont_sample_a.bed"), "ont")
    assert pb_col == 3
    assert ont_col == 10


def test_prepare_series_detailed_auto_detects_col():
    region = Region("chr5", 495, 565)
    named = {"PB": str(DATA_DIR / "pb_sample_a.bed")}
    series, raw, used_cols = prepare_series_detailed(
        named,
        region,
        window_size=2,
        min_points_for_smooth=2,
        drop_first_col=False,
        percent_col=None,
        platform="pb",
    )

    assert len(series) == 1
    assert len(raw["PB"]) == 4
    assert used_cols["PB"] == 3


def test_draw_series_writes_prefix_outputs(tmp_path):
    region = Region("chr1", 100, 200)
    series = [("A", [120, 140, 160], [20, 30, 40])]

    out_prefix = tmp_path / "plot_out"
    draw_series(series, region, out_path=str(out_prefix))

    assert (tmp_path / "plot_out.png").exists()
    assert (tmp_path / "plot_out.svg").exists()


def test_draw_series_creates_missing_output_directory(tmp_path):
    region = Region("chr1", 100, 200)
    series = [("A", [120, 140, 160], [20, 30, 40])]

    out_path = tmp_path / "nested" / "folder" / "plot.png"
    draw_series(series, region, out_path=str(out_path))

    assert out_path.exists()
