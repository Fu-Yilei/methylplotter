from pathlib import Path

from methylplotter.core import draw_series, prepare_series
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


def test_prepare_series_pb_uses_percent_col_8():
    region = Region("chr5", 495, 565)
    named = {"PB": str(DATA_DIR / "pb_sample_a.bed")}

    series, raw = prepare_series(
        named,
        region,
        window_size=2,
        min_points_for_smooth=2,
        drop_first_col=False,
        percent_col=8,
    )

    assert len(raw["PB"]) == 4
    assert len(series[0][1]) == 3


def test_draw_series_writes_prefix_outputs(tmp_path):
    region = Region("chr1", 100, 200)
    series = [("A", [120, 140, 160], [20, 30, 40])]

    out_prefix = tmp_path / "plot_out"
    draw_series(series, region, out_path=str(out_prefix))

    assert (tmp_path / "plot_out.png").exists()
    assert (tmp_path / "plot_out.svg").exists()
