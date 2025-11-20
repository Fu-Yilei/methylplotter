import os
import sys
from src.ont_bed_parsing import draw_series, prepare_series
from src.parse_args import parse_args  # assuming you defined parse_args in args.py


def main(argv):
    args = parse_args(argv)
    gene_info = args.gene.split(":")

    gene_chr = gene_info[0]
    gene_start_ = int(gene_info[1].split("-")[0])
    gene_end_ = int(gene_info[1].split("-")[1])
    gene_name_ = gene_info[2]
    if args.line:
        line = [i for i in args.line.split(",")]
        

    if args.region:
        region_info = args.region.split(":")
        region_chr = region_info[0]
        region_start_ = int(region_info[1].split("-")[0])
        region_end_ = int(region_info[1].split("-")[1])
        region = (region_chr, region_start_, region_end_)
    else:
        region = (gene_chr, gene_start_-500, gene_end_+500)

    if args.platform == "ont":
        named_bed_paths = {name: path for name, path in zip(args.sample, args.bed)}
        series, raw_region = prepare_series(
            named_bed_paths,
            region,
            window_size=20,
            min_points_for_smooth=5,
            drop_first_col=False,
            percent_col=10
        )
        draw_series(
            series,
            region,
            annotate_spans=[(gene_name_, gene_start_, gene_end_)],
            annotate_vlines=[(line[0], int(line[1]))] if args.line else None,
            out_path = f"{args.output}"
        )
    elif args.platform == "pb":
        named_bed_paths = {name: path for name, path in zip(args.sample, args.bed)}
        series, raw_region = prepare_series(
            named_bed_paths,
            region,
            window_size=20,
            min_points_for_smooth=5,
            drop_first_col=False,
            percent_col=8
        )
        draw_series(
            series,
            region,
            annotate_spans=[(gene_name_, gene_start_, gene_end_)],
            annotate_vlines=[(line[0], int(line[1]))] if args.line else None,
            out_path = f"{args.output}"
        )

if __name__ == "__main__":
    main(sys.argv[1:])