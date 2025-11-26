#!/usr/bin/env python

import argparse
import polars as pl

def parse_args(args=None):
    description = "Convert parquet file into a TSV."
    epilog = "Example usage: python parquet_to_tsv.py --help"

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-i",
        "--input",
        help="Path to input parquet file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Path to output TSV file.",
    )
    parser.add_argument(
        "-c",
        "--columns",
        help="Columns to write out to TSV.",
        nargs='+'
    )
    parser.add_argument(
        "-a",
        "--header",
        help="Whether to include header in output.",
        type=bool,
        default=False
    )
    parser.add_argument('--version', action='version', version='1.0.0')
    return parser.parse_args(args)



def main(args=None):
    args = parse_args(args)

    # stream the parquet file and write to TS
    (
        pl.scan_parquet(args.input).select(args.columns) if args.columns else pl.scan_parquet(args.input)
            .sink_csv(args.output, separator='\t', include_header=args.header)
    )

if __name__ == "__main__":
    main()
