#!/usr/bin/env python

import argparse
import polars as pl

def parse_args(args=None):
    description = "Convert TSV file into a compressed parquet."
    epilog = "Example usage: python tsv_to_parquet.py --help"

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-i",
        "--input",
        help="Path to input TSV file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Path to output parquet file.",
    )
    parser.add_argument('--version', action='version', version='1.0.0')
    return parser.parse_args(args)



def main(args=None):
    args = parse_args(args)

    # stream the TSV file and write to parquet
    (
        pl.scan_csv(args.input, separator='\t')
            .sink_parquet(args.output, compression='zstd')
    )

if __name__ == "__main__":
    main()
