#!/usr/bin/env python

import argparse
import glob

import polars as pl

def parse_args(args=None):
    description = "Combine multiple parquet files having the same structure into one."
    epilog = "Example usage: python combine_anis.py --help"

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-i",
        "--input",
        help="Path to input parquet files.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Path to combined TSV file.",
    )
    parser.add_argument('--version', action='version', version='1.0.0')
    return parser.parse_args(args)



def main(args=None):
    args = parse_args(args)

    df_lst = []

    for file in glob.glob(args.input):
        df = (
            pl.read_parquet(file)
                .with_columns([
                    pl.col('query').cast(pl.String),
                    pl.col('reference').cast(pl.String),
                    pl.col('tani').cast(pl.Float32),
                    pl.col('gani').cast(pl.Float32),
                    pl.col('ani').cast(pl.Float32),
                    pl.col('qcov').cast(pl.Float32),
                    pl.col('rcov').cast(pl.Float32)
                ])
        )
        df_lst.append(df)

    pl.concat(df_lst).write_csv(args.output, separator='\t')

if __name__ == "__main__":
    main()
