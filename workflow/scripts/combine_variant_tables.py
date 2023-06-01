#!/usr/bin/env python3

import pandas as pd
import pathlib
import shutil

def main(args):
    df = pd.DataFrame()
    for filepath in args.input:
        sample_name = filepath.stem
        sample_df = pd.read_csv(filepath, sep='\t')
        sample_df['Sample'] = sample_name
        df = pd.concat([df, sample_df])
    df_grouped = df.groupby(['Sample'] + args.fields).size().reset_index(name='Count')
    if args.mqc is not None:
        shutil.copy(args.mqc, args.output)
        df_wide = pd.pivot(data=df_grouped, index='Sample', columns='FILTER', values='Count').reset_index().fillna(0)
        df_wide.to_csv(args.output, mode='a', sep='\t', index=False)
    else:
        df_grouped.to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input",
                        type=pathlib.Path,
                        required=True,
                        metavar="PATH",
                        nargs="+")
    parser.add_argument("-o", "--output",
                        type=str,
                        required=True,
                        metavar="PATH")
    parser.add_argument("-f", "--fields",
                        type=str,
                        required=True,
                        metavar="FIELD",
                        nargs="+")
    parser.add_argument("--mqc",
                        help="Append data to this MultiQC header",
                        type=pathlib.Path)
        
    args = parser.parse_args()

    main(args)