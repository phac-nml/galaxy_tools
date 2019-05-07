#!/usr/bin/env python

# Import dependancies needed
import argparse

import pandas as pd

# Define the main function:


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f',
        '--filename',
        required=True,
        help='Specify your biohansel tsv or other tabular separated input')
    parser.add_argument(
        '-o',
        '--output',
        default='output.csv',
        help='Specify output name')
    args = parser.parse_args()
    tsv_file = args.filename
    out_name = args.output

    df_input = pd.read_csv(tsv_file, sep='\t')

    df_no_comma = df_input.replace(',', '/', regex=True)
    df = qc_shortener(df_no_comma)
    df.to_csv(out_name, index=False)

# Shorten QC results:


def qc_shortener(df):
    for i, row in df.iterrows():
        message = row['qc_message']
        try:
            if len(message) > 150:
                df.at[i,'qc_message'] = message[0:150]
                df.at[i, 'qc_message_2'] = message[150:]
                if len(message) > 300:
                    df.at[i, 'qc_message_3'] = message[300:450]
        except TypeError:
            pass
    return df


if __name__ == '__main__':
    main()
