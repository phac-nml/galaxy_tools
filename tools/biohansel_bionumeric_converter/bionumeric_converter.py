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
        help='Specify your tsv input')
    parser.add_argument(
        '-o',
        '--output',
        default='output.csv',
        help='Specify output name')
    args = parser.parse_args()
    tsv_file = args.filename
    out_name = args.output

    no_comma_tsv = comma_remover(tsv_file)
    df = qc_shortener(no_comma_tsv)
    df.to_csv(out_name, index=False)

# Remove comma function:


def comma_remover(tsv_file):
    # Create a table from the tsv file as an input into the dataframe.
    df = pd.read_csv(tsv_file, sep='\t')
    # Change all commas to / in the QC message
    no_comma_tsv = df.replace(',', '/', regex=True)
    return no_comma_tsv

# Shorten QC results:


def qc_shortener(df):
    for count in df.index:
        message = str(df.at[count, 'qc_message'])
        if len(message) > 150:
            results = message.find('|')
            new_message = "Truncated after first '|' : " + message[0:results]
            df['qc_message'] = df['qc_message'].replace(message, new_message)
    return df


if __name__ == '__main__':
    main()
