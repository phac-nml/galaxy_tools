# Import dependancies needed
import pandas as pd
import numpy as np
import argparse

#Def the main function
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f',
        '--filename',
        required=True,       
        help= 'Specify your tsv input')
    parser.add_argument(
        '-o',
        '--output',
        default='output.csv',
        help='Specify output name')
    args = parser.parse_args()
    tsv_file = args.filename
    out_name = args.output

    no_comma_tsv = comma_remover(tsv_file)
    short_qc = qc_shortener(no_comma_tsv)
    tsv_to_csv(short_qc, out_name)

# Remove comma function
def comma_remover(tsv_file):
    # Create a table from the tsv file as an input into the dataframe.
    df = pd.read_csv(tsv_file,sep ='\t')
    # Change all commas to / in the QC message
    no_comma_tsv = df.replace(',', '/', regex=True)
    return no_comma_tsv

# Shorten QC results    
def qc_shortener(df):
    count = 0
    while count != len(df.index):
        message = str(df.at[count,'qc_message'])
        if len(message) > 150:
            results = message.find('|')
            new_message = "Truncated after first '|' : " + message[0:results]
            df['qc_message'] = df['qc_message'].replace(message, new_message)
        count += 1
    return df

# Converter function
def tsv_to_csv(df, out_name):
    # Write balues to a comma-separated values file
    df.to_csv(out_name,index=False)

if __name__ == '__main__':
    main()
