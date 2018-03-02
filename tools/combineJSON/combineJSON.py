#!/usr/bin/env python
import json
import argparse
import sys


def init_parser():
    parser = argparse.ArgumentParser(
        prog="combineJSON",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Combine JSON data arrays into a single array")
    parser.add_argument('-i',
                        nargs='*',
                        help="Input JSON files to be combined")
    parser.add_argument('-o',
                        help='Output file name')
    return parser


parser = init_parser()
args = parser.parse_args()
input_files = args.i
json_file = []

if input_files is None or len(input_files) < 2:
    print('Not enough input files. '
          'Please use -i filename.txt filename1.txt '
          'to combine JSON data files')
    sys.exit(0)


for file_path in input_files:
    try:
        # Attempt to open each input file, parse as JSON
        with open(file_path, 'r') as curr_file:
            file_data = curr_file.read()
            parsed_json_file = json.loads(file_data)
            # Append each valid JSON data array
            for entry in parsed_json_file:
                json_file.append(entry)
    except Exception as e:
        print("Help! I can't parse this file {}. "
              "Are you sure this is a valid JSON file?"
              .format(file_path))
        raise(e)

if args.o:
    with open(args.o, 'w') as out_json:
        json.dump(json_file, out_json)
else:
    print(json.dumps(json_file))
