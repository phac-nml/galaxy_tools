#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Version 1.01 - bugs kindly corrected by Jan van Haarst
# Modified by Matthew Gopez October 13th, 2017
# Rewritten by Matthew Gopez May 25th, 2020

import os
import subprocess
import argparse
import shutil
from pathlib import Path


PERL_OUT_FILES = ['stats.txt', 'sorted_contigs.fa', 'histogram_bins.dat.png',
                  'summed_contig_lengths.dat.png', 'histogram_bins.dat',
                  'summed_contig_lengths.dat']


def init_parser():
    """Create argument parser and return parser obj."""
    parser = argparse.ArgumentParser(description="usage: %prog [options]")

    parser.add_argument(
        "-d",
        "--working-dir",
        dest="working_dir",
        required=True)

    parser.add_argument(
        "-t",
        "--type",
        dest="file_type",
        required=True)

    parser.add_argument(
        "-b",
        "--bucket",
        dest="bucket",
        action='store_true')

    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        required=True)

    parser.add_argument(
        "-s",
        "--stats",
        dest="stats",
        required=True)

    parser.add_argument(
        "-sc",
        "--sorted-contigs",
        dest="sorted_contigs",
        required=True)

    parser.add_argument(
        "-hpng",
        "--histogram-png",
        dest="histogram_png",
        required=True)

    parser.add_argument(
        "-spng",
        "--summed-contigs-png",
        dest="summed_contigs_png",
        required=True)

    parser.add_argument(
        "-hd",
        "--histogram-data",
        dest="histogram_data",
        required=True)

    parser.add_argument(
        "-scd",
        "--summed-config-data",
        dest="summed_contig_data",
        required=True)

    return parser


def exec_fasta_summary(input_data, file_type, bucket, working_dir):
    """Execute fasta_summary.pl script with user arguments."""
    script_dir = Path(__file__).parent.absolute()

    if bucket:
        bucket_arg = '-b'
    else:
        bucket_arg = ''

    cli_command = '{}/fasta_summary.pl -i {} -t {} {} -o {} > /dev/null'.format(
        script_dir, input_data, file_type, bucket_arg, working_dir)

    try:
        subprocess.check_output(
            cli_command,
            stderr=subprocess.STDOUT,
            shell=True,
            universal_newlines=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError('Error running assembly_stats.py!\nReturn Code: {}\nOutput: {}'.format(
            exc.returncode, exc.output))


def main():
    """This is where the magic happens. (not really)

    1. Gets command line arguments.
    2. Grabs the user's desired parameters for running the perl script.
    3. Ensures the directories are in place.
    4. Executes fasta_summary.pl
    5. Move the out files from the perl script to the desired location the user specified.

    """
    parser = init_parser()
    args = parser.parse_args()

    working_dir = args.working_dir

    out_file_names = [args.stats, args.sorted_contigs, args.histogram_png,
                      args.summed_contigs_png, args.histogram_data, args.summed_contig_data]

    # Ensure working directory is created.
    Path(working_dir).mkdir(parents=True, exist_ok=True)

    # Execute Perl Script
    exec_fasta_summary(args.input, args.file_type, args.bucket, working_dir)

    # Rename out files to desired file names
    for perl_out_file, dest_file in zip(PERL_OUT_FILES, out_file_names):
        shutil.move(os.path.join(working_dir, perl_out_file),
                    dest_file)


if __name__ == "__main__":
    main()
