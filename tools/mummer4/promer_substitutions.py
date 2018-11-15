#!/usr/bin/env python

import os
import subprocess
import sys

from Bio import SeqIO
from Bio.SeqIO import FastaIO

"""
# =============================================================================
Copyright Government of Canada 2018
Written by: Eric Marinier, Public Health Agency of Canada,
    National Microbiology Laboratory
Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this file except in compliance with the License. You may obtain a copy of the
License at:
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
# =============================================================================
"""

__version__ = '0.2.0'

"""
# =============================================================================
GLOBALS
# =============================================================================
"""

FASTA_DIRECTORY = "fasta"

REF_START = 2
REF_END = 3

HEADER_ROW = "[P1]\t[SUB]\t[SUB]\t[P2]\t[BUFF]\t[DIST]\
    \t[R]\t[Q]\t[FRM]\t[FRM]\t[TAG]\t[TAG]\n"

"""
# =============================================================================
REORIENT FILE
# =============================================================================
"""
def reorient_file(fasta_location, start, end):

    record = list(SeqIO.parse(fasta_location, "fasta"))[0]

    # reversed
    if start > end:
        record.seq = record.seq[(end - 1):start].reverse_complement()

    # same orientation
    else:
        record.seq = record.seq[(start - 1):end]

    SeqIO.write(record, fasta_location, "fasta")


"""
# =============================================================================
PROMER
# =============================================================================
"""
def promer(reference_location, query_location):

    # promer
    subprocess.check_output(
        ['promer', reference_location, query_location],
        universal_newlines=True)

    # filter
    output = subprocess.check_output(
        ['delta-filter', '-grq', 'out.delta'], universal_newlines=True)

    filter_file = open("out.filter", "w")
    filter_file.write(output)
    filter_file.close()


"""
# =============================================================================
MAIN
# =============================================================================
"""

# File locations from arguments:
reference_location = sys.argv[1]
query_location = sys.argv[2]

# Make directories:
if not os.path.exists(FASTA_DIRECTORY):
    os.mkdir(FASTA_DIRECTORY)

# Read query FASTA:
query = list(SeqIO.parse(query_location, "fasta"))
fasta_locations = []

# Prepare final output:
snps_file = open("snps.tab", "w")
snps_file.write(HEADER_ROW)

# Split the query FASTA file into multiple files, each with one record:
# (This is done to work around a bug in promer.)
for record in query:

    output_name = str(record.id) + ".fasta"
    output_location = os.path.join(FASTA_DIRECTORY, output_name)

    fasta_locations.append(output_location)
    SeqIO.write(record, output_location, "fasta")

# Run promer on each (new) FASTA file:
for fasta_location in fasta_locations:

    promer(reference_location, fasta_location)

    # show-coords
    output = subprocess.check_output(
        ['show-coords', '-THrcl', 'out.filter'], universal_newlines=True)

    if not output:
        continue

    ref_start = int(output.split()[REF_START])
    ref_end = int(output.split()[REF_END])

    reorient_file(fasta_location, ref_start, ref_end)
    promer(reference_location, fasta_location)

    # show-coords
    output = subprocess.check_output(
        ['show-coords', '-THrcl', 'out.filter'], universal_newlines=True)

    # show snps
    output = str(subprocess.check_output(
        ['show-snps', '-T', '-q', 'out.filter'], universal_newlines=True))
    output = output.split("\n")[4:]
    output = "\n".join(output)

    snps_file.write(output)

snps_file.close()

# Write all FASTA files into one output:
output_location = "contigs.fasta"
records = []

for fasta_location in fasta_locations:

    record = list(SeqIO.parse(fasta_location, "fasta"))[0]
    records.append(record)

contigs_file = open(output_location, 'w')
fasta_writer = FastaIO.FastaWriter(contigs_file, wrap=None)
fasta_writer.write_file(records)
contigs_file.close()
