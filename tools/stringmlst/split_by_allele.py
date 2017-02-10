#!/usr/bin/env python
import os
import getopt
import sys
from Bio import SeqIO

ERROR_MSG = "Error could not parse out allele name and number from '%s'"


def split_allele_file(alleles, profiles):

    writers = {}

    handle = open(alleles, "rU")
    for record in SeqIO.parse(handle, "fasta"):

        seqid = record.id

        # split out the alelle name from the version number
        # attempting to split based on '-' first, if that fails, then '_'
        result = seqid.split('_')

        if len(result) != 2:
            result = seqid.split('-')
            if len(result) == 2:
                newid = '_'.join(result)
                record.id = newid
            else:
                print ERROR_MSG % seqid
                exit(0)

        name, num = result

        # if writer exist, then write to that current fasta file
        if name in writers:
            SeqIO.write(record, writers[name], "fasta")
        else:
            # new allele found, create new writer and add the first record
            file_name = name + '.fasta'
            output_fh = open(file_name, "w")
            SeqIO.write(record, output_fh, "fasta")
            writers[name] = output_fh

    handle.close()

    # create config file based on the alleles found
    with open('config.txt', 'w') as cfile:
        cfile.write("[loci]\n")
        for name, writer in writers.iteritems():
            path = os.path.realpath(writer.name)
            cfile.write("%s\t%s\n" % (name, path))
        cfile.write("[profile]\n")
        cfile.write("profile\t%s\n" % profiles)

    return


alleles = None
profiles = None

"""Input arguments"""
options, remainder = getopt.getopt(sys.argv[1:], '', [
 'alleles=',
 'profiles='
])

for opt, arg in options:
    if opt in ('--alleles'):
        alleles = arg
    elif opt in ('--profiles'):
        profiles = arg

if alleles and profiles:
    split_allele_file(alleles, profiles)
