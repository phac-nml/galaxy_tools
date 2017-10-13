#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Version 1.01 - bugs kindly corrected by Jan van Haarst
# Modified by Matthew Gopez October 13th, 2017

import logging
import os
import subprocess
import sys


log = logging.getLogger(__name__)

assert sys.version_info[:2] >= (2, 4)


def stop_err(msg):
    sys.stderr.write('%s\n' % msg)
    sys.exit()


def __main__():

    # Parse Command Line

    working_dir = sys.argv[2]
    type = sys.argv[3]
    bucket = sys.argv[4]
    input = sys.argv[5]
    stats = sys.argv[6]
    sortedcontigs = sys.argv[7]
    histogrampng = sys.argv[8]
    summedcontigspng = sys.argv[9]
    histogramdata = sys.argv[10]
    summedcontigdata = sys.argv[11]
    try:  # for test - needs this done
        os.makedirs(working_dir)
    except Exception, e:
        stop_err('Error running assembly_stats_txt.py ' + str(e))

    cmdline = '%s/fasta_summary.pl -i %s -t %s %s -o %s > /dev/null' \
        % (os.path.dirname(sys.argv[0]), input, type, bucket,
           working_dir)
    try:
        proc = subprocess.Popen(args=cmdline, shell=True,
                                stderr=subprocess.PIPE)
        returncode = proc.wait()

        # get stderr, allowing for case where it's very large

        stderr = ''
        buffsize = 1048576
        try:
            while True:
                stderr += proc.stderr.read(buffsize)
                if not stderr or len(stderr) % buffsize != 0:
                    break
        except OverflowError:
            pass
        if returncode != 0:
            raise Exception
    except Exception, e:
        stop_err('Error running assembly_stats.py ' + str(e))

    stats_path = os.path.join(working_dir, 'stats.txt')
    sorted_contigs_path = os.path.join(working_dir, 'sorted_contigs.fa')
    histogram_png_path = os.path.join(working_dir,
                                      'histogram_bins.dat.png')
    summed_contigs_path = os.path.join(working_dir,
                                       'summed_contig_lengths.dat.png')
    histogram_data_path = os.path.join(working_dir, 'histogram_bins.dat')
    summed_contigs_data_path = os.path.join(working_dir,
                                            'summed_contig_lengths.dat')

    out = open(stats, 'w')
    for line in open(stats_path):
        out.write('%s' % line)
    out.close()

    out = open(sortedcontigs, 'w')
    for line in open(sorted_contigs_path):
        out.write('%s' % line)
    out.close()

    out = open(histogrampng, 'w')
    for line in open(histogram_png_path):
        out.write('%s' % line)
    out.close()

    out = open(summedcontigspng, 'w')
    for line in open(summed_contigs_path):
        out.write('%s' % line)
    out.close()

    out = open(histogramdata, 'w')
    for line in open(histogram_data_path):
        out.write('%s' % line)
    out.close()

    out = open(summedcontigdata, 'w')
    for line in open(summed_contigs_data_path):
        out.write('%s' % line)
    out.close()


if __name__ == '__main__':
    __main__()
