#!/usr/bin/env python
#Version 1.01 - bugs kindly corrected by Jan van Haarst
import pkg_resources
import logging, os, string, sys, tempfile, glob, shutil, types, urllib
import shlex, subprocess
from optparse import OptionParser, OptionGroup
from stat import *


log = logging.getLogger( __name__ )

assert sys.version_info[:2] >= ( 2, 4 )

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()

def __main__():
    #Parse Command Line
    s = 'assembly_stats_txt.py:  argv = %s\n' % (sys.argv)
    argcnt = len(sys.argv)
    html_file = sys.argv[1]
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
    try: # for test - needs this done
        os.makedirs(working_dir)
    except Exception, e:
        stop_err( 'Error running assembly_stats_txt.py ' + str( e ) )
  
    
    cmdline = '%s/fasta_summary.pl -i %s -t %s %s -o %s > /dev/null' % (os.path.dirname(sys.argv[0]),input, type, bucket, working_dir)
    try:
        proc = subprocess.Popen( args=cmdline, shell=True, stderr=subprocess.PIPE )
        returncode = proc.wait()
        # get stderr, allowing for case where it's very large
        stderr = ''
        buffsize = 1048576
        try:
            while True:
                stderr += proc.stderr.read( buffsize )
                if not stderr or len( stderr ) % buffsize != 0:
                    break
        except OverflowError:
            pass
        if returncode != 0:
            raise Exception, stderr
    except Exception, e:
        stop_err( 'Error running assembly_stats.py ' + str( e ) )
 
    stats_path = os.path.join(working_dir,'stats.txt')
    sorted_contigs_path = os.path.join(working_dir,'sorted_contigs.fa')
    histogram_png_path = os.path.join(working_dir,'histogram_bins.dat.png')
    summed_contigs_path = os.path.join(working_dir,'summed_contig_lengths.dat.png')
    histogram_data_path =  os.path.join(working_dir,'histogram_bins.dat')
    summed_contigs_data_path = os.path.join(working_dir,'summed_contig_lengths.dat')

    out = open(stats,'w')
    for line in open( stats_path ):
        out.write( "%s" % (line) )
    out.close()

    out = open(sortedcontigs,'w')
    for line in open(sorted_contigs_path ):
        out.write( "%s" % (line) )
    out.close()

    out = open(histogrampng,'w')
    for line in open(histogram_png_path ):
        out.write( "%s" % (line) )
    out.close()

    out = open(summedcontigspng,'w')
    for line in open(summed_contigs_path ):
        out.write( "%s" % (line) )
    out.close()


    out = open(histogramdata,'w')
    for line in open(histogram_data_path ):
        out.write( "%s" % (line) )
    out.close()

    out = open(summedcontigdata,'w')
    for line in open(summed_contigs_data_path ):
        out.write( "%s" % (line) )
    out.close()









#    rval = ['<html><head><title>Assembly stats Galaxy Composite Dataset </title></head><p/>']
#    rval.append('<div>%s<p/></div>' % (cmdline) )
#    rval.append('<div>This composite dataset is composed of the following files:<p/><ul>')
#    rval.append( '<li><a href="%s" type="text/plain">%s </a>%s</li>' % (stats_path,'stats.txt','stats.txt' ) )
#    rval.append( '<li><a href="%s" type="text/plain">%s </a>%s</li>' % (sorted_contigs_path,'sorted_contigs.fa','sorted_contigs.fa' ) )
#    rval.append( '<li><a href="%s" type="image/png">%s </a>%s</li>' % (histogram_png_path,'histogram_bins.dat.png','histogram_bins.dat.png' ) )
#    rval.append( '<li><a href="%s" type="image/png">%s </a>%s</li>' % (summed_contigs_path,'summed_contig_lengths.dat.png','summed_contig_lengths.dat.png' ) )
#    rval.append( '<li><a href="%s" type="text/plain">%s </a>%s</li>' % (histogram_data_path,'histogram_bins.dat','histogram_bins.dat' ) )	
#    rval.append( '<li><a href="%s" type="text/plain">%s </a>%s</li>' % (summed_contigs_data_path,'summed_contig_lengths.dat','summed_contig_lengths.dat' ) )


#	
#    rval.append( '</ul></div></html>' )
#    f = file(html_file,'w')
#    f.write("\n".join( rval ))
#    f.write('\n')
#    f.close()

if __name__ == "__main__": __main__()
