#!/usr/bin/env perl

#==============================================================================================

# Script to output statistsics and histograms for reads and contigs/isotigs


# Outputs include:
#    Mean, N50, StdDev or reads or contig lengths,
#    Mean and Modal read or contig lengths.
#    Number of reads or contigs > 1 kb in length 
#    Summed contig length (by number of contigs, in sorted order)
#    Histogram of read or contig lengths,
#    Graph of sums of read-lengths
#    File of reads or contigs sorted by read or contig length
#    Test for mono/di-nucelotide repeats
#    Randomly selected reads or contigs


# Needs gnuplot installed to create the histograms:
#   On Fedora/Redhat linux: sudo yum install gnuplot 
#   On Ubuntu/Debian: sudo apt-get install gnuplot

#  Uses a linux pipe to call gnu-plot directly, rather than as a separate shell script.

# Original written by Sujai Kumar, 2008-09-05 University of Edinburgh
# Modified by Stephen: 29-Apr-2009:
# Last changed by Stephen: 9-Aug-2010


# Usage: fasta_summary.pl -i infile.fasta  -o process_reads -t read OR contig OR isotig (to use 'read' or 'contig' or 'isotig' in the output table & graphs. Isotig is for 'runAssembly -cdna ...' output file '454Isotigs.fna') [-r 1 to indicate count simple nucleotide repeats] [-n number of random reads to output] [-c cutoff_length] [-l 1 to indicate output the longest read] [-f (s or t or w) for spacer, tab or wiki format table output.]

# Note: The parameters above in the [] are optional.

# eg: fasta_summary.pl -i myfile.fasta  -o process_reads -t read
# Where:
#  -i reads or contigs as input, in fasta format.
#  -o output_dir (created if it doesn't exist)
#  -t read, contig or isotig

# Gives back
# - N50
# - num of contigs > 1 kb
# - num of contigs
# - Read or Contig Histogram and graphs.
# - Summed contig length (by number of contigs, in sorted order)

#==============================================================================================


use strict;
use warnings;
use Getopt::Long;

my $infile;
my $output_dir;
my $type='read'; # Defaults to 'read' at present
my $repeats=1;
my $num_random_reads_to_output=0;
my $cutoff_length=-1; # -1 means won't check this cutoff
my $longest_read=-1; # -1 mean's don't output the sequence for the longest read.
my $doCommify=1; # Outputs statistics numbers in format: 9,999,999
my $format="t";  # "s"=spaces between columns, "t"=tabs between columns, "w"=wiki '||' and '|'.
my $bucket1=0; # For optional exact length histogram distribution as asked for by JH.

if ($#ARGV==-1) {die "
 Usage:  

    fasta_summary.pl -i infile.fasta  -o output_dir  -t ( read | contig | isotig ) [ -r 0 ] [ -n num_reads ] [ -c cutoff_length ] [ -l 1 ] [ -d 0 ] [ -f (w | t ) ] [ -bucket1 ]

 where:

    -i or -infile  infile.fasta  :  input fatsa file of raeds, contigs or isotigs, 

    -o or -output_dir output_directory : directory to put output stats and graphs into.

    -t or -type (read or contig or isotig) : for displaying the graph title, where type is 'read' or 'contig' or 'isotig'.

    -r or -repeats 0 or 1        :  1=count number of reads that contain over 70% simple mono-nucleotide and di-nucleotide repeat bases; 0=don't count.

    -n or -number num_reads      : For outputting specified number of randomly selected reads or contigs.

    -c or -cutoff cutoff_length  : Give a number of reads to do extra analysis (calculating again the number of reads and number of bases in reads above this length)

    -l or -longest 0 or 1        : 1=Output the longest read;  0= don't output the longest read

    -d or -doCommify 0 or 1      : Output numbers formatted with commas to make easier to read: 0=no commas, default=1

    -f or -format w or t         : w=wiki_format (ie. table with || and | for column dividers), t=tabs between column symbols for the wiki pages, default is spaces between columns.

    -b or -bucket1               : To also output histogram file of exact read lengths (ie. bucket size of 1)


  eg: For 454 reads:      fasta_summary.pl  -i RawReads.fna -o read_stats -t read
      For 454 isotigs:    fasta_summary.pl  -i 454Isotigs.fna -o isotig_stats -t isotig

";}

GetOptions (
	"infile=s"      => \$infile,
	"output_dir=s"  => \$output_dir,
	"type=s"        => \$type, ## type is 'read' or 'contig' or 'isotig' - for displaying the graph title
	"repeats=i"     => \$repeats, # To count simple repeats
	"number=i"      => \$num_random_reads_to_output, # For outputting specified number of random reads
	"cutoff=i"      => \$cutoff_length, # Give a number of reads to do extra analysis (calculating again the number of reads and number of bases in reads above this length)
	"longest=i"     => \$longest_read, # Output the longest read.
	"doCommify=i"   => \$doCommify, # Output numbers formatted with commas to make easier to read: 0=no commas, default=1
	"format=s"      => \$format, # "w"=wiki_format (ie. table with || and | for column dividers), "t"=tabs between column symbols for the wiki pages, default is spaces between columns.
        "bucket1"       => \$bucket1, # To also output histogram file of exact read lengths (ie. bucket size of 1)
);
if ($#ARGV>-1) {die "Unused options specified: @ARGV\n";}
if ( (! defined $infile) || ($infile eq '') ) {die "\nPlease give input fasta file, preceeded by the '-i' option\n\n";}
if ( (! defined $output_dir) || ($output_dir eq '') ) {die "Please give output_directory, preceeded by the '-o' option\n\n";}
if ( (! defined $type) || (($type ne 'contig') && ($type ne 'read') && ($type ne 'isotig')) ) {die "ERROR: On commandline: -t type must be 'contig' or 'read' or 'isotig'\n\n";}


my ($L,$M,$R, $Lh,$Mh,$Rh, $Lhnewline,$Mhnotab);
if    ($format eq 's') {($L,$M,$R, $Lh,$Mh,$Rh, $Lhnewline,$Mhnotab)=('','  ','', '','  ','', "\n",'');}
elsif ($format eq 't') {($L,$M,$R, $Lh,$Mh,$Rh, $Lhnewline,$Mhnotab)=("\t","\t",'', "","\t",'', "\n",'');} # There is correctly a tab for the $L, but not the $Lh.
elsif ($format eq 'w') {($L,$M,$R, $Lh,$Mh,$Rh, $Lhnewline,$Mhnotab)=('| ',' | ',' |', '|| ',' || ',' ||', '|| ',' || ');}
else {die "\nInvalid output format code: '$format'. Should be 's', 't' or 'w'.\n\n";}

### create output_dir if it doesn't exist
if (-d $output_dir) {
	print STDERR "  Directory '$output_dir' exists, so the existing fasta_summary.pl output files will be overwritten\n";
} else {
	mkdir $output_dir;
	print STDERR "  Directory '$output_dir' created\n";
}

my $gc_count = 0;

#--------------- Read in contigs from fasta file -------------------

open INFILE, "<$infile" or die "Failed to open file: '$infile' : $!";
open STATS, ">$output_dir/stats.txt" or die "Failed to open $output_dir/stats.txt: '' : $!";

my $header = <INFILE>;
if (! defined($header) ) {print "\n** ERROR: First line of input fasta file is undefined - so file must be empty **\n\n"; print STATS "No sequences found\n"; exit 1;}
if ($header!~/^>/) {print "\nERROR: First line of input fasta file should start with '>' BUT first line is: '$header'\n"; print STATS "No sequences found\n"; exit 1;}

my $seq = "";
my @sequences;

while (<INFILE>) {
	if (/^>/) {
		push @sequences, [length($seq), $header, $seq];
		$header = $_;
		$seq = "";
	} else {
		chomp;
		$gc_count += tr/gcGC/gcGC/;
		$seq .= $_;
	}
}
push @sequences, [length($seq), $header, $seq];
close INFILE;
if ($#sequences==-1) {print "\nERROR: There are zero sequences in the input file: $infile\n\n"; print STATS "No sequences found\n"; exit 1;}


my $all_contigs_length=0;
foreach (@sequences) {$all_contigs_length += $_->[0];}
if ($all_contigs_length==0) {print "\nERROR: Length of all contigs is zero\n\n"; exit 2;}

# Find number and number of bases in reads greater than the optional cut-off length given at command-line.
my $num_reads_above_cutoff=0;
my $num_of_bases_in_reads_above_cutoff=0;
if ($cutoff_length>0) 
  {
  foreach (@sequences)
    {
    if ($_->[0]>=$cutoff_length) {$num_of_bases_in_reads_above_cutoff+= $_->[0]; $num_reads_above_cutoff++;}
    }
  }


#--------------- Gather Plots Data, Find N50, Print sorted contig file -------------------

my $summed_contig_length = 0;
my @summed_contig_data; # <-- For graph of summed length (in number of bases) versus number of contigs.
my @summed_contig_data_contigLens; # <-- Added by SJBridgett to get graph of summed contig length versus min. contig length included (ie. X-axis is sort of inverse of above)

my $contig1k_count = 0;
my $contig1k_length = 0;

open SORTED, ">$output_dir/sorted_contigs.fa" or die $!;

# top row in stats file
#print STATS "N50\nMax contig size\nNumber of bases in contigs\nNumber of contigs\nNumber of contigs >=1kb\nNumber of contigs in N50\nNumber of bases in contigs >=1kb\nGC Content of contigs\n";

my $N50size=-1;
my $N50_contigs = 0;

my @sorted_by_contig_length = sort {$b->[0] <=> $a->[0]} @sequences;

### variables and initialization for histogram (stored in @bins)
my $max = $sorted_by_contig_length[0][0];
my $mean= $all_contigs_length/($#sequences+1); # <-- Added by Stephen Bridgett.  Note: as $# gives the highest index number, so add 1 as arrays are zero-based.

# Calculate standard deviation
my $sumsquares = 0;
foreach (@sequences) {$sumsquares += ($_->[0] - $mean) ** 2;} # <-- Taken from John's "mean_fasta_length.pl" script.
my $stddev = ( $sumsquares/($#sequences+1) ) ** 0.5;

my $min = 0;
# Aim for approximately 100 bins, so 

my $bin_size=1;
my $min_max_range=$max - $min;
# $bin_size = ($min_max_range)/(99); # (99 is 100-1) so   1000/100
if    ($min_max_range>=100000000) {$bin_size=1000000;}
elsif ($min_max_range>=10000000) {$bin_size=100000;}
elsif ($min_max_range>=1000000) {$bin_size=10000;}
elsif ($min_max_range>=100000) {$bin_size=1000;}
elsif ($min_max_range>=10000) {$bin_size=100;}
else  {$bin_size=10;}  # elsif  ($min_max_range>=1000) {$bin_size=10;}
#elsif ($min_max_range>=100) {$bin_size=1;}
#elsif ($min_max_range>=10) {}
#elsif ($min_max_range>=1) {}
# WAS: my $bin_size = ($type eq 'contig') ? 1000 : 10;

my @bins;
$#bins = int(($min_max_range)/$bin_size) + 1; # <-- Set the bins array size.
foreach (@bins) {$_ = 0};

foreach (@sorted_by_contig_length) {

	my $curr_contig_length = $_->[0];
	push @summed_contig_data_contigLens, $curr_contig_length; # <-- added by Stephen.
	
	$bins[int(($curr_contig_length + 1 - $min)/$bin_size)]++;
	
	$summed_contig_length += $curr_contig_length;
	push @summed_contig_data, $summed_contig_length;

	### sorted contigs file
	print SORTED $_->[1] . $_->[2] . "\n";

	if ($curr_contig_length >= 1000) {
		$contig1k_count++;
		$contig1k_length += $curr_contig_length;
	}
	
	$N50_contigs++ unless ($N50size>-1); # Was unless $N50_found
	
	if ($summed_contig_length > ($all_contigs_length / 2) and $N50size == -1) {
		$N50size = $curr_contig_length;
	}
}


if ($bucket1!=0) 
  {
=pod
  # This firsdt method works and agress with the second, but the lengths are in reverse order, at the @sorted_by_contig_length array was sorted with longest contig first.
  open BUCKET1, ">$output_dir/lengths_hist1.txt" or die "Failed to open file '$output_dir/lengths_hist1.txt' : $!\n";
  print BUCKET1 "Length\tFrequency\n";
  my $len=-1;
  my $count=0;
  foreach (@sorted_by_contig_length)
    {  
    if ( $len != $_->[0] ) {if ($len>-1) {print BUCKET1 "$len\t$count\n";} $len=$_->[0]; $count=0;}
    $count++;
    }
  if ($len>-1) {print BUCKET1 "$len\t$count\n";} # Print length of final length grouping.
  close BUCKET1;
=cut
  open BUCKET1, ">$output_dir/lengths_hist1_with_zeros.txt" or die "Failed to open file '$output_dir/lengths_hist1_with_zeros.txt' : $!\n";
  print BUCKET1 "Length\tFrequency\n";
  my @bucket=(); # To check the result by using array.
  foreach (@sequences)
    {
    my $len=$_->[0];
    if (defined $bucket[$len]) {$bucket[$len]++;}
    else {$bucket[$len]=1;}
    }
  for (my $i=0; $i<=$#bucket; $i++)  
# for (my $i=$#bucket; $i>=0; $i--) # <-- for reverse order
   {
   if (defined $bucket[$i]) {print BUCKET1 "$i\t$bucket[$i]\n";} 
   else {print BUCKET1 "$i\t0\n";} # Can uncomment this later if don't want zeros in the output.
   }
  close BUCKET1;
  }


my $type_plural=$type.'s';
print STATS  $Lh."Statistics for $type lengths:".$Mhnotab.$Rh."\n";
print STATS  $L."Min $type length:".$M.&commify_if($sorted_by_contig_length[$#sequences][0],$doCommify).$R."\n";
print STATS  $L."Max $type length:".$M.&commify_if($max,$doCommify).$R."\n";
printf STATS $L."Mean %s length:".$M."%.2f".$R."\n", $type,$mean; # <-- Added by Stephen Bridgett, April 2009.
printf STATS $L."Standard deviation of %s length:".$M."%.2f".$R."\n", $type,$stddev;  ## <-- Added by Stephen Bridgett, May 2009.
print STATS  $L."Median $type length:".$M.&commify_if($sorted_by_contig_length[int($#sequences/2)][0],$doCommify).$R."\n";
print STATS  $L."N50 $type length:".$M.&commify_if($N50size,$doCommify).$R."\n";

print STATS  $Lhnewline."Statistics for numbers of $type_plural:".$Mhnotab.$Rh."\n";
print STATS  $L."Number of $type_plural:".$M.&commify_if($#sequences+1,$doCommify).$R."\n";
print STATS  $L."Number of $type_plural >=1kb:".$M.&commify_if($contig1k_count,$doCommify).$R."\n";
print STATS  $L."Number of $type_plural in N50:".$M.&commify_if($N50_contigs,$doCommify).$R."\n";

print STATS  $Lhnewline."Statistics for bases in the $type_plural:".$Mhnotab.$Rh."\n";
print STATS  $L."Number of bases in all $type_plural:".$M.&commify_if($all_contigs_length,$doCommify).$R."\n";
print STATS  $L."Number of bases in $type_plural >=1kb:".$M.&commify_if($contig1k_length,$doCommify).$R."\n";
printf STATS $L."GC Content of %s:".$M."%.2f %%".$R."\n", $type_plural,(100*$gc_count/$all_contigs_length);

if ($cutoff_length>0) 
  {
  print STATS $Lhnewline."Statistics for $type_plural >= $cutoff_length bp in length:".$Mhnotab.$Rh."\n";
  print STATS $L."Number of $type_plural >= $cutoff_length bp:".$M.&commify($num_reads_above_cutoff,$doCommify).$R."\n";
  print STATS $L."\tNumber of bases in $type_plural >= $cutoff_length bp:".$M.&commify($num_of_bases_in_reads_above_cutoff,$doCommify).$R."\n";
  }

if ($repeats==1) {&countRepeats();}

print STATS "\n";


# Output random selection of reads if requested on commandline:
my $fastaLineLen=60; # <-- The line length used for 454 sffinfo output, but could use a value read from input file (but be careful not to read a short line)
if ($num_random_reads_to_output>0)
  {
  my @randlist;
  if ($num_random_reads_to_output<($#sequences+1))
    {    
    print STATS "\nSome randomly selected reads:\n\n";
    @randlist= &getListOfRandomNumbers($#sequences, $num_random_reads_to_output); # Don't use ($#sequences + 1), just ($#sequences) otherwise would be outside the array.
    }
  else # Just print all the sequences:
    {
    print STATS "\nAll ".($#sequences+1)." reads:\n\n";
    for (my $i=0;$i<=$#sequences;$i++) {push @randlist,$i;}
    }
  &print_sequences(\@randlist)
  }    


# Print the longest read:
if ($longest_read>0)
  {
  my $length_of_longest_read=-1;
  my @longest_read=();
  my $i=0;
  foreach (@sequences)
    {
    if ($_->[0]>$length_of_longest_read) {$length_of_longest_read=$_->[0]; $longest_read[0]=$i;}
    $i++;
    }
  if ($length_of_longest_read>0) {print STATS "\nLongest read:\n"}
  &print_sequences(\@longest_read);
  }


=pod
print STATS "\n$type\tSummed\nlength\tlength\n"; # <-- Added by Stephen Bridgett, but better to produce a graph instead.

my $i=0;
foreach (@summed_contig_data) {
#	print STATS $sorted_by_contig_length[$i]->[0]."\t".$summed_contig_data_contigLens[$i]."\t".$_."\t".$summed_contig_data[$i]."\n";
	print STATS $sorted_by_contig_length[$i]->[0]."\t".$_."\n";
	$i++;
}
=cut

open  SUMMED, ">$output_dir/summed_contig_lengths.dat" or die $!;
print SUMMED join "\n",@summed_contig_data;
close SUMMED;

open  HISTOGRAMBINS, ">$output_dir/histogram_bins.dat" or die $!;
my $bin_size_counter = 0;
foreach (@bins) {
	print HISTOGRAMBINS eval($bin_size_counter++ * $bin_size + $bin_size/2) . "\t$_\n";
}
close HISTOGRAMBINS;


# Graph of cumulative (summed) number of reads on y-axis, versus length of read (decending order) on x-axis
open  SUMREAD_READLEN, ">$output_dir/sum_reads_vs_read_len.dat" or die $!;
#my $read_counter= 0;
my $read_counter= $#sorted_by_contig_length+1;

foreach (@sorted_by_contig_length) {
#	$read_counter++;
	$read_counter--;
	print SUMREAD_READLEN "$_->[0]\n"; # $read_counter
}
close SUMREAD_READLEN;




my $properType=ucfirst($type); # Makes the first letter an upper case letter, ie. 'Config' or 'Read'
#if ($type eq 'contig') 
#  {
  # print the outcome of the gnu_plot as may have a write permissions error sometimes.
  my $YHistogramScaleType = ($type eq 'read') ? '' : 'log y';    # Not using log scale for reads, just for contig/isotigs.
  &plot_graph('histogram', "$output_dir/histogram_bins.dat",  "Histogram of $type lengths", "$properType length", "Number of $type_plural", '0.9', $YHistogramScaleType);
  &plot_graph('line',      "$output_dir/summed_contig_lengths.dat", "Summed $type lengths", "Number of $type_plural", "Summed $type length in bases", '0.9', '');
  &plot_graph('line',      "$output_dir/sum_reads_vs_read_len.dat", "X-axis gives the Number of $type_plural that are greater than the $properType-length given on the Y-axis", "$properType length", "Cummulative number of $type_plural",  '0.9', '');

=pod
  #  print `gnuplot_histogram.sh $output_dir/histogram_bins.dat`;

  &plot_graph("$output_dir/summed_contig_lengths.dat", 'Summed contig lengths', 'Number of contigs', 'Summed contig length in bases', '0.9', '');
  #  print `gnuplot_summedcontigs.sh $output_dir/summed_contig_lengths.dat`;

  &plot_graph('line', "$output_dir/sum_reads_vs_read_len.dat", 'X-axis gives the Number of contigs that are greater than the Contig-length given on the Y-axis', 'Contig length', 'Cummulative number of contigs',  '0.9', '');
  #  print `gnuplot_sum_contig_vs_contig_len.sh $output_dir/sum_reads_vs_read_len.dat`;
  }
elsif ($type eq 'read')
  {

  #  print `gnuplot_readshistogram_logY.sh $output_dir/histogram_bins.dat`; # There's also optionally a "...._linearY.sh"

  &plot_graph('line', "$output_dir/summed_contig_lengths.dat",'Summed read lengths', 'Number of reads', 'Summed read length in bases', '0.9', '');
  #  print `gnuplot_summedreads.sh $output_dir/summed_contig_lengths.dat`;

  &plot_graph('line', "$output_dir/sum_reads_vs_read_len.dat", 'X-axis gives the Number of reads that are greater than the Read-length given on the Y-axis', 'Read length', 'Cummulative number of reads', '0.9', '');
  #  print `gnuplot_sum_reads_vs_read_len.sh $output_dir/sum_reads_vs_read_len.dat`;  
  }
else  {die "\n** ERROR: Invalid type='$type' **\n\n";}
=cut

close SORTED;
close STATS;


# Use pipe to plot directly with gnuplot, rather than calling a separate shell script:
# http://www.vioan.ro/wp/2008/09/30/calling-gnuplot-from-perl/
# http://forums.devshed.com/perl-programming-6/plotting-with-gnuplot-within-perl-script-549682.html
# Another option is the perl module: GnuplotIF: http://lug.fh-swf.de/perl/GnuplotIF.html OR: http://lug.fh-swf.de/perl/

# PlPlot: Perl: http://search.cpan.org/~dhunt/PDL-Graphics-PLplot-0.47/plplot.pd
#               http://plplot.sourceforge.net/
# dislin: http://www.mps.mpg.de/dislin/overview.html
# MathGL: http://mathgl.sourceforge.net/index.html


sub plot_graph
{
# Plots a histogram or xy-line graph
# Parameters: GraphType (histogram/line) DataFile, Title, X-label, Y-label, Y-range
# Graphfile should end with '.png'
# The $yloglinear is 'log y' for log, or '' for linear
my ($graphtype, $datafile, $title,$xlabel,$ylabel,$yrange,$yloglinear)=@_;  # yrange for reads: 0.1, and for contigs: 0.9

my $graphstyle='';
if ($graphtype eq 'histogram') {$graphstyle="plot \"$datafile\" using 1:2 with boxes";}
elsif ($graphtype eq 'line') {$graphstyle="plot \"$datafile\" using 1 with lines";}
else {die "\n** ERROR: Invalid graphtype='$graphtype'\n\n";}
my $yloglinearscale= ($yloglinear eq '') ? '' : "set $yloglinear";

# To capture any errors that are normally sent from gnuplot to stderr, could use: open3 pipe interface: 
# http://www.clintoneast.com/articles/perl-open3-example.php
# http://hell.org.ua/Docs/oreilly/perl2/prog/ch16_03.htm
# But the following should be fine, as the stderr will display when running the script anyway.
# If needed a simpler way would be to sent the output to a file using eg: open (GNUPLOT, "|gnuplot > gnu_out.txt 2>&1")  or die .... The read the resulting file.
open (GNUPLOT, "|gnuplot") or die "\n**ERROR: Failed to open gnuplot : $!\n\n **";
print GNUPLOT <<ENDPLOT;
set terminal png
set output "$datafile.png"
set nokey
$yloglinearscale
set xlabel "$xlabel"
set ylabel "$ylabel"
set yrange [$yrange:]
set title "$title"
$graphstyle
ENDPLOT
close(GNUPLOT);
if ($? != 0) {print "\n** WARNING: GNUplot pipe returned non-zero status: '$?'\n\n";} # $? is the status returned by the last pipe close (or backtick or system operator)
if (! -e "$datafile.png") {die "\n** ERROR: Failed to create '$datafile.png'**\n\n";}

=pod
#PNG
set term png small xFFFFFF
set output "$file.png"
set size 1 ,1
set nokey
set data style line
set xlabel "frequency" font "VeraMono,10"
set title "Fast Fourier Transform" font "VeraMono,10"
set grid xtics ytics
set xtics 100
plot "$file" using 1:2 w lines 1
=cut

#WAS PREVIOUSLY AS .sh script
=pod
# The 'gnuplot_readshistogram_logY.sh' is:
set terminal png
set output "$1.png"
set log y
set xlabel "Read length"
set ylabel "Frequency"
set yrange [0.9:]
set title "Histogram of read lengths"
plot "$1" using 1:2 with boxes
=cut
}


# Was previously a separate .sh file:
=pod
#!/bin/sh
gnuplot << EOF
set terminal png
set output "$1.png"
set xlabel "Number of contigs"
set ylabel "Summed contig length in bases"
set yrange [0.9:]
set title "Summed contig lengths"
plot "$1" using 1 with lines
EOF
=cut




# Added function to count number of simple dinucleotide repeats:
sub countRepeats 
  {
  # To count the number of sequences that contain mostly repeats.
  # This would be faster if called a C program on the file.

# Common simple repeats are listed here: http://www.bioinfo.de/isb/2005/05/0041/
# Dinucleotide
#   AT/TA
#   AC/TG
#   AG/TC
#   CG/GC
# Trinucleotide	
#   AAT/TTA
#   CTA/GAT
#   ATG/TAC
#   ACT/TGA
#   CTC/GAG
#   AGG/TCC
#   CAG/GTC
#   AAG/TTC
#   ATA/TAT
#   CAA/GTT
#   AGC/TCG
#   ACA/TGT
#   ACG/TGC
#   AGA/TCT
#   ACC/TGG
#   Other
# Tetranucleotide	
#   AAAT
#   AAAC
#   CACG
#   AACA
#   AATA
#   AAGA
#   TGAA
#   AAAG
#   ACAT
#   AATG
#   AGCC
#   Other
# Pentanucleotide	
#   AAAAC
#   AATTG
#   GCTAA
#   ATAAT
#   AAAAT
#   AAACA
#   ATATA
#   TTGCC
#   Other

  # I also add mono-nucleotide repeats: - ie. just all T's, or A's, etc
  # Just consider the dinucleotide repeats for now:
  my ($ATseq,$CGseq,$ACseq,$TGseq,$AGseq,$TCseq)=(0,0,0,0,0,0);
  my ($AAseq,$TTseq,$CCseq,$GGseq)=(0,0,0,0);
  foreach (@sequences)
    {
    my $seq_len=$_->[0];
    my $seq=$_->[2]; # This copy might be slow, maybe should just stick with using the reference.
    my $mnt=0.35*$seq_len;  # Mononucleotide threshold:  HERE 0.35 also means 70%;  0.4 means 80% dinucleotide repeats, as really one base so 0.5 = 100%
    my $dnt=0.35*$seq_len;  # Dinucleotide threshold:  0.35 means 70%;  0.4 means 80% dinucleotide repeats, as two bases so 0.5 = 100%
    my ($AT,$CG,$AC,$TG,$AG,$TC)=(0,0,0,0,0,0);
    my ($AA,$TT,$CC,$GG)=(0,0,0,0);
    # See: http://www.allinterview.com/showanswers/76719.html

    # AT/TA seems most common repeat so process it first to save time:
    while ($seq=~/AT/g) {$AT++;} if ($AT>$dnt) {$ATseq++; next;} # AT is same as TA.   If has 80% AT's then won't have 80% AC etc.
    while ($seq=~/CG/g) {$CG++;}  if ($CG>$dnt) {$CGseq++; next;} # CG is same as GC.
    # AC,TG
    while ($seq=~/AC/g) {$AC++;} if ($AC>$dnt) {$ACseq++; next;}
    while ($seq=~/TG/g) {$TG++;} if ($TG>$dnt) {$TGseq++; next;}
    # AG/TC
    while ($seq=~/AG/g) {$AG++;} if ($AG>$dnt) {$AGseq++; next;}
    while ($seq=~/TC/g) {$TC++;} if ($TC>$dnt) {$TCseq++; next;}

    # Added my simple mononucleotde repeat count:
    while ($seq=~/AA/g) {$AA++;} if ($AA>$mnt) {$AAseq++; next;}
    while ($seq=~/TT/g) {$TT++;} if ($TT>$mnt) {$TTseq++; next;}
    while ($seq=~/CC/g) {$CC++;} if ($CC>$mnt) {$CCseq++; next;}
    while ($seq=~/GG/g) {$GG++;} if ($GG>$mnt) {$GGseq++; next;}
    }

  my $num_seq=($#sequences+1);
  my $total_din_repeats_seq= $ACseq+$TGseq+$ATseq+$AGseq+$TCseq+$CGseq;
  my $percent_din_repeats=100*$total_din_repeats_seq/$num_seq;
  print STATS "\nSimple Dinucleotide repeats:\n";
  printf STATS "\tNumber of %s with over 70%% dinucleotode repeats:\t%.2f %% (%d %s)\n", $type_plural, $percent_din_repeats, $total_din_repeats_seq, $type_plural;
  printf STATS "\tAT:\t%.2f %% (%d %s)\n", (100*$ATseq/$num_seq),$ATseq,$type_plural;
  printf STATS "\tCG:\t%.2f %% (%d %s)\n", (100*$CGseq/$num_seq),$CGseq,$type_plural;
  printf STATS "\tAC:\t%.2f %% (%d %s)\n", (100*$ACseq/$num_seq),$ACseq,$type_plural;
  printf STATS "\tTG:\t%.2f %% (%d %s)\n", (100*$TGseq/$num_seq),$TGseq,$type_plural;
  printf STATS "\tAG:\t%.2f %% (%d %s)\n", (100*$AGseq/$num_seq),$AGseq,$type_plural;
  printf STATS "\tTC:\t%.2f %% (%d %s)\n", (100*$TCseq/$num_seq),$TCseq,$type_plural;

  my $total_mono_repeats_seq= $AAseq+$TTseq+$CCseq+$GGseq;
  my $percent_mono_repeats=100*$total_mono_repeats_seq/$num_seq;
  print STATS "\nSimple mononucleotide repeats:\n";
  printf STATS "\tNumber of %s with over 50%% mononucleotode repeats:\t%.2f %% (%d %s)\n", $type_plural, $percent_mono_repeats, $total_mono_repeats_seq, $type_plural;
  printf STATS "\tAA:\t%.2f %% (%d %s)\n", (100*$AAseq/$num_seq),$AAseq,$type_plural;
  printf STATS "\tTT:\t%.2f %% (%d %s)\n", (100*$TTseq/$num_seq),$TTseq,$type_plural;
  printf STATS "\tCC:\t%.2f %% (%d %s)\n", (100*$CCseq/$num_seq),$CCseq,$type_plural;
  printf STATS "\tGG:\t%.2f %% (%d %s)\n", (100*$GGseq/$num_seq),$GGseq,$type_plural;

  return $percent_din_repeats+$percent_mono_repeats;
  }


sub commify_if
  {
  # If doCommify is >0 then converts output to commas.   
  # Formats '1234567890.01' with commas as "1,234,567,890.01
  # Based on: http://www.perlmonks.org/?node_id=110137
  my ($number,$doCommify)=@_;
  if ($doCommify > 0) {$number =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;}
  return $number;
  }


#--------------- Produce ordered list of random numbers -------------------
# This is copied from: my_random_contigs.pl

sub getListOfRandomNumbers
{
# Use: @list=getListOfRandomNumbers(200,20); to return sorted list of 20 numbers in range from 0 to 200 inclusive.
my %list2=();
my $i=0;
my $MaxNumber=$_[0];
my $NumToPick=$_[1];
while ($i<$NumToPick)
  {
  my $intRand = int(rand($MaxNumber+1)); # For Zero-based perl-arrays. The +1 means this generates random integers between 0 and $MaxNumber. (See: http://perldoc.perl.org/functions/rand.html )
  if ($intRand>$MaxNumber) {$intRand=$MaxNumber} # Just to be extra sure that don't exceed $MaxNumber.
  if ( !exists($list2{$intRand}) ) {$list2{$intRand}=1; $i++;}
  }
#foreach my $key(keys %list2) {print "$key ";} 
# Sort the list of numbers:
#my @SortedList2 = sort { $a <=> $b } keys(%list2);
#return @SortedList2;
return (sort { $a <=> $b } keys(%list2));
#print "Sorted list of ".$NumToPick." random numbers:\n";
#foreach my $num(@SortedList2) {print "$num\n";}
#print "\n\n";
}


sub print_sequences
  {
  # Print the sequences wrapping sequences using index array, at line length of '$fastaLineLen' characters:
  # Uses the global '@sequences' array.
  my $sequence_indexes_list=$_[0]; # This is an array reference, not the array itself.
  foreach my $num(@{$sequence_indexes_list})
    {
#print "$num (max=$#sequences)\n";
    print STATS $sequences[$num]->[1]; # Prints the header, no "\n" needed after it.
    my $pos=0;
    my $seqlen=$sequences[$num]->[0];
    while ($pos<$seqlen)
      {
      print STATS substr($sequences[$num]->[2],$pos,$fastaLineLen)."\n";
      $pos+=$fastaLineLen;
      }
    print STATS "\n";
    }
  }


=pod
# Some test runs for mono-nucleotides and dinucelotides:
>FUOMOGO01AQV42DUMMYA length=339 xy=0189_0676 region=1 run=R_2009_04_23_17_54_06_
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>FUOMOGO01AQV42DUMMYB length=339 xy=0189_0676 region=1 run=R_2009_04_23_17_54_06_
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>FUOMOGO01AQV42DUMMYC length=339 xy=0189_0676 region=1 run=R_2009_04_23_17_54_06_
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
>FUOMOGO01AQV42DUMMYD length=339 xy=0189_0676 region=1 run=R_2009_04_23_17_54_06_
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG

>FUOMOGO01AQV42 length=339 xy=0189_0676 region=1 run=R_2009_04_23_17_54_06_
TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG
TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG
TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG
TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG
TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG
TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT
>FUOMOGO01AUK0D length=214 xy=0231_0843 region=1 run=R_2009_04_23_17_54_06_
ACACACACACACACACACACACACACACACACACACACACACACACACACACACACACAC
ACACACACACACACACACACACACACACACACACACACACACACACACACACACACACAC
ACACACACACACACACACACACACACACACACACACACACACACACACACACACACACAC
ACACACACACACACACACACACGACGACGACGAC
>FUOMOGO01AUB7C length=64 xy=0228_1718 region=1 run=R_2009_04_23_17_54_06_
ATATATATATATATATATATATATATATATATATATATATATATATATATAGTACGTACG
TACG
>FUOMOGO01AU00B length=213 xy=0236_1097 region=1 run=R_2009_04_23_17_54_06_
ACACACACACACACACACACACACACACACACACACACACACACACACACACACACACAC
ACACACACACACACACACACACACACACACACACACACACACACACACACACACACACAC
ACACACACACACACACACACACACACACACACACACACACACACACACACACACACACAC
ACACACACACACACACACACGACGACGACGACG
>FUOMOGO01ATYRT length=169 xy=0224_0695 region=1 run=R_2009_04_23_17_54_06_
TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG
TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG
TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT
>FUOMOGO01ARMLN length=400 xy=0197_2201 region=1 run=R_2009_04_23_17_54_06_
TATATATATATATATATATATATATATATATATATATATATATATATATATATATATATA
TATAGTAGTAGTAGTATATATATATATATATATATATATATATATATATATATATATATA
TATATATATATATATATATATATATATATATATATATATATATATATATATATATATATA
TATATATATATATATATATATATATATATATATATATATATATATATATATATATATATA
TATATATATATATATATATATATATATATATATATATATATATATATATATATATATATA
TATATATATATATATATATATATATATATATATATATATATATATATATATATATATATA
TATATATATATATATATATATATATATATATATATATATA
>FUOMOGO01AVGRX length=44 xy=0241_1051 region=1 run=R_2009_04_23_17_54_06_
TATATATATATATATATATATATATATATATATATATATATATA
>FUOMOGO01ASZ6K length=315 xy=0213_0922 region=1 run=R_2009_04_23_17_54_06_
TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG
TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG
TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG
TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG
TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG
TGTGTGTGTGTGTGT
>FUOMOGO01ARSZF length=65 xy=0199_2281 region=1 run=R_2009_04_23_17_54_06_
TATATATATATATATATATATATATATATATATATATATATATATATATATAGTACGTAC
GTACG
>FUOMOGO01AYV8U length=49 xy=0280_1324 region=1 run=R_2009_04_23_17_54_06_
ATATATATATATATATATATATATATATATATATATATATATATATATA
>FUOMOGO01AYV9X length=40 xy=0280_1363 region=1 run=R_2009_04_23_17_54_06_
TATATATATATATATATATATATATATATATATATATATA
>FUOMOGO01AUX4M length=40 xy=0235_1460 region=1 run=R_2009_04_23_17_54_06_
TATATATATATATATATATATATATATATATATATATATA
>FUOMOGO01AWOTU length=54 xy=0255_0800 region=1 run=R_2009_04_23_17_54_06_
ATATATATATATATATATATATATATATATATATATATATATATATATATAGTA
>FUOMOGO01A11TC length=66 xy=0316_1054 region=1 run=R_2009_04_23_17_54_06_
ATATATATATATATATATATATATATATATATATATATATATATATATATATAGTACGTA
CGTACG
>FUOMOGO01ASRJP length=401 xy=0210_2019 region=1 run=R_2009_04_23_17_54_06_
TATATATATATATATATATATATATATATATATATATATATATATATATATATAGTATAT
AGTAGTAGTAGTATATATATATATATATATATATATATATATATATATATATATATATAT
ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT
ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT
ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT
ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT
ATATATATATATATATATATATATATATATATATATATATA
>FUOMOGO01AU1ZH length=67 xy=0236_2363 region=1 run=R_2009_04_23_17_54_06_
TATATATATATATATATATATATATATATATATATATATATATATATATATATAGTACGT
ACGTACG
=cut
