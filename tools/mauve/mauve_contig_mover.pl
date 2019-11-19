#!/usr/bin/perl

use strict;
#use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Copy;
use File::Basename;
#use Archive::Zip;


my ($output, $reference_gbk, $reference_dat, $draft_fasta, $draft_dat, $alignment_file, $fasta_file, $html_file, $help, $best_alignment, 
	@files, @sorted_files, $num_of_alignments, $mauve_cmd, $out, $best_alignment_file, $best_fasta_file);

Getopt::Long::Configure('bundling');
GetOptions(
	'r|reference=s' => \$reference_dat,
	'd|draft=s'	=> \$draft_dat,
	'o|output=s'	=> \$output,
	'a|alignment=s'	=> \$alignment_file,
	'f|fasta=s'	=> \$fasta_file, 
	'l|html=s'	=> \$html_file,
	'h|help'	=> \$help
);
pod2usage(1) if $help;

#Format the fasta file. Some fastas were not working, so we're going to format all incoming fastas
my $draft_temp = "temporary.fasta";

open my $in, '<', $draft_dat or die "Could not open draft file: $?";
open my $out, '>', $draft_temp or die "Could not open file for writing: $?";

my $first_header = <$in>;
chomp $first_header;
print $out "$first_header\n";
while (my $line = <$in>)
{
   chomp $line;

   if ($line =~ />/)
   {
      print $out "\n$line\n";
   }
   else
   {
      print $out $line;
   }
}



#progressiveMauve checks the file extension of inputs, and did not like .dat files passed in.

#First get the file format for the extention.
my $format = get_format($reference_dat);
die "Input reference file isn't a properly formatted fasta or genbank file!\n" if $format eq "na";

#So here we create symbolic links to the .dat files using the proper file extensions.
$reference_gbk = "reference.".$format;
$draft_fasta = "draft.fasta";

system("ln -s $reference_dat $reference_gbk");
system("ln -s $draft_temp $draft_fasta");

#First, run mauve with the given inputs
$mauve_cmd = "java -Djava.awt.headless=true -Xmx500m -cp \$path2jar org.gel.mauve.contigs.ContigOrderer -output ".$output." -ref ".$reference_gbk." -draft ".$draft_fasta;

$out = system("$mauve_cmd");


#Get all of the alignmentX folders and pick the latest one
opendir(DIR, $output) or die "Can't opendir $output";
@files = readdir(DIR);
@sorted_files = sort @files;
$best_alignment = $sorted_files[@sorted_files -1];
$num_of_alignments = $best_alignment;
$num_of_alignments  =~ s/[^\d.]//g;


#Now let's give galaxy the right outputs. Is there a better way to do this?
#First extract the names from the paths	
my ($a_name, $a_path, $a_suffix) = fileparse($best_alignment);
my ($d_name, $d_path, $d_suffix) = fileparse($draft_fasta, ".fasta");


#Now we want some files (alignment and the final fasta file) to show up 
#in the history. So we copy these files to send them to galaxy.	
$best_alignment_file = $output."/".$best_alignment."/".$a_name;
$best_fasta_file = $output."/".$best_alignment."/".$d_name.".fasta";

#Now copy them to the galaxy locations
copy($best_alignment_file, $alignment_file) or die "$best_alignment_file Copying alignment failed: $!";
copy($best_fasta_file, $fasta_file) or die "Copying fasta file failed: $!";


#Let us write the html file!
open my $html_out, ">", $html_file;
printf $html_out "<!DOCTYPE html>
<html>
<style type=\"text/css\">

body {
	font-family: sans-serif;
	color: #000;
	}

table {
	margin-left: 3em;
	text-align: center;
	}
th {
	text-align:center;
	background-color: #000080;
	color: #FFF;
	padding: 0.4em;
	}
td {
	font-family: monospace;
	text-align: left;
	background-color: #EEE;
	color: #000;
	padding: 0.4em;
	}
h2 {
	color: #800000;
	padding-bottom: 0;
	margin-bottom: 0;
	clear: left;
	}
</style></head>

  <body>


    <h2 id=\"M0\">Mauve Output Summary</h2><br>
 <ul><li>Number of alignments performed: $num_of_alignments (last alignment is usually the best)</li><li>The fasta and alignment files from the last alignment are shown in the history</li><li>To download the complete Mauve output, download the contents of this file</li><li>The contig orders of each alignment are shown below: </li></ul><br>";


my %summary;
my $contig_count;

#generate data html
foreach my $folder (@sorted_files)
{
	my $start = "Ordered Contigs";
        my $stop = "Contigs with conflicting ordering information";
        my $header = "type      label   contig  strand  left_end        right_end";
        my $in_section = 0;
        $contig_count = 0;
        my @alignment_info;
        next if ($folder =~ m/^\./);
	my $file_path = $output."/".$folder."/".$d_name."_contigs.tab";
	open my $curr_file, "<", $file_path;

	#Go through the lines
	while(<$curr_file>)
	{
		#find correct section. Note flip flop operators didn't work here for me 
		if(/$start/) {$in_section=1;}
		elsif(/$stop/) {$in_section=0;}

		next if(/$start/ || /$stop/ || /$header/);

		#Store the line content
		if ($in_section and /\S/)
		{
			$contig_count++;
			my @tmp = split('\t', $_);
			my @columns = ($tmp[1], $tmp[3]);
			push(@alignment_info, [@columns]);
		}
	}
	close $curr_file;
	$summary{$folder} = [@alignment_info];
}


#print out the headers

printf $html_out "<table border=\"1\"><tbody><tr><th>Alignment</th>\n";

for my $a (sort keys %summary)
{
	my $tmp = $a;
	$tmp =~ s/[^\d.]//g;
        printf $html_out "<th colspan=\"2\">$tmp</th>\n";
}

printf $html_out "</tr>";

#print out the data
for my $i (0 .. $contig_count-1)
{
        printf $html_out "<tr><td></td>\n";
        for my $alignment (sort keys %summary)
        {
                printf $html_out "<td>". $summary{$alignment}[$i][0]."</td>\n<td>".$summary{$alignment}[$i][1]. "</td>\n";
        }
        printf $html_out "</tr>\n";
}

printf $html_out "</tbody></table></body></html>";

#close all the things
closedir(DIR);
close $html_out;

exit($out);


sub get_format
{
	my $file = shift;
	my $format;

	open my $in, '<', $file or die "Could not open file for reading. $!";

	my $line = <$in>;

	if ($line =~ /LOCUS/)
	{
		$format = "gbk";
	}
	elsif ($line =~ /^>/)
	{
		$format = "fasta";
	}
	else
	{
		$format = "na";
	}
	close $in;

	return $format;
}

__END__

=head1 NAME

	mauve_contig_mover.pl - A wrapper for galaxy to run Mauve Contig Mover

=head1 SYNOPSIS

	mauve_contig_mover.pl -r <reference> -d <draft> -o <output> -a <alignment file output> -f <fasta file output> -l <html file output> -h <help>

=head1 OPTIONS

=over 8

=item B<-r> B<--reference>

The input reference strain in either a fasta or genbank format

=item B<-d> B<--draft>

The input draft genome in fasta format

=item B<-o> B<--output>

The output folder created by Mauve

=item B<-a> B<--alignment>

The best output alignment produced by Mauve

=item B<-f> B<--fasta>

The best output fasta file produced by Mauve

=item B<-l> B<--html>

The html file containing all of the output files produced by Mauve

=item B<-h> B<--help>

Print a help message and exits

=back

=head1 DESCRIPTION

B<mauve_contig_mover> is a galaxy wrapper for Mauve Contig Mover. This script runs the command line version of the Mauve Contig Mover

=cut
