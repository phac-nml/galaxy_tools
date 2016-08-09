#!/usr/bin/env perl
## A wrapper script to call spades.py and collect its output
use strict;
use warnings;
use File::Temp qw/ tempfile tempdir /;
use File::Copy;
use Getopt::Long;

# Parse arguments
my ($out_contigs_file,
    $out_paths_file,
    $out_log_file,
    $new_name,
    @sysargs) = @ARGV;


my $output_dir = 'output_dir';

# Create log handle
open my $log, '>', $out_log_file or die "Cannot write to $out_log_file: $?\n";

# Run program
runSpades(@sysargs);
collectOutput($new_name);
print $log "Done\n";
close $log;
exit 0;

# Run spades
sub runSpades {
    my $cmd = join(" ", @_) . " -o $output_dir";
    my $return_code = system($cmd);
    if ($return_code) {
	print $log "Failed with code $return_code\nCommand $cmd\nMessage: $?\n";
	die "Failed with code $return_code\nCommand $cmd\nMessage: $?\n";
    }
    return 0;
}

# Collect output
sub collectOutput{
    my ($new_name) = @_;
    
    # To do: check that the files are there
    # Collects output
    if ( not -e "$output_dir/transcripts.fasta") {
        die "Could not find transcripts.fasta file\n";
    }
    if ( not -e "$output_dir/transcripts.paths") {
        die "Could not find transcripts.paths file\n";
    }

    #if a new name is given for the contigs, change them before moving them
    if ( $new_name ne 'NODE') {
        renameContigs($new_name);
    }
    else {
        move "$output_dir/transcripts.fasta", $out_contigs_file;
        move "$output_dir/transcripts.paths", $out_paths_file;        
    }

    

    open LOG, '<', "$output_dir/spades.log" 
	or die "Cannot open log file $output_dir/spades.log: $?";
    print $log $_ while (<LOG>);
    return 0;
}

#Change name in contig and fastg file
sub renameContigs{
    my ($name) = @_;

    open my $in, '<',"$output_dir/transcripts.fasta" or die $!;
    open my $out,'>', $out_contigs_file;

    while ( my $line = <$in>) {
        #remove the NODE_ so we can rebuilt the display_id with our contig name with the contig number.
        #also move the remainder of the length
        if ( $line =~ />NODE_(\d+)_(.+)/) {
            $line = ">$name" . "_$1 $2\n";
        }
        print $out $line;
    }
    close $in;
    close $out;
    

    open $in, '<',"$output_dir/transcripts.paths" or die $!;
    open $out,'>', $out_paths_file;

    while ( my $line = <$in>) {
        #remove the NODE_ so we can rebuilt the display_id with our contig name with the contig number.
        #also move the remainder of the length
        if ( $line =~ />NODE_(\d+)_(.+)/) {
            $line = ">$name" . "_$1 $2\n";
        }
        print $out $line;
    }
    close $in;
    close $out;

}
