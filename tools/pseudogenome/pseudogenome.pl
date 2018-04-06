#!/usr/bin/env perl
use strict;
use warnings;
use autodie qw(:all);
use Bio::SeqIO;
use Readonly;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
Getopt::Long::Configure('bundling');

=head1 NAME

nml_pseudogenome.pl - To create a single pseudo genome out of multiple contigs provided in a single fasta file. Contig are combined in order of appearances in file

=head1 SYNOPSIS

nml_pseudogenome.pl -i F<file_name.fna> -n 100 -c X -o F<filename.fna>

=head1 OPTIONS

=over

=item B<-i>, B<--input>

Multiple fasta file

=item B<-n>, B<--number>

Number of filler base pairs to be added, default : 10

=item B<-c>, B<--chars>

Character to be used as the 'glue' between contigs, default : 'N'

=item B<--id>

Name of fasta file to be used default: pseudogenome

=item B<-o>, B<--output>

Output file name, default : Same as input

=item B<-s>, B<--stitch>

Add the stitch pattern between contigs only

=item B<-h>, B<--help>

Print this help

=item EXAMPLE

nml_pseudogenome.pl -i multiple_fasta.fna -n 100 -c X -o pseudo.fna

nml_pseudogenome.pl -i another_multiple.fna

=back

=head1 DESCRIPTION

To create a single pseudo genome out of multiple contigs provided in a single fasta file. Contig are combined in order of appearances in file.

=cut

# Nonsub perlcode

Readonly my $DEFAULT_NUM_CHAR => 10;
Readonly my $stitch_pattern => 'NNNNNCACACACTTAATTAATTAAGTGTGTGNNNNN';
Readonly my $DEFAULT_CHAR => 'N';
my ( $input,$id, $number, $char, $output,$stitch, $help );

GetOptions(
    'i|input=s'    => \$input,
    'n|number=s' => \$number,
    'c|char=s'   => \$char,
    'o|output=s'   => \$output,
    'h|help'       => \$help,
    's|stitch' => \$stitch,
    'id=s' => \$id
);
($input,$id,$number,$char,$output) =check_inputs( $input, $number, $char,$output, $help,$stitch );



my $in = Bio::SeqIO->new(-file=>$input,-format=>'fasta');

my $sequence;

#go thru every sequence and append to main sequence
while (my $seq = $in->next_seq()) {
    if ($stitch) {
        $sequence .= $seq->seq . $stitch_pattern;
    }
    else {
        $sequence .= $seq->seq . ($char x $number );        
    }

}

my $main = Bio::Seq->new(-display_id=>$id,-seq=>$sequence);

my $out = Bio::SeqIO->new(-file => ">$output" ,-format=>'fasta');
$out->write_seq($main);

exit;

=begin HTML

=head2 check_inputs

     Title   : check_inputs
     Usage   : check_inputs($fasta,$num,$filler,$out_to,$usage);
     Function: check arguments provided by the user to see if they are usable and more or less correct
     Returns : Return 1 if all is correct,otherwise 0
     Args    : $query: Query that we are looking for in the database.  Could be accession number or locus_tag
               $db: Name of database we are looking for using the query provided
               $format: Ensure that format was given by user and is valid format
               $usage: If true, return usage description
     Throws  : none

=cut

sub check_inputs {
    my ( $fasta, $num, $filler, $out_to, $usage,$use_stitch ) = @_;

    if ( $help || !( $fasta || $num || $filler || $out_to ) ) {
	pod2usage();
        exit;
    }

    if ( !($fasta) || !( -e $fasta ) ) {
        print STDERR "Error: Input file not given or does not exist\n";
	pod2usage();
        exit;
    }

    if ($use_stitch) {
        print "Using stitch pattern\n";
        
    }
    else {
        if ( !$num ) {
            $num = $DEFAULT_NUM_CHAR;
            print STDERR "Number of character not given, using $num\n";
        }
        elsif ( !( $num =~ /^\d+$/xms ) ) {
            print STDERR "Error: Number of character was not a number\n";
	pod2usage();
            exit;
        }
        
        if ( !$filler ) {
            $filler = $DEFAULT_CHAR;
            print STDERR "No filler character given, using 'N'\n";
        }
        
    }

    if ( !($out_to) ) {
        $out_to = fileparse($fasta) . ".pseudogenome";
        print
          "Output file was not given. Result will be written to '$out_to'\n";
    }
    if ( ! $id) {
        $id = 'pseudogenome';
    }

    return ( $fasta,$id, $num, $filler, $out_to );
}

=end HTML

=head1 SEE ALSO

No related files.

=head1 AUTHOR

Philip Mabon, <philip.mabon@canada.ca>

=head1 BUGS

None reported.

=head1 COPYRIGHT & LICENSE 

Copyright (C) 2018 by Public Health Agency of Canada

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 DEVELOPER PAGE

No developer documentation.

=cut
