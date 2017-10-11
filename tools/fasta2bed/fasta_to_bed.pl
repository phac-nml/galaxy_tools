#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;
use Bio::Seq;
use Getopt::Long;

my ($fasta_file, $out);

GetOptions(
	"i|input=s" => \$fasta_file, #contigs fasta file
	"o|output=s" => \$out
);

my $seqio_object= Bio::SeqIO->new (-format =>'fasta', -file=>$fasta_file);
open(my $out_fh, ">", $out) || die "Could write to file '$out'\n";

while (my $seq_object = $seqio_object->next_seq()) {
	my $seq_id = $seq_object->display_id();
	my $length = $seq_object->length();
	print $out_fh $seq_id . "\t" . "1" . "\t" . $length . "\t" . $seq_id . "\n";

}


