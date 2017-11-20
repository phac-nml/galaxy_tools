#!/usr/bin/env perl

use strict;
use warnings;

use Bio::DB::Sam;
use Bio::SeqIO;
use Getopt::Long;

my ($bam, $ref, $output);
my $overlap_cutoff = 65;
my $identity_cutoff = 75;
GetOptions (
    "bam=s"             => \$bam,
    "overlap_cutoff=s"  => \$overlap_cutoff,
    "identity_cutoff=s" => \$identity_cutoff,
    "ref=s"             => \$ref,
    "output=s"          => \$output
);

if ( not defined $bam or not defined $ref or not defined $output or $overlap_cutoff < 0 or $overlap_cutoff > 100 or $identity_cutoff < 0 or $identity_cutoff > 100 ) {
    die("Usage: ./bam2mappingstats.pl --bam example.bam --overlap_cutoff 65 --identity_cutoff 75 --ref ref.fasta --output output.csv\nError in command line arguments\n");
}

my $seqIO = Bio::SeqIO->new(-file=>$ref, -format=>"fasta");
my $seq   = $seqIO->next_seq();
my $refid = $seq->display_id;

my $sam = Bio::DB::Sam->new(
    -bam       => "$bam",
    -fasta     => "$ref",
    -autoindex => 1,
);

my @stats;
my @nt_stats;

my $caller = sub {
    my ($seqid,$pos,$p) = @_;
    my $refbase = $sam->segment($seqid,$pos,$pos)->dna;
    $stats[$pos]{'wt'} = uc($refbase);

    for my $pileup (@$p) {
        my $b = $pileup->alignment;
        next if $pileup->is_refskip;

        my $identity = $b->query->length / $b->l_qseq * 100;

        if ($identity >= $overlap_cutoff) {
            my $matches = 0;
            if ( $b->has_tag('NM') ) {
                $matches = $b->query->length - $b->get_tag_values('NM');
            }
            else {
                for (my $i=0; $i < scalar(@{$b->cigar_array}); $i++) {
                    if ( @{@{$b->cigar_array}[$i]}[0] eq 'M' || @{@{$b->cigar_array}[$i]}[0] eq '=' ) {
                        $matches += @{@{$b->cigar_array}[$i]}[1];
                    }
                }
            }

            my $percent_identity = $matches / $b->query->length * 100;

            if ($percent_identity >= $identity_cutoff) {

                if ( not defined $stats[$pos] ) {
                    $stats[$pos] = {};
                }

                if ( $pileup->indel > 0 ) {
                    $stats[$pos]{'ins'} += 1;
                }
                elsif ( $pileup->indel < 0 ) {
                    for ( my $i=$pos; $i<$pos-$pileup->indel; $i++ ) {
                       $stats[$i]{'del'} += 1;
                    }
                }
                else {
                    my $qbase = substr($b->qseq,$pileup->qpos,1);
                    if ( $qbase =~ /[nN]/ ) {
                        $stats[$pos]{'ns'} += 1;
                    }
                    elsif ( lc($refbase) ne lc($qbase) ) {
                        if ( ( $qbase =~ /[aA]/ and $refbase =~ /[gG]/ ) or
                            ( $qbase =~ /[cC]/ and $refbase =~ /[tT]/ ) ) {
                            $stats[$pos]{'tsit'} += 1;
                            $nt_stats[$pos]{uc($qbase)} += 1;
                        }
                        else {
                            $stats[$pos]{'tver'} += 1;
                            $nt_stats[$pos]{uc($qbase)} += 1;
                        }
                    }
                    else {
                        $stats[$pos]{'m'} += 1;
                        $nt_stats[$pos]{uc($refbase)} += 1;
                    }
                }
            }
        }
    }
};

$sam->max_pileup_cnt(1000000);
$sam->pileup("$refid",$caller);

open my $fh, ">", $output;

print $fh "#position,wildtype,matches,mismatches,ins,del,a,c,g,t,snp_rate,ins_rate,del_rate\n";
for ( my $pos=1; $pos<scalar(@stats); $pos++ ) {
    my $m          = $stats[$pos]{'m'}    || 0;
    my $tsit       = $stats[$pos]{'tsit'} || 0;
    my $tver       = $stats[$pos]{'tver'} || 0;
    my $ns         = $stats[$pos]{'ns'}   || 0;
    my $ins        = $stats[$pos]{'ins'}  || 0;
    my $del        = $stats[$pos]{'del'}  || 0;
    my $quality    = $stats[$pos]{'qual'} || 0;
    my $mismatch_rate = 0;
    my $ins_rate = 0;
    my $del_rate = 0;
    if (($tsit+$tver+$m) > 0) {
    	$mismatch_rate = ($tsit+$tver)/($tsit+$tver+$m);
        $ins_rate = $ins/($tsit+$tver+$m);
        $del_rate = $del/($tsit+$tver+$m);
    }
    my $wildtype   = $stats[$pos]{'wt'}   || substr($seq->seq,$pos-1,1);
    my $a          = $nt_stats[$pos]{'A'} || 0;
    my $c          = $nt_stats[$pos]{'C'} || 0;
    my $g          = $nt_stats[$pos]{'G'} || 0;
    my $t          = $nt_stats[$pos]{'T'} || 0;

    printf $fh ("%i,%s,%i,%i,%i,%i,%i,%i,%i,%i,%.10f,%.10f,%.10f\n", $pos, $wildtype, $m, $tsit+$tver, $ins, $del, $a, $c, $g, $t, $mismatch_rate, $ins_rate, $del_rate);
}

close $fh;
