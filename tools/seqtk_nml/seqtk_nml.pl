#!/usr/bin/env perl
package seqtk_nml;
use warnings;
use strict;
use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage;
use File::Copy;
__PACKAGE__->run unless caller;


my $rv;

sub get_parameters {
    my ($fastaref, $type, $coverage, $length,$log);
    my ($for,$rev,$out_for,$out_rev);
    #determine if our input are as sub arguments or getopt::long
    if ( @_ && $_[0] eq __PACKAGE__ ) {
        Getopt::Long::Configure('bundling');
        GetOptions(
            'ref=s' => \$fastaref,
            'type=s' => \$type,
            'forward=s'   => \$for,
            'reverse=s'   => \$rev,
            'out_forward=s'   => \$out_for,
            'out_reverse=s'   => \$out_rev,
            'log=s'=> \$log,
            'cov=s' => \$coverage
        );
    }

    if ( !$for || !( -e $for ) ) {
        print "ERROR: Was not given or could not find fastq file: '$for'\n";
        pod2usage( -verbose => 1 );
    }

    if ( !$out_for ) {
        print "ERROR: Was not given output file path for fastq\n";
        pod2usage( -verbose => 1 );
    }

    if ( $type eq 'paired') {
        if ( !$rev || !( -e $rev ) ) {
            print "ERROR: Was not given or could not find reverse fastq file: '$rev'\n";
            pod2usage( -verbose => 1 );
        }
        
        if ( !$out_rev ) {
            print "ERROR: Was not given output file path for reverse fastq\n";
            pod2usage( -verbose => 1 );
        }        
    }

    if ( !$coverage ) {
        print "ERROR: Was not given a coverage number\n";
        pod2usage( -verbose => 1 );
    }

    if ( $coverage <=0  ) {
        print "ERROR: Was given a coverage less than 0\n";
        pod2usage( -verbose => 1 );
    }

    
    
    if ( !$log ) {
        print "ERROR: Was not given a log file\n";
        pod2usage( -verbose => 1 );
    }

    
    return ($fastaref, $type, $coverage, $length, $log,$for,$rev,$out_for,$out_rev);
}


sub run {
    my ($fastaref, $type, $coverage, $length, $log,$for,$rev,$out_for,$out_rev) = get_parameters(@_);
    my $subsample_size;


    #open log fh here
    open my $log_fh,">" ,"$log";

    
    my @in_fastqs;
    my @out_fastqs;
    
    
    if ($type eq "single"){
	$in_fastqs[0] = $for;
        $out_fastqs[0] = $out_for;
    }elsif ($type eq 'paired' ) {
	$in_fastqs[0] = $for;
	$in_fastqs[1] = $rev;
        $out_fastqs[0] = $out_for;
        $out_fastqs[1] = $out_rev;
    }
    else {
        die "Given unknown read type of '$type'";
    }

    #get total read lengths from all fastq files given
    my $total= get_total_length(@in_fastqs);

    
    if (!($coverage)){
	$coverage = 50;
    }

    my $seq_in  = Bio::SeqIO->new(
        -format => 'fasta',
        -file   => $fastaref,
    );

    while ( my $seq = $seq_in->next_seq()) {
        $length += $seq->length();
    }


    print $log_fh "Downsampling to coverage of: $coverage\n";
    print $log_fh "Total number of Basepairs: $total\n";
    print $log_fh "Length of Reference: $length\n";
    
    my $rawcoverage = $total/$length;
    printf $log_fh "Raw Coverage: %.3f\n",$rawcoverage;

    if($rawcoverage > $coverage){
	#need to downsample
	#calculate $subsample_size
	$subsample_size = $coverage/$rawcoverage;
	printf $log_fh "subsample: %.3f",$subsample_size;

	foreach my $fastq (@in_fastqs){
            my $out = shift @out_fastqs;
            #seed always set to 42 for reproducibility 
            my $seqCommand = "seqtk sample -s42 $fastq $subsample_size > $out";
            $rv = system($seqCommand);
            #need to bit shift 8 bit because seqtk exit code for some reason are greater then standard 0-255 values that most unix application expect
            die "Error when running '$seqCommand' command" if $rv >>8; 
	}	
    } else {
	#no sampling needed, just copy the fastq's to the output
	print "Subsampling not required\n";
	foreach my $fastq (@in_fastqs){
            my $out = shift @out_fastqs;
            copy($fastq,$out) || die "Not able to copy '$fastq' to '$out' with error $!";
	}
    }



}


sub get_total_length {
    my (@files) = @_;
    my $total;
    foreach my $fastq( @files) {

        open my $in, "<",$fastq || die "Could not open file '$fastq'";
        #skip first 3 lines
        for ( 0..2) {
            my $line = <$in>;
        }

        while ( <$in>) {
            chomp;
            $total+=length($_);
            #skip first 3 lines
            for ( 0..2) {
                my $line = <$in>;
            }            
        }
        close $in;
        
        
    }
    
    return $total;
    
}

1;


=head1 NAME



seqtk_nml.pl - Down sample fastq(s) if raw coverage is above user provided level


=head1 SYNOPSIS

     seqtk_nml.pl --ref reference.fasta --forward first_R1.fastq --reverse --reverse_R2.fastq --out_forward answer_R1.fastq --out_reverse answer_R2.fastq --log log-file


=head1 OPTIONS

=over

=item

=item B<--ref>

Reference fasta file that we getting the expected length [Required]


=item B<--cov>

Coverage desired i.e 50.0


=item B<--forward>

Forward fastq read file. [Required]

=item B<--reverse>

Reverse fastq read file. Can be optional


=item B<--out_forward>

Downsampled forward fastq read file. [Required]

=item B<--out_reverse>

Downsampled reverse fastq read file. Can be optional


=item B<--log>

Log file that indicate what has happen. [Required]

=item B<--type>

Indicate to application if we are receiving one or two fastq files [Required] ['paired','single']



=back

=head1 DESCRIPTION


Downsample fastq(s) reads based on the raw coverage from reference fasta file. Needed when we have too much data to run correctly in downstream analyses tools. i.e spades , snvphyl , etc..


=cut





=back


=head1 SYNOPSIS



=head1 DESCRIPTION



=cut

