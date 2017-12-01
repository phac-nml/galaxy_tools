#!/usr/bin/env perl
package nml_fastqc_stats;
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use autodie;
use File::Basename;
use Bio::SeqIO;
use File::Temp qw/ tempdir /;
__PACKAGE__->run unless caller;



sub run {
    #grabbing arguments either from command line or from another module/script
    my ( $rawdatas, $fastqs,$num_bps,$sample,$out) = get_parameters(@_);

    print_header($out,scalar @$fastqs);
    my @results;
    foreach my $rawdata ( @$rawdatas) {
        my %results = %{ parse_FastQC($rawdata) } ;
        push @results,\%results;
    }


    my %results;
    
    
    if ( scalar @$rawdatas ==1  ) {
        %results = % { $results[0]};
        
    }
    else {
        #need to combine the info into a SINGLE results file

        my ($r1,$r2) = ($results[0],$results[1]);

        
        foreach my $key( keys %$r1) {
            if ( exists $r2->{$key}) {
                my ($value1,$value2) = ($r1->{$key},$r2->{$key});

                #check to see if data that could/should be different
                my %diff_ignore = ('Filename'=>1);
                my %combine = ('duplicate_lvl'=>1,'dist_length'=>1,'overre'=>1,'Total Sequences'=>1,'Sequence length'=>1,'%GC'=>1);
                if (exists $combine{$key} ) {
                    if ( $key eq 'duplicate_lvl') {
                        $results{'Duplicate level R1'} = $value1;
                        $results{'Duplicate level R2'} = $value2;
                    }
                    elsif ( $key eq 'dist_length') {
                        my %dist=%{$value1};
                        foreach my $key( keys %$value2) {
                            #adding the number of reads to either existing range or new
                            if ( exists $dist{$key}) {
                                $dist{$key}{'count'}+=$value2->{$key}{'count'};
                                
                            }
                            else {
                                $dist{$key} = $value2->{$key};
                            }
                        }

                        $results{$key} = \%dist;
                        
                    }
                    elsif ( $key eq 'Total Sequences') {
                        $results{$key}= $value1 + $value2;
                    }
                    elsif ( $key eq 'overre') {
                        $results{$key}=$value1 + $value2;
                    }
                    elsif ( $key eq '%GC') {
                        $results{$key}=($value1 + $value2)/2;
                    }
                    elsif ( $key eq 'Sequence length') {
                        if ( $value1 eq $value2) {
                            $results{$key}= $value1;
                        }
                        else {
                            #figure out the range by splitting all values into array then sort
                            my @number = split /-/,$value1;
                            map { push @number,$_ }  split'-',$value2;
                            @number = sort {$a <=> $b } @number;
                            $results{$key} = $number[0] . '-' . $number[$#number];
                        }
                    }
                    
                }
                #both the same we are good
                elsif ( $value1 eq $value2) {
                    $results{$key}=$value1;
                }
                elsif ( ! exists $diff_ignore{$key}) {
                    die "Have different value between fastqc files for '$key'\n";
                }

                 
                
            }
            else {
                die "Cannot find key : '$key' in reverse fastqc file\n";
            }
        }
        
        
    }

    #add sample name given by user
    $results{'filename'} = $sample;
    
    #perform coverage, total base pair and SE or PE
    if ( $fastqs) {
        my $result = determine_stats($fastqs,$num_bps);
        if ( $result) {
            foreach ( keys %$result) {
                if ( exists $results{$_}) {
                    die "Have two functions with key: '$_'\n";
                }
                else {
                    $results{$_}= $result->{$_};
                    }
            }
        }
    }

    print_csv($out,\%results);

    
    
    return 1;
}

sub determine_stats {
    my ($fastqs,$num_bps) = @_;
    
    
    my %results;
    my $number = scalar @$fastqs;
    my $status;
    
    if ( $number ==1) {
        $status='SE';
    }
    elsif ( $number ==2) {
        $status='PE';
    }
    else {
        $status='N/A';
    }

    $results{'pe_status'} = $status;

    #determine total base pair

    my $all_fastq = join (' ' , map { "\"$_\"" } @$fastqs);

    #one liner from : http://onetipperday.blogspot.ca/2012/05/simple-way-to-get-reads-length.html with some modification by Philip Mabon
    my $total=`cat $all_fastq | perl -ne '\$s=<>;<>;<>;chomp(\$s);print length(\$s)."\n";' | sort | uniq -c | perl -e 'while(my \$line=<>){chomp \$line; \$line =~s/^\\s+//;(\$l,\$n)=split/\\s+/,\$line; \$t+= \$l*\$n;} print "\$t\n"'`;
    chomp $total;
    
    $results{'total_bp'} = $total;
    
    if ($total) {

        $results{'coverage'} = sprintf "%.2f", ($total/$num_bps);
       
	#reference name stuff
	#my $name = basename($reference);
        #$name =~ s/\.fasta//;
        
        $results{'reference'} = $num_bps;
        
        
    }
    
    return \%results;
    
}


sub print_csv {
    my ($out_fh,$results) = @_;
    my %results = %{$results};

    my @line;
    
    #Name
    my $name = $results{'filename'} || 'N/A';
    $name =~ s/.fastq//;
    push @line,$name;

    #Indicate if PE,SE or multiple
    my $status= $results{'pe_status'} || 'N/A';
    push @line,$status;

    #encoding of fastq file
    push @line, $results{'Encoding'} || 'N/A';
    
    #number of reads found
    my $reads = $results{'Total Sequences'} || 'N/A';

    push @line,$reads || 'N/A';

    #number of Total Base Pairs
    push @line, $results{'total_bp'} || 'N/A';

    #sequence read length range
    push @line, $results{'Sequence length'} || 'N/A';
    
    #most abundant read length

    my ($most_abundant) = sort {$results{'dist_length'}{$b}{'count'} <=> $results{'dist_length'}{$a}{'count'} } keys %{$results{'dist_length'}};
    my ($most_count) = $results{'dist_length'}{$most_abundant}{'count'} || 'N/A';
    
    push @line, $most_abundant;
    push @line, $most_count;

    #coverage against reference
    push @line, $results{'coverage'} || 'N/A';
    push @line, $results{'reference'} || 'N/A';

    
    #duplicate level
    if ( $results{'pe_status'} eq 'SE') {
        push @line, $results{'duplicate_lvl'} || 'N/A';
    }
    elsif (  $results{'pe_status'} eq 'PE') {
        push @line, $results{'Duplicate level R1'} || 'N/A';
        push @line, $results{'Duplicate level R2'} || 'N/A';
    }



    #determine percentage of reads that are over-represented sequences
    push @line, exists $results{'overre'} ? $results{'overre'} : 'N/A';
    
    print $out_fh join("\t",@line) . "\n";
    
    return;
}
    
sub print_header {
    my ($out_fh,$fastqs) = @_;

    my @headers;
    if ( $fastqs==2) {
        @headers = ('Name','SE/PE','Encoding' , '# of Reads', 'Total # Base Pairs', 'Sequence length range','Most abundant read length','# of reads for abundant','Estimated Coverage','Reference length','Duplicate % R1','Duplicate % R2','# of Overrepresented sequences');  
    }
    else {
        @headers = ('Name','SE/PE','Encoding' , '# of Reads', 'Total # Base Pairs', 'Sequence length range','Most abundant read length','# of reads for abundant','Estimated Coverage','Reference length','Duplicate %','# of Overrepresented sequences');  
    } 

    print $out_fh join("\t",@headers) . "\n";
    

    return;
}
    


sub parse_FastQC {
    my ($file) = @_;

    my %results;
    
    my %todo = (
        'Per base sequence quality' => \&per_base_quality,
        'Per sequence quality scores' => \&per_seq_quality,
        'Per base sequence content' => \&per_seq_content,
        'Per base GC content' => \&per_base_gc_content,
        'Per sequence GC content' => \&per_seq_gc_content,
        'Per base N content' => \&per_base_n_content,
        'Sequence Length Distribution' => \&seq_length_dis,
        'Sequence Duplication Levels' => \&seq_dupli_lvl,
        'Overrepresented sequences' => \&overrepresented_seq,
        'Kmer Content' => \&kmer_content,
        'Basic Statistics' => \&basic_stats
   
    );
    
    
    open my $in , '<', $file;

    #get header
    my $version = <$in>;

    #create functional code
    my $next_sec = next_section($in);
    
    my %sections;
    
    while ( my ($sec_name,$fh_lines) = $next_sec->()) {
        if (! ($sec_name || $fh_lines)) {
            last;
        }
        

        #see if we are at a beginning of new section
        if ($sec_name =~ /^>>(.*)\s+(warn|fail|pass)$/) {
            my ($name,$status) = ($1,$2);
            $sections{$name}= $status;
            
            if ( exists $todo{$name}) {
                my $result =$todo{$name}->($fh_lines);
            
                if ( $result) {
                    #combine results together
                    foreach ( keys %$result) {
                        if ( exists $results{$_}) {
                            die "Have two functions with key: '$_'\n";
                        }
                        else {
                            $results{$_}= $result->{$_};
                        }
                    }
                }
            }
        }
    }

    
    return \%results;
}


sub basic_stats {
    my ($lines) = @_;
    my $header = <$lines>;

    my %stats;
    
    while (<$lines>) {
        chomp;
        my @data = split/\t/;
        my $value = pop @data;
        if ( $value eq '>>END_MODULE') {
            last;
        }
        my $key = join(' ',@data);
        $stats{$key}=$value;
    }
    
    return \%stats;
}



sub per_base_quality {
	my ($lines) = @_;
	my %results;

	return \%results;
}

sub per_seq_quality {
	my ($lines) = @_;
	my %results;

	return \%results;
}

sub per_seq_content {
	my ($lines) = @_;
	my %results;

	return \%results;
}
sub per_base_gc_content {
	my ($lines) = @_;
	my %results;

	return \%results;
}

sub per_seq_gc_content {
	my ($lines) = @_;
	my %results;

	return \%results;
}

sub per_base_n_content {
	my ($lines) = @_;
	my %results;

	return \%results;
}

sub seq_length_dis {
    my ($lines) = @_;
    my %results;
    
    my $header = <$lines>;
    my %lengths;
    while (<$lines>) {
        chomp;
        if ( $_ =~ /^>>/) {
            next;
        }
        
        my ($key,$value) = split/\s+/;
        if ( $key =~ /(\d+)-(\d+)/) {
            $lengths{$1}{'count'} =$value;
            $lengths{$1}{'key'} =$key;
            
        }
        else {
            $lengths{$key}{'count'} =$value;
                $lengths{$key}{'key'} =$key;
        }
        
    }
        
    
    return {'dist_length' => \%lengths};
}
sub seq_dupli_lvl {
    my ($lines) = @_;

    my $perc_dupl = <$lines>;
    my $perc;
    if ( $perc_dupl =~ /^\#Total Duplicate Percentage\s+(.*)$/ ) {
        $perc = $1;
    }
    elsif ($perc_dupl =~ /^\#Total Deduplicated Percentage\s+(.*)$/ ) {
        $perc = $1;
        #version 0.11.2 and above, it indicates as Deduplicate insetad of Duplicated
        #we want to know % of Duplicate sequences instead
        $perc = $perc - 100.00;
    }

    my $header = <$lines>;
    #ignoring the results of the graph for now

    $perc = sprintf "%.2f",$perc;

    return {'duplicate_lvl' => $perc};
}

sub overrepresented_seq {
    my ($lines) = @_;
    
    my %results;
    
    my $header = <$lines>;
    my %seqs;
    my $total;
    
    while ( <$lines>) {
        chomp;
        if ( $_=~ />>END_MODULE/) {
            last;
        }
        my @data = split/\s+/;
        
        $seqs{$data[0]} =$data[2];
        $total+=$data[2];
    
    }

    if ( ! $total) {
        $results{'overre'}=0;
        
    }
    else {
        $results{'overre'} = sprintf "%.2f",$total;        
    }

    

    return \%results;
}

sub kmer_content {
	my ($lines) = @_;
	my %results;

	return \%results;
}



sub next_section {
    my ($in)=@_;
    my $lines;
    
    return sub {
        local $/ = ">>END_MODULE\n";    
        $lines= <$in>;
        my ($name,$sec);
        {
            if ($lines) {
                local $/ = "\n";
                open $sec ,'<',\$lines;
                $name = <$sec>;
                chomp $name;
            }
            
        }
        
        
        return ($name,$sec);
    }
}


sub get_parameters {
    my ($out,$sample,@rawdatas,@fastq,$ref,$rawdataSE,$rawdataPE_R1,$rawdataPE_R2,$fastqSE,$fastqPE_R1,$fastqPE_R2,$galaxy,$num_bps);
    #determine if our input are as sub arguments or getopt::long
    if ( @_ && $_[0] eq __PACKAGE__ ) {
        # Get command line options
        GetOptions(
            'o|out=s'   => \$out,
            'fastq_se=s' => \$fastqSE,
            'fastq_pe_1=s' => \$fastqPE_R1,
            'fastq_pe_2=s' => \$fastqPE_R2,
            'g_rawdata_se=s' => \$rawdataSE,
            'g_rawdata_pe_1=s' => \$rawdataPE_R1,
            'g_rawdata_pe_2=s' => \$rawdataPE_R2,
            'sample=s' => \$sample,
            'ref=s' => \$ref,
	    'num_bps=s' => \$num_bps
        );
    }
    else {
        die "NYI\n";
        #( $file, $out ) = @_;
    }

    if ( !$sample) {
        $sample = "Unknown sample";
    }
    
    if ( $fastqSE) {
        if ( ! (-e $fastqSE) ) {
            print "ERROR: Was given a fastq SE file but could not find it: '$fastqSE'\n";
            pod2usage( -verbose => 1 );
        }
        else {
            push @fastq,$fastqSE;
        }
        if ( !$rawdataSE || !( -e $rawdataSE)) {
            print "ERROR: Was not given or could not rawdata file: '$rawdataSE'\n";
            pod2usage( -verbose => 1 );
        }
        else {
            push @rawdatas,$rawdataSE;
        }
        
    }
    elsif (  !$fastqPE_R1 || ! (-e $fastqPE_R1)  || !$fastqPE_R2 || ! (-e $fastqPE_R2) ) {
        print "ERROR: Was given a fastqPE R1 or R2 file but could not find it: '$fastqPE_R1' , '$fastqPE_R2'\n";
        pod2usage( -verbose => 1 );

    }
    else {
        push @fastq,$fastqPE_R1;
        push @fastq,$fastqPE_R2;

        
        for my $rawdata ( ($rawdataPE_R1,$rawdataPE_R2)) {
            if (!$rawdata || !(-e $rawdata)) {
                print "ERROR: Was not given or could not rawdata file: '$rawdata'\n";
                pod2usage( -verbose => 1 );                
            }
            else {
                push @rawdatas,$rawdata;
            }
        }        
    }

    if ( $ref && !( -e $ref ) ) {
        print "ERROR: Was given a reference file but could not find it: '$ref'\n";
        pod2usage( -verbose => 1 );
    }

    if ($ref && $num_bps)
    {
	print "ERROR: Was given both a reference file and number of base pairs. One or the other please.";
	pod2usage( -verbose => 1 );
    }

    if ($ref)
    {
	$num_bps = 0;
        my $in = Bio::SeqIO->new(-format=>'fasta' , -file => $ref);
        while ( my $seq = $in->next_seq()) {
            $num_bps += $seq->length();
        }

	if ($num_bps == 0)
	{
	    print "ERROR: number of base pairs read from reference file is 0. Please check validity of reference file.";
	    pod2usage( -verbose=> 1 );
	}
    }

    $out = _set_out_fh($out);
    
    return (\@rawdatas,\@fastq,$num_bps,$sample,$out);
}



sub _set_out_fh {
    my ($output) = @_;
    my $out_fh;
    
    if ( defined $output && ref $output && ref $output eq 'GLOB' ) {
        $out_fh = $output;
    }
    elsif ( defined $output ) {
        open( $out_fh, '>', $output );
    }
    else {
        $out_fh = \*STDOUT;
    }

    return $out_fh;
}


    
1;

=head1 NAME

nml_fastqc_stats.pl.pl - (Galaxy only script) Parse fastqc runs into a csv summary report


=head1 SYNOPSIS

     nml_fastqc_stats.pl.pl -z fastqc.txt
     nml_fastqc_stats.pl.pl -z fastqc.txt -o results.csv --ref NC_33333.fa
     nml_fastqc_stats.pl.pl -z fastqc.txt -o results.csv --ref NC_33333.fa --fastq ya_R1.fastq --fatq ya_R2.fastq

=head1 OPTIONS

=over


=item B<-z>

FastQC txt rawdata file


=item B<-o> B<--out>

Output csv file

=item  B<--fastq>

Path to fastq file(s) . Used for providing total base pairs and coverage information


=item  B<--ref>

Path to reference file. Needed for providing coverage information

=back

=head1 DESCRIPTION




=cut
