#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use Getopt::Long;

#quick and dirty script to combine a list of assembly stats tab files into a simple csv file where each row is one strain


my ($files,$output) = prepare_inputs();


my @strains = sort { $a cmp $b } keys %{$files};


#get first file so we can determine the header
my $first = shift @strains;
my $top_header;
my $second_header;



open my $out,'>',$output;
process($first,$files->{$first},$out,1);


foreach my $name( @strains) {
    process($name,$files->{$name},$out);
}



close $out;


exit;

sub process {
    my ($name,$file,$out,$header) = @_;

    my @header = ("Strain");
    my @values = ($name);
    
    open my $in,'<',$file;
    while ( <$in>) {
        chomp;

        if (length $_ ==0) {
            next;
        }

        #if we hit this section, we are done reading this file since the rest we do not care about
        if ( $_ =~ /Simple Din.*repeats/) {
            last;
        }
        
        
        my ($key,$value) = split /:/;

        #trim out the tabs
        $key =~ s/\t//g;
        $value =~ s/\t//g;

        if ( $value) {
            push @header,$key;
            push @values,$value;
        }

    }


    close $in;

    #check to see if we are printing out the header
    if ( $header) {
        print $out join ("\t",@header) . "\n";
    }
    print $out join ("\t",@values) . "\n";
    
    return;
}


sub prepare_inputs {

    my ($output,%files);
    


    if (!GetOptions('stats=s' => \%files,
                    'output=s' => \$output
                )){
        
        die "Invalid options given\n";
    }
    
    
    if ( scalar keys %files == 0){
        die "No files given\n";
    }

    return (\%files,$output);
}
