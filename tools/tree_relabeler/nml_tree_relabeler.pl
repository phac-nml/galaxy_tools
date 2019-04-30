#!/usr/bin/env perl

use strict;
use warnings;
use Bio::TreeIO;
use Bio::Tree::Tree;
use IO::String;
use Getopt::Long;
use Pod::Usage;

my ($treefile, $tabfile, $delim, $outfile, $template, $help, $replace, $tabline, $man );

GetOptions(
    'i|treefile=s'          => \$treefile,
    't|tabfile=s'           => \$tabfile,
    'o|outfile=s'           => \$outfile,
    'd|delim:s'             => \$delim,
    'p|print-template=s'    => \$template,
    'r|replace'             => \$replace,
    'h|help'                => \$help,
    'm|man'                 => \$man
);

if ($help){
    pod2usage(-verbose => 99,
              -sections => "SYNOPSIS|OPTIONS AND ARGUMENTS|DESCRIPTION|DIAGNOSTICS");
} elsif ($man){
    pod2usage(-verbose => 2);
} elsif (!$treefile) {
    pod2usage(-msg => "**Tree file is required.**\n",
              -exitval => 2,
              -verbose => 1,
              -output  => \*STDERR);
} elsif ( !(($template) || ($tabfile && $outfile)) ){
    pod2usage(-msg => "**Either select a template file, or a tab file and outfile**\n",
              -exitval => 2,
              -verbose => 1,
              -output  => \*STDERR);
} elsif ( $treefile && $template) {
    generate_template();
} elsif ( $treefile && $tabfile && $outfile ){
    relabel_tree();
}

sub generate_template
{
    check_treefile();
    print_template();
}

sub relabel_tree
{
    my %new_labels;
    check_delimiter();
    %new_labels = %{load_tabfile()};
    check_treefile();
    write_treefile(\%new_labels);
}

sub print_template
{
    open my $tempout, '>', $template or die "Could not open file: $!";
    print $tempout "#label\n";

    my $treein = Bio::TreeIO->new(
        -format => "newick",
        -file => "$treefile"
        );

    while(my $t = $treein->next_tree)
    {
        my @nodes = $t->get_nodes;

        if (scalar(@nodes) <= 1)
        {
            print STDERR "Tree is not in newick format.\n";
            exit(2);
        }
        else
        {
            foreach my $node ($t->get_leaf_nodes)
            {
                print $tempout $node->id,"\n";
            }
        }
    }
    close ($tempout);
}

sub check_delimiter
{
    if (!($delim)){
        $delim = ' ';
    } else{
        my $delimlen = length $delim;
        # delim length less than 1 indicates empty string.
        if ($delimlen < 1){
            $delim = ' ';
        }elsif ($delim =~ /[\(\)\;\,\:]/){
            print STDERR "Delimiters cannot be Newick reserved characters '(),;:'.\n";
            exit(1);
        }
    }
}

sub load_tabfile
{
    my %new_labels;
    if (!(-e $tabfile)){
        # exit if error in tab file
        print STDERR "Error opening tab file.\n";
        exit(1);
    }

    open (my $tabin, '<', $tabfile);

    # append the > to the front of the outfile string
    # $outfile = '>'.$outfile;

    # go through tab file to add new labels to a hash file
    while ($tabline = <$tabin>){
        # skip the first row if it starts with a #
        next if $tabline =~ s/^#//;
        $tabline =~ s/\r//g;
        chomp $tabline;

        if ($tabline =~ /[\(\)\;\,\:]/){
            print STDERR "New labels cannot contain Newick reserved characters '(),;:'.\n";
            exit(1);
        }

        my @splits = split("\t", $tabline);

        # Check that the tab file has more than one column 
        my $num_cols = scalar @splits;
        if ($num_cols <= 1){
            # exit if one column or less; no new info to add to tree/error with tab layout.
            print STDERR "Tab file does not contain new labels.\n";
            exit(1);
        }
        # set the hash label to the first value in a row, is the original tip label.
        my $label = $splits[0];
        # If user chose find and replace instead, get rid of the first value.
        shift @splits if ($replace);
        # join all values from @split into one string with delim separating them
        my $new_info = join($delim, @splits);
        # add the new info to the hash
        $new_labels{$label} = $new_info;
    }

    close ($tabin);

    return \%new_labels;
}

sub check_treefile
{
    if (!(-e $treefile)){
        # exit if error in tree file
        print STDERR "Error opening tree file.\n";
        exit(1);
    }


    # open tree file to check format
    if (!(-e $treefile)){
        # exit if error in tab file
        print STDERR "Error opening tree file.\n";
        exit(1);
    }
    if (-z $treefile){
        print STDERR "Tree file is empty.\n";
        exit(1);
    }
    my $linecount = 0;
    open (my $treein, '<', $treefile);
    my $line = <$treein>;
    my $nextline = <$treein>;
    if (defined $nextline){
        print STDERR "Tree is not in newick format. More than one line\n";
        exit(2);
    }
    close ($treein);

    my @chars = split("",$line);
    my $lastelem = @chars-2;
    my $bracketcount = 0;
    # look for non-spaces at end of line
    while($chars[$lastelem] eq " "){
        $lastelem--;
    }
    if ($chars[$lastelem] ne ";"){
        # newick formats end in ;
        print STDERR "Tree is not in newick format. Does not end in ; \n";
        exit(2);
    }

    foreach my $char (@chars){
        if ($char eq ")"){
            if ($bracketcount == 0){
                # There is a ) before a (
                print STDERR "Tree is not in newick format. Missing a (\n";
                exit(2);
            }else {
                $bracketcount--;
            }
        }elsif ($char eq "("){
            $bracketcount++;
        }
    }
    if ($bracketcount != 0){
        # There were not equal number of ( and )
        print STDERR "Tree is not in newick format. Brackets do not match. \n";
        exit(2);
    }
}

sub write_treefile
{
    my $temp = shift;
    my %new_labels = %{$temp};

    # open tree file as a tree
    my $treeinput = Bio::TreeIO->new(
        -file => "$treefile"
        );

    #create output tree file
    my $treeoutput = Bio::TreeIO->new(
        -format => "newick",
        -file => ">$outfile"
        );

    while(my $t = $treeinput->next_tree ){
        # check if tree is valid
        if (scalar($t->get_nodes)<=1){
            # if tree has only 1 or less nodes, then it's not in correct format
            print STDERR "Tree is not in newick format.\n";
            exit(2);
        } else{
            foreach my $label(keys %new_labels){
                my $tip = $t->find_node(-id =>$label);
                if ($tip){
                    # if found, change to the new label
                    $tip->id($new_labels{$label});
                    print $label," found and changed.\n"
                } else{
                    # if label from tab file is not found, notify the user
                    print $label," not found.\n";
                }
            }
            $treeoutput->write_tree($t);
        }
    }
}

exit;

=head1 NAME

nml_tree_relabeler.pl - Changes the tip labels on a newick formatted tree

=head1 VERSION

This documentation refers to nml_tree_relabeler.pl version 0.0.2.

=head1 SYNOPSIS

    nml_tree_relabeler.pl -i treefile [-t tabfile -o outfile (-d delim) (-r) | -p template]

=head1 OPTIONS AND ARGUMENTS

=over

=item B<-i>, B<--treefile>

The name of the tree file containing the tree to adjust the tip laels. Only accepts trees in newick format. (required)

=item B<-t>, B<--tabfile>

The name of the tab delimited file containing current tip labels and the info to be replaced/added tothe labels. The first column must contain the current tree labels. Must not contain one of the Newick reserved characters '(),:;' (required option)

=item B<-o>, B<--out>

The output file. (required option)

=item B<-d>, B<--delim>

The character to use to divide the information of the labels. Must not be one of the Newick reserved characters '(),:;' (optional)

=item B<-r>, B<--replace>

Replace the tip names. This option will replace the tree tip names with the specified labels, instead of adding them to the tip name.

=item B<-p>, B<--print-template>

The name of the output template file. Prints out a template for the tabfile.(required option)

=item B<-h>, B<--help>

To display help message

=item B<-m>, B<--man>

To display manual

=back

=head1 DESCRIPTION

=over 

nml_tree_relabeler takes a newick format tree file to modify tip labels and a tab-delimited file containing current tip labels and additional information to add to the tips in 2 or more columns. Header row of the tab delimited file must start with a '#'. An example is below:

 #label outbreak    year    location
 orgs1  outbreak1   year1   location1
 orgs2  outbreak2   year2   location2

and so on.

The information in the tab file is inserted into the tree file so the new information will appear on the tip labels.

Alternatively, nml_tree_relabeler can print out the tip names to a tab-delimited template file.

=back

=head1 DIAGNOSTICS

=over

=item B<Tree file, tab file, and output file are required>

Use the proper command line arguments (-i, -t, -o respectively) to add the filenames of the tree file, tab file, and output file.

=item B<Tree file is required>

Use the -i command line argument to add the tree file.

=itemB<Either select a template file, or a tab file and outfile>

Use the proper command line arguments to either add a template file (-p) to print a tab template, or to add a tab file and an output file (-t, -o respectively) to relabel a tree.

=item B<Label not found>

A warning that a label provided in the tab file was not found in the tree file. Relabeling continues.

=item B<Error opening tab/tree file>

An error occured while opening the tab/tree file, please check path/file.

=item B<Tree is not in newick format>

The tree file does not appear to be in newick format. Please check file and convert if necessary.

=item B<Tab file does not contain new labels>

The tab file only contains one column and therefore does not have any additional information to add to the tree. Please check the tab file.

=item B<Delimiter/tabfile cannot contain Newick reserved characters '(),;:' >

The tab file or delimiter selected contains one of the characters used in the Newick format. This will cause an error when trying to read the tree. Please modify your tab file or select a new delimiter.

=back

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over

=item use Bio::TreeIO

=item use Bio::Tree::Tree

=back

=head1 INCOMPATIBILITIES

This script only works for NEWICK formatted trees. All other tree formats are not compatible. 

=head1 AUTHOR

Jen Cabral, <jencabral@gmail.com>

=head1 BUGS AND LIMITATIONS

There are no known bugs in this module.

Please report problems to Jen Cabral, <jencabral@gmail.com>

=head1 COPYRIGHT & LICENSE 

Copyright (C) 2015 by NML

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>

=cut