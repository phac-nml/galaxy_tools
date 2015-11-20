#!/usr/bin/perl

use strict;
use warnings;
use File::Copy;
use File::Basename;
use Getopt::Long;
use Pod::Usage;

my ($html_file, $html_path, @args, %information, $folder);

Getopt::Long::Configure('bundling');
GetOptions(
	'i|input=s'	=> \@args,
	'h|html=s'	=> \$html_file,
	'p|path=s'	=> \$html_path
	);

pod2usage(1) unless @args && $html_file && $html_path;

#At this point, the output directory does not exist yet. So we have to make it
mkdir $html_path or die "Could not make directory $!";

#Now make a folder for all our files
my $data_folder = $html_path."/Bundled_Collection";
mkdir $data_folder or die "Could not make directory $!";

#Go through each list item
foreach my $entry (@args)
{
	#Get key and value. Remove any spaces
	my ($info, $file) = split /=/, $entry;
	my ($name, $ext) = split /,/, $info;
	$name=~s/ /_/g;
	my $full_name = $name.".".$ext;

	#We store this for later to make the html file
	$information{$name}{$ext}=1;

	#copy each file to its directory

	my $file_path = $data_folder."/".$full_name;

	copy($file,$file_path) or die "Could not copy $file to $file_path: $!";
}

#Write out the html file
open my $out, ">", $html_file or die "Could not open html file: $!";

my $num_keys = scalar(keys %information);
my $num_vals = scalar(@args);

printf $out "<!DOCTYPE html>
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

<h2 id=\"M0\">Bundle Collection Summary</h2><br><br>

Number of keys: $num_keys<br>
Number of values: $num_vals<br><br> 

<table border=\"1\"><tr><th>File name</th><th>File type</th></tr>";

foreach my $key (sort(keys %information))
{

	foreach my $val (keys %{$information{$key}} )
	{
		printf $out "<tr><td>$key</td><td>$val</td></tr>";
	}
}

printf $out "</table></body></html>";

close $out;

__END__

=head1 name
	
		bundle_collection.pl - Downloads a collection from Galaxy

=head1 SYNOPSIS

	bundle_collection.pl -h html_file -p output_path -o "key=value"

=back
=cut