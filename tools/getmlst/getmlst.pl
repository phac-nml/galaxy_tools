#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;
use File::Copy;

my ($mlst_db, $mlst_defs, $species) = @ARGV;

$species =~ s/__pd__/#/ig;

my $command = "getmlst.py --species '$species'";

my $rv = system($command);

if ($rv == 0)
{
	#need to find output files in the dir
	my $cur_dir = getcwd();

	foreach my $file (<$cur_dir/*>)
	{
		if ($file =~ /\.fasta$/)
		{
			move($file, $mlst_db);
		}
		elsif ($file =~/\.txt$/)
		{
			move($file, $mlst_defs);
		}
	}
}
$rv = $rv >> 8;
exit $rv;
