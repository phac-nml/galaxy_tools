#!/usr/bin/env perl


use strict;
use warnings;
use Cwd;
use File::Copy;

#The first 3 arguments should be in this format:
#  bam_output scores_output pileup_output ...


my ($bam_results, $scores, $pileup, $job_type, $txt_results, $genes_results, $fullgenes_results, $name, $databases);
my ($allele_results,$allele_type);


$bam_results = $ARGV[0];
shift(@ARGV);
$scores = $ARGV[0];
shift(@ARGV);
$pileup = $ARGV[0];
shift(@ARGV);

#Now that we shifted the first 4 arguments, get the rest depending on the job type. There are three:
#I pass a letter to tell us which job type was selected.
$job_type = $ARGV[0];
shift(@ARGV);

#If m, mlst only:  we only have one mlst output file
if($job_type eq "m")
{
   $txt_results = $ARGV[0];
   shift(@ARGV);
   $allele_results = $ARGV[0];
   shift(@ARGV);
   $allele_type = $ARGV[0];
   shift(@ARGV);
}
#If g, genedb only: we have two outputs: genes and fullgenes. We also get the name of the database
elsif($job_type eq "g")
{
   $genes_results = $ARGV[0];
   shift(@ARGV);
   $fullgenes_results = $ARGV[0];
   shift(@ARGV);
   $databases = $ARGV[0];
   shift (@ARGV);
}
#If b, both mlst and genedb: We will have three output files and the database name.
else
{
   $txt_results = $ARGV[0];
   shift(@ARGV);
   $genes_results = $ARGV[0];
   shift(@ARGV);
   $fullgenes_results = $ARGV[0];
   shift(@ARGV);
   $databases = $ARGV[0];
   shift(@ARGV);
}

#After we get the output files/database name, now we get the name of the reads
#This allows SRST2 to give meaningful output instead of just printing 'dataset_xxx' as the sample name
my $filename = $ARGV[0];
shift(@ARGV);

#This index offset is used to determine where the 'genedb' option is located.
#If only a single-end input, the genedb will be 3 positions into the arguments:
#  -input_se sample.fastq --genedb database.fasta
my $index_offset = 3;

print @ARGV;
#change the file extensions of the input reads so srst can use them

#Usually the file name looks like this: sample_S1_L001_R1_001.fastq
#If we use this file name, it confuses srst2 a lot. So we just extract
#the sample name to use as input file name.
my @file_name = split /_/, $filename;
$name = "temp_file_name";

my ($for_read, $rev_read, $sing_read, $database);
if ($ARGV[0] eq "--input_pe")
{
        #Increment index offset if we paired-end inputs
        $index_offset++;

	$for_read = $name."_1.dat";
	$rev_read = $name."_2.dat";

	symlink($ARGV[1], $for_read);
	symlink($ARGV[2], $rev_read);

	$ARGV[1] = $for_read;
	$ARGV[2] = $rev_read;
}
else
{
        $sing_read = $name.".dat";

        symlink($ARGV[1], $sing_read);

        $ARGV[1] = $sing_read;

}
#If we are running a job to include genedb, use the database name for input file name
if ($job_type eq 'g' | $job_type eq 'b')
{
  my @db_names = split /,/, $databases;
  my $num_db = @db_names;
  my %names_hash = ();
  # loop through dbs to replace spaces with _ and check for duplicates
  for (my $i = 0; $i < $num_db; $i++){
    $db_names[$i]=~s/ /_/g;
    if( exists($names_hash{$db_names[$i]}) ) {
      print STDERR  "More than one database with the same name";
      exit(1);
    }else{
      $names_hash{$db_names[$i]}=$db_names[$i];
    }
  }


  foreach my $db_name (@db_names){
    (my $base = $db_name) =~ s/\.[^.]+$//;
    $database = $base.".dat";

    symlink($ARGV[$index_offset], $database);
    $ARGV[$index_offset] = $database;

    $index_offset++;
  }
}

for (my $i  =0; $i< @ARGV; $i++){
  if (index($ARGV[$i], "maxins") != -1){
    my ($maxins, $minins);
    my @b2args = split(' ', $ARGV[$i]);
    for (my $j = 0; $j < @b2args; $j++){
      if (index($b2args[$j], "maxins") != -1){
        $maxins = $b2args[$j+1];
      }
      if (index($b2args[$j], "minins") != -1){
        $minins = $b2args[$j+1];
      }
    }
    if ($maxins - $minins < 0){
      print STDERR  "--minins cannot be greater than --maxins";
      exit(1);
    }
  }
}



my $command = "srst2 @ARGV";

my $exit_code = system($command);

my $cur_dir = getcwd();
# make arrays for using multiple custom databases (creates multiple output files - need to be concatenated)
my (@genefiles, @bamfiles, @pileupfiles, @fullgenefiles, @scoresfiles);

# go through files in the output directory to move/concatenate them as required.

foreach my $file (<$cur_dir/*>)
{
    print $file, "\n";
    #Will cause problems if any files have 'mlst' in them
  if ($file =~ /mlst/)
  {
    move($file, $txt_results);
  }
  elsif ($file =~ /\.bam$/)
  {
    push @bamfiles, $file;
  }
  elsif ($file =~ /\.scores$/)
  {
    push @scoresfiles, $file;
  }
  elsif ($file =~ /\.pileup$/)
  {
    push @pileupfiles, $file;
  }
  elsif ($file =~ /__fullgenes__/)
  {
    push @fullgenefiles, $file;
  }
  elsif ($file =~ /__genes__/)
  {
    push @genefiles, $file;
  }
  elsif ($file =~ /all_consensus_alleles.fasta$/ && $allele_type eq 'all') {
    move($file,$allele_results);
  }
  elsif ($file =~ /new_consensus_alleles.fasta$/ && $allele_type eq 'new'){
    move($file,$allele_results);
  }
}


my ($cmd, $temp_head, $temp_full, $cat_header, $final_bam, @headers );

# create new concatenated bam file with all bam files.
if (@bamfiles > 1){
  my $counter = 0;
  $cat_header = "cat_header";
  while ($counter < @bamfiles) {
    $headers[$counter] = "bam_header".$counter;
    # make a header file for each bam results file
    my $cmd = "samtools view -H $bamfiles[$counter] > $headers[$counter]";
    system($cmd);
    if ($counter >= 1){
      # only keep the @hd and @pg from first file because the final concatenated file can only have one of each (doesn't matter location)
      $temp_head="cut_head".$counter;
      # cut off first row and last row of each file (the @HD and @PG)
      $cmd = "tail -n +2 $headers[$counter] | head -n -1 > $temp_head";
      system($cmd);
      unlink $headers[$counter];
      # replace the old header with the new cut header in the array
      $headers[$counter] = $temp_head;
    }
    $counter++;
  }
  # combine all header files
  $cmd = "cat ".join(" ",@headers)." > $cat_header";
  system($cmd);

  $final_bam = "final_bam";
  # concatenate all the bam files *must include the concatenated header file created above
  $cmd = "samtools cat -h $cat_header -o $final_bam ".join(" ",@bamfiles)." ";
  system($cmd);

  # sort the bam file so it can be indexed
  $cmd = "samtools sort $final_bam 'last_bam'";
  system($cmd);

  # move bam file to where Galaxy expects it.
  $cmd = "mv 'last_bam.bam' $bam_results";
  system($cmd);
} else {
  # only one bam file, don't need to concatenate
  move($bamfiles[0], $bam_results);
}

# concatenate all pileup files
if (@pileupfiles > 1){
  $cmd = "cat ".join(" ",@pileupfiles)." > $pileup";
  system($cmd);
} else {
  move($pileupfiles[0], $pileup);
}

# perform find-replace to restore original user-specified file names
foreach my $gene (@genefiles){
  my $data = read_file($gene);
  $data =~ s/temp_file_name/$file_name[0]/g;
  write_file($gene, $data);
}

foreach my $gene (@fullgenefiles){
  my $data = read_file($gene);
  $data =~ s/temp_file_name/$file_name[0]/g;
  write_file($gene, $data);
}

# concatenate gene files with a space separating each file
if (@genefiles > 1){
  my $join = join(" <(echo) ", @genefiles);
  my @args = ("bash", "-c", "cat $join > $genes_results");
  system(@args);
} else {
  # only one gene results file
  move($genefiles[0], $genes_results);
}

# concatenate full gene results files but only keep header in first file.
if (@fullgenefiles >1){
  for (my $i= 1; $i < @fullgenefiles; $i++){
    # go through all files but the first one to remove headers
    # create a temp file to save the file after header is removed
    $temp_full = "temp_full".$i;
    $cmd = "tail -n +2 $fullgenefiles[$i] > $temp_full";
    system($cmd);
    unlink $fullgenefiles[$i];
    $fullgenefiles[$i] = $temp_full;
  }
  $cmd = "cat ". join(" ",@fullgenefiles)." > $fullgenes_results";
  system($cmd);
} else{
  # only one full gene results file
  move($fullgenefiles[0], $fullgenes_results);
}

# concatenate full gene results files but only keep header in first file.
if (@scoresfiles >1){
  for (my $i= 1; $i < @scoresfiles; $i++){
    # go through all files but the first one to remove headers
    # create a temp file to save the file after header is removed
    $temp_full = "temp_full".$i;
    $cmd = "tail -n +2 $scoresfiles[$i] > $temp_full";
    system($cmd);
    unlink $scoresfiles[$i];
    $scoresfiles[$i] = $temp_full;
  }
  $cmd = "cat ". join(" ",@scoresfiles)." > $scores";
  system($cmd);
} else{
  # only one scores file
  move($scoresfiles[0], $scores);
}

# cleanup srst2 output and temp files
foreach my $file (@fullgenefiles){
  unlink $file;
}
foreach my $file (@genefiles){
  unlink $file;
}
foreach my $file (@pileupfiles){
  unlink $file;
}
foreach my $file (@bamfiles){
  unlink $file;
}
foreach my $file (@headers){
  unlink $file;
}
foreach my $file (@scoresfiles){
  unlink $file;
}
unlink $temp_head;
unlink $temp_full;
unlink $cat_header;
unlink $final_bam;

#get rid of symlinks
if ($for_read)
{
	unlink $for_read;
}
if ($rev_read)
{
	unlink $rev_read;
}
if ($sing_read)
{
         unlink $sing_read;
}
if ($database)
{
         unlink $database;
}
$exit_code = $exit_code >> 8;
exit $exit_code;

sub read_file {
    my ($filename) = @_;

    open my $in, '<:encoding(UTF-8)', $filename or die "Could not open '$filename' for reading $!";
    local $/ = undef;
    my $all = <$in>;
    close $in;

    return $all;
}

sub write_file {
    my ($filename, $content) = @_;

    open my $out, '>:encoding(UTF-8)', $filename or die "Could not open '$filename' for writing $!";;
    print $out $content;
    close $out;

    return;
}
