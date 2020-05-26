#!/bin/bash

#wrapper for abacas
##our arguments will come in as an array

args=("$@")

#to keep things readable assign vars
input_ref="${args[0]}"
input_query="${args[1]}"
ordered_contigs="${args[2]}"
use_bin="${args[3]}"
non="${args[4]}"
append_bin="${args[5]}"
out_bin="${args[6]}"
out_crunch="${args[7]}"
out_gaps="${args[8]}"
out_fasta="${args[9]}"
out_tab="${args[10]}"
out_unused="${args[11]}"
out_multi="${args[12]}"
out_binmulti="${args[13]}"
out_nonpseudo="${args[14]}"

##ok lets do up the optional tags.....I have four thus far
##generate ordered multifasta file

if [ "$ordered_contigs" == "yes" ]; then
	options="-m"
fi

if [ "$use_bin" == "yes" ]; then
        options="$options -b"
fi

if [ "$non" == "yes" ]; then
        options="$options -N"
fi

if [ "$append_bin" == "yes" ]; then
        options="$options -a"
fi

options="$options -o ab_out"

script_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)"

#run abacas
eval "perl $script_path/abacas.1.1.pl $options -r $input_ref -q $input_query -p nucmer"
echo "perl $script_path/abacas.1.1.pl $options -r $input_ref -q $input_query -p nucmer"
#ok now we have to remove the nucmer files to cleanup a bit...
rm nucmer.delta
rm nucmer.filtered.delta
rm nucmer.tiling

#name the datafiles properly
mv ab_out.bin "$out_bin"
if [ "$ordered_contigs" == "yes" ]; then
	mv ab_out.contigs.fas "$out_multi"
fi

if [ "$use_bin" == "yes" ]; then
	mv ab_out.contigsInbin.fas "$out_binmulti"
fi

mv ab_out.crunch "$out_crunch"
mv ab_out.fasta "$out_fasta"
mv ab_out.gaps "$out_gaps"

if [ "$non" == "yes" ]; then
	mv ab_out.NoNs.fasta "$out_nonpseudo"
fi

mv ab_out.tab "$out_tab"
mv unused_contigs.out "$out_unused"

#end script for now
