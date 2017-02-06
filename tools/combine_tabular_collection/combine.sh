#!/bin/bash

output=$1
shift

i=1

for var in "$@"
do 
	if [[ -s $var ]] ; then
		( head -q -n 1 $var ) > $output
		break
	fi
	i=$[i+1]
done

if [ $i -le "$#" ]
	then
		( tail -q -n +2  $@ )>> $output
else
	exit 5
fi


