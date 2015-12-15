#/bin/bash

name=$1
shift

spolpred $@

sed -i s/^.*\t/$name/ output.txt

exit 0
