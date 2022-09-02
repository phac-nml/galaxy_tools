#!/bin/bash

if grep -q ">NODE" $1; then
    sed -r "s/>NODE(_[0-9]+)_(.*)/>$1\1 \2/g" $1 >$2
fi
if grep -q ">contig" $1; then
    sed -r "s/>contig([0-9]+) (.*)/>$1\1 \2/g" $1 >$2
fi