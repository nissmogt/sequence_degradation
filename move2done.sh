#!/bin/bash
while read line
do
mv aln/${line}.fa done/
done < "$1"
