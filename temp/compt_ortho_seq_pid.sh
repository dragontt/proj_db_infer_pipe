#! /bin/bash
# Compute sequence percent identity of one query gene's sequence to those of all ortholog genes. 

DIR_IN=$1
DIR_OUT=$2

for FILE in $DIR_IN/*
do
	FN_OUT=$DIR_OUT/$(basename $FILE)
	clustalo --infile $FILE --outfile ${FN_OUT}.out --seqtype dna --distmat-out ${FN_OUT}.pim --full --percent-id --force
done
