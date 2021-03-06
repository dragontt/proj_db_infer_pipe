#! /bin/bash

# Compute DBD percent identity of one query motif's DBD to all motif DBDs in the database. 
# 1. Compile the fasta files and make a temp file of paired query DBD to individual databse DBD;
# 2. Compute the percent idenity of those paired DBDs;
# 3. Create a ordered list of aaids of the query motif DBD to all database motif DBDs.

QUERY_MOTIF=$1
DIR_PROJ=$2
DIR_QUERY_DBDS=$3
DIR_DATABASE_DBDS=$4
DIR_OUTPUT=$5
IS_CISBP=$6

# create temp directories
# DIR_TEMP=$DIR_PROJ/temp/${QUERY_MOTIF}_dbd_pid
DIR_TEMP=/tmp/${QUERY_MOTIF}_dbd_pid
mkdir -p $DIR_TEMP

# create files of paired query DBD to individual datatbase DBD
python $DIR_PROJ/scripts/generate_paired_fasta.py -m $QUERY_MOTIF -f1 $DIR_QUERY_DBDS -f2 $DIR_DATABASE_DBDS -o $DIR_TEMP -c $IS_CISBP

# align and compute percent identiy of each dbd pair
for filename_full in $DIR_TEMP/*.fasta; do
	filename_partial=$(basename $filename_full)
	# echo $filename_partial
	$HOME/usr/clustalo --infile ${filename_full} --outfile $DIR_TEMP/${QUERY_MOTIF}_${filename_partial}.out --seqtype protein --distmat-out $DIR_TEMP/${QUERY_MOTIF}_${filename_partial}.pim --full --percent-id --force
done

# parse all percent identity matrix files, and output a ordered list of aadis
python $DIR_PROJ/scripts/parse_paired_aaid.py -i $DIR_TEMP -o $DIR_OUTPUT/$QUERY_MOTIF

rm -r $DIR_TEMP
