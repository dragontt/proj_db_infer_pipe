#!/usr/bin/bash

function compute_paired_pid {
	QUERY_MOTIF=$1
	DIR_QUERY_DBDS=$2
	DIR_DATABASE_DBDS=$3
	DIR_OUTPUT=$4

	DIR_PROJ=/home/mblab/ykang/proj_db_infer_pipe

	# create temp directories
	# DIR_TEMP=$DIR_PROJ/temp/${QUERY_MOTIF}_dbd_pid
	DIR_TEMP=/tmp/${QUERY_MOTIF}_dbd_pid
	mkdir -p $DIR_TEMP

	# create files of paired query DBD to individual datatbase DBD
	python $DIR_PROJ/scripts/generate_paired_fasta.py -m $QUERY_MOTIF -f1 $DIR_QUERY_DBDS -f2 $DIR_DATABASE_DBDS -o $DIR_TEMP

	# align and compute percent identiy of each dbd pair
	for filename_full in $DIR_TEMP/*.fasta; do
		filename_partial=$(basename $filename_full)
		# echo $filename_partial
		$HOME/usr/clustalo --infile ${filename_full} --outfile $DIR_TEMP/${QUERY_MOTIF}_${filename_partial}.out --seqtype protein --distmat-out $DIR_TEMP/${QUERY_MOTIF}_${filename_partial}.pim --full --percent-id --force
	done

	# parse all percent identity matrix files, and output a ordered list of aadis
	python $DIR_PROJ/scripts/parse_paired_aaid.py -i $DIR_TEMP -o $DIR_OUTPUT/$QUERY_MOTIF

	rm -r $DIR_TEMP
}

DIR_RES=$HOME/proj_db_infer_pipe/resources/cisbp_1.01
# DIR_OUT=$HOME/proj_db_infer_pipe/output/fly_cisbp_known_motif_dbd_pid
# DIR_OUT=$HOME/proj_db_infer_pipe/output/yeast_cisbp_known_motif_dbd_pid
DIR_OUT=$HOME/proj_db_infer_pipe/output/cisbp_known_motif_dbd_pid
# list=$DIR_RES/cisbp_motifs_Drosophila_melanogaster.txt
# list=$DIR_RES/cisbp_motifs_Saccharomyces_cerevisiae.txt
list=$DIR_RES/cisbp_motifs_all_species.txt

counter=0
while read line; do
	counter=$[$counter +1]
	if [ -e $DIR_OUT/$line ]; then
		echo "Processing $line ... Existed ... $counter ... Done"
	else
		echo -n "Processing $line ... New ... $counter"
		# compute_paired_pid $line $DIR_RES/cisbp.dbd.fasta.Drosophila_melanogaster $DIR_RES/cisbp.dbd.fasta.Drosophila_melanogaster $DIR_OUT
		# compute_paired_pid $line $DIR_RES/cisbp.dbd.fasta.Saccharomyces_cerevisiae $DIR_RES/cisbp.dbd.fasta.Saccharomyces_cerevisiae $DIR_OUT
		compute_paired_pid $line $DIR_RES/cisbp.dbd.fasta $DIR_RES/cisbp.dbd.fasta $DIR_OUT
		echo " ... Done"
	fi
done < $list

