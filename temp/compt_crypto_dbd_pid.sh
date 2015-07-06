#! /usr/bin/bash

function compt_pairwise_dbd_pid {
	QUERY_MOTIF=$1
	DIR_QUERY_DBDS=$2
	DIR_DATABASE_DBDS=$3

	DIR_PROJ=$HOME/proj_db_infer_pipe

	# create temp directories
	# DIR_TEMP=$DIR_PROJ/temp/${QUERY_MOTIF}_dbd_pid
	DIR_TEMP=/tmp/${QUERY_MOTIF}_dbd_pid
	mkdir -p $DIR_TEMP

	# create files of paired query DBD to individual datatbase DBD
	python $DIR_PROJ/temp/generate_paired_fasta.py -m $QUERY_MOTIF -f1 $DIR_QUERY_DBDS -f2 $DIR_DATABASE_DBDS -o $DIR_TEMP

	# align and compute percent identiy of each dbd pair
	for filename_full in $DIR_TEMP/*.fasta; do
		filename_partial=$(basename $filename_full)
		# echo $filename_partial
		$HOME/usr/clustalo --infile ${filename_full} --outfile $DIR_TEMP/${QUERY_MOTIF}_${filename_partial}.out --seqtype protein --distmat-out $DIR_TEMP/${QUERY_MOTIF}_${filename_partial}.pim --full --percent-id --force
	done

	# parse all percent identity matrix files, and output a ordered list of aadis
	python $DIR_PROJ/scripts/parse_paired_aaid.py -i $DIR_TEMP -o $DIR_PROJ/resources/crypto_aa_seq/clustal_dbd_pairwise/$QUERY_MOTIF

	rm -r $DIR_TEMP
}

tf_list=$1

DIR_RES=$HOME/proj_db_infer_pipe/resources/crypto_aa_seq
list=$DIR_RES/${tf_list}
# list=$DIR_RES/crypto_network_orig/regulator.gids
counter=0
while read line; do
	counter=$[$counter +1]
	if [ -e $DIR_RES/clustal_dbd_pairwise/$line ]; then
		echo "Processing $line ... Existed ... $counter ... Done"
	else
		echo -n "Processing $line ... New ... $counter"
		compt_pairwise_dbd_pid $line $DIR_RES/crNeoH99.peptides.lit.regulator.dbd.fasta $DIR_RES/crNeoH99.peptides.lit.regulator.dbd.fasta
		echo " ... Done"
	fi
done < $list
