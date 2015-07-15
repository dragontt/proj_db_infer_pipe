#! /bin/bash
# bash align_fimo.sh resources/scertf_names.txt resources/scertf_pwms output/scertf Saccharomyces_cerevisiae

# Input variables
IN_TF_LIST=$1 	# list of tf names
IN_TF_PWM=$2	# directory of tf pwm
OUT_FIMO=$3		# directory of fimo alignment output 
# OUT_RANK=$4		# direcotry of the ranked lists of fimo output

# Algin pfms individually 
PROJ_DIR=$HOME/proj_db_infer_pipe

mkdir -p $PROJ_DIR/$OUT_FIMO
# mkdir -p $PROJ_DIR/$OUT_RANK

counter=0

while read -a line
do
	counter=$[$counter +1]
	motif=${line[0]}
	if [ -e $PROJ_DIR/$IN_TF_PWM/$motif ]
		then
		echo  "*** Processing $motif ... $counter"
		cd $HOME/usr/meme/bin
		./fimo -o $PROJ_DIR/$OUT_FIMO/$motif --thresh 5e-3 --bgfile $PROJ_DIR/resources/cisbp_all_species_bg_freq/Drosophila_melanogaster.bmf $PROJ_DIR/$IN_TF_PWM/$motif $PROJ_DIR/resources/fly_promoter_seq/rsat_dmel_upstream_-2000_+200.filtered.fasta 
		#./fimo -o $PROJ_DIR/$OUT_FIMO/$motif --thresh 5e-3 --bgfile $PROJ_DIR/resources/cisbp_all_species_bg_freq/Drosophila_melanogaster.bmf $PROJ_DIR/$IN_TF_PWM/$motif $PROJ_DIR/resources/crypto_promoter_seq/crNeoH99.promoters.upstream600.fasta
		sed ' 1d ' $PROJ_DIR/$OUT_FIMO/$motif/fimo.txt | cut -f 1,2,7 > $PROJ_DIR/$OUT_FIMO/$motif/temp.txt
		ruby $PROJ_DIR/scripts/estimate_affinity.rb -i $PROJ_DIR/$OUT_FIMO/$motif/temp.txt > $PROJ_DIR/$OUT_FIMO/$motif.summary
		# mv $PROJ_DIR/$OUT_FIMO/$motif/fimo.txt $PROJ_DIR/$OUT_FIMO/$motif.fimo
		rm -r $PROJ_DIR/$OUT_FIMO/$motif
		echo "*** Done"
	else
		echo  "*** No PWM for $motif ... $counter *** Done"
	fi
done < $PROJ_DIR/$IN_TF_LIST

# Compute the rankded lists
# echo "*** Computing rankings ... "
# python $PROJ_DIR/scripts/compt_fimo_rank.py -i $PROJ_DIR/$OUT_FIMO -t $PROJ_DIR/resources/crypto_network_orig/gids -o $PROJ_DIR/$OUT_RANK
# echo "*** Done"

echo "*** ALL DONE! ***"
