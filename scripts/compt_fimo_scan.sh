#! /bin/bash
# bash align_fimo.sh resources/scertf_names.txt resources/scertf_pwms output/scertf Saccharomyces_cerevisiae

# Input variables
IN_TF_LIST=$1 	# list of tf names
IN_TF_PWM=$2	# directory of tf pwm
IN_PROMOTERS=$3	# promoter sequence file (e.g. yeast_promoter_seq/s_cerevisiae.promoters.fasta, fly_promoter_seq/rsat_dmel_upstream_-2000_+200.filtered.fasta)
IN_BACKGROUND=$4	# background frequency file (e.g. Saccharomyces_cerevisiae.bmf, Drosophila_melanogaster.bmf)
OUT_FIMO=$5		# directory of fimo alignment output 

counter=0

while read -a line
do
	counter=$[$counter +1]
	motif=${line[0]}
	if [ -e $IN_TF_PWM/$motif ]
		then
		echo  "*** Processing $motif ... $counter"
		fimo -o $OUT_FIMO/$motif --thresh 5e-3 --bgfile $IN_BACKGROUND $IN_TF_PWM/$motif resources/$IN_PROMOTERS
		sed ' 1d ' $OUT_FIMO/$motif/fimo.txt | cut -f 1,2,7 > $OUT_FIMO/$motif/temp.txt
		ruby $HOME/proj_db_infer_pipe/scripts/estimate_affinity.rb -i $OUT_FIMO/$motif/temp.txt > $OUT_FIMO/$motif.summary
		sed ' 1d ' $OUT_FIMO/$motif/fimo.txt | sort -t $'\t' -k 1,1 | cut -f 1,2,3,4,5 > $OUT_FIMO/$motif.loci
		rm -r $OUT_FIMO/$motif
		echo "*** Done"
	else
		echo  "*** No PWM for $motif ... $counter *** Done"
	fi
done < $IN_TF_LIST

echo "*** ALL DONE! ***"
