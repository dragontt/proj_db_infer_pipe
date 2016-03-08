#! /bin/bash
# Input variables
IN_TF_LIST=$1 	# list of tf names
IN_TF_PWM=$2	# directory of tf pwm
IN_PROMOTERS=$3	# promoter sequence file (e.g. yeast_promoter_seq/s_cerevisiae.promoters.fasta, fly_promoter_seq/rsat_dmel_upstream_-2000_+200.filtered.fasta)
IN_BACKGROUND=$4	# background frequency file (e.g. Saccharomyces_cerevisiae.bmf, Drosophila_melanogaster.bmf)
OUT_FIMO=$5		# directory of fimo alignment output 

counter=0

# while read -a line
# do
# 	counter=$[$counter +1]
# 	motif=${line[0]}
# 	if [ -e $IN_TF_PWM/$motif ]
# 		then
# 		echo  "*** Processing $motif ... $counter"
# 		fimo -o $OUT_FIMO/$motif --thresh 5e-3 --bgfile $IN_BACKGROUND $IN_TF_PWM/$motif $IN_PROMOTERS
# 		sed ' 1d ' $OUT_FIMO/$motif/fimo.txt | cut -f 1,2,7 > $OUT_FIMO/$motif/temp.txt
# 		ruby $HOME/proj_db_infer_pipe/scripts/estimate_affinity.rb -i $OUT_FIMO/$motif/temp.txt > $OUT_FIMO/$motif.summary
# 		sed ' 1d ' $OUT_FIMO/$motif/fimo.txt | sort -t $'\t' -k 1,1 | cut -f 1,2,3,4,5 > $OUT_FIMO/$motif.loci
# 		# rm -r $OUT_FIMO/$motif
# 		echo "*** Done"
# 	else
# 		echo  "*** No PWM for $motif ... $counter *** Done"
# 	fi
# done < $IN_TF_LIST

while read -a line
do
	## process fimo scan
	# counter=$[$counter +1]
	# motif=${line[0]}
	# echo  "*** Processing $motif ... $counter"
	# fimo -o $OUT_FIMO/$motif --thresh 5e-3 --bgfile $IN_BACKGROUND $IN_TF_PWM/$motif $IN_PROMOTERS
	
	# sed ' 1d ' $OUT_FIMO/$motif/fimo.txt | cut -f 1,2,7 > $OUT_FIMO/$motif/temp.txt
	# ruby $HOME/proj_db_infer_pipe/scripts/estimate_affinity.rb -i $OUT_FIMO/$motif/temp.txt > $OUT_FIMO/${motif}.summary

	# rm $OUT_FIMO/$motif/cisml.css
	# rm $OUT_FIMO/$motif/cisml.xml
	# rm $OUT_FIMO/$motif/fimo.gff
	# rm $OUT_FIMO/$motif/fimo.html
	# rm $OUT_FIMO/$motif/fimo-to-html.xsl
	# rm $OUT_FIMO/$motif/fimo.xml

	# post process masked fimo scores
	counter=$[$counter +1]
	motif=${line[0]}
	echo  "*** Processing $motif ... $counter"

	python $HOME/proj_db_infer_pipe/scripts/mask_fimo.py -m mask1 -i $OUT_FIMO/$motif/fimo.txt -a $HOME/proj_db_infer_pipe/resources/fly_promoter_seq/dmel_-2000_+200.rsat_promoters.cons_elem_mask.txt -o $OUT_FIMO/$motif/fimo_mask1.txt
	sed ' 1d ' $OUT_FIMO/$motif/fimo_mask1.txt | cut -f 1,2,7 > $OUT_FIMO/$motif/temp_mask1.txt
	ruby $HOME/proj_db_infer_pipe/scripts/estimate_affinity.rb -i $OUT_FIMO/$motif/temp_mask1.txt > $OUT_FIMO/${motif}.summary_mask1

	python $HOME/proj_db_infer_pipe/scripts/mask_fimo.py -m mask2 -i $OUT_FIMO/$motif/fimo.txt -a $HOME/proj_db_infer_pipe/resources/fly_promoter_seq/dmel_-2000_+200.rsat_promoters.cons_elem_mask.txt -o $OUT_FIMO/$motif/fimo_mask2.txt
	sed ' 1d ' $OUT_FIMO/$motif/fimo_mask2.txt | cut -f 1,2,7 > $OUT_FIMO/$motif/temp_mask2.txt
	ruby $HOME/proj_db_infer_pipe/scripts/estimate_affinity.rb -i $OUT_FIMO/$motif/temp_mask2.txt > $OUT_FIMO/${motif}.summary_mask2

	python $HOME/proj_db_infer_pipe/scripts/mask_fimo.py -m mask3 -i $OUT_FIMO/$motif/fimo.txt -a $HOME/proj_db_infer_pipe/resources/fly_promoter_seq/dmel_-2000_+200.rsat_promoters.cons_elem_mask.txt -o $OUT_FIMO/$motif/fimo_mask3.txt
	sed ' 1d ' $OUT_FIMO/$motif/fimo_mask3.txt | cut -f 1,2,7 > $OUT_FIMO/$motif/temp_mask3.txt
	ruby $HOME/proj_db_infer_pipe/scripts/estimate_affinity.rb -i $OUT_FIMO/$motif/temp_mask3.txt > $OUT_FIMO/${motif}.summary_mask3

	python $HOME/proj_db_infer_pipe/scripts/mask_fimo.py -m mask3 -i $OUT_FIMO/$motif/fimo.txt -a $HOME/proj_db_infer_pipe/resources/fly_promoter_seq/dmel_-2000_+200.rsat_promoters.cons_threshold_0.05_mask.txt -o $OUT_FIMO/$motif/fimo_mask3_cons_thd_0.05.txt
	sed ' 1d ' $OUT_FIMO/$motif/fimo_mask3_cons_thd_0.05.txt | cut -f 1,2,7 > $OUT_FIMO/$motif/temp_mask3_cons_thd_0.05.txt
	ruby $HOME/proj_db_infer_pipe/scripts/estimate_affinity.rb -i $OUT_FIMO/$motif/temp_mask3_cons_thd_0.05.txt > $OUT_FIMO/${motif}.summary_mask3_cons_thd_0.05

	echo "*** Done"
done < $IN_TF_LIST

echo "*** ALL DONE! ***"
