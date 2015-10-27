#!/bin/bash
# SPECIES=yeast
# SPECIES_FULL=Saccharomyces_cerevisiae
# SPECIES=fly
# SPECIES_FULL=Drosophila_melanogaster
SPECIES_FULL=all_species

DIR_PROJ=/home/mblab/ykang/proj_db_infer_pipe
DIR_PFM=$DIR_PROJ/resources/cisbp_1.01/cisbp_all_pfms
FN_TARGET_PFM=$DIR_PROJ/resources/cisbp_1.01/cisbp_pfms.${SPECIES_FULL}.meme
FN_QUERY_LIST=$DIR_PROJ/resources/cisbp_1.01/cisbp_motifs_${SPECIES_FULL}.txt
# FN_BACKGROUND=$DIR_PROJ/resources/cisbp_all_species_bg_freq/${SPECIES_FULL}.bmf

## specific species, e.g. yeast, fly
# DIR_OUTPUT=$DIR_PROJ/output/${SPECIES}_cisbp_known_motif_pwm_similarity

## all species with direct evidnce in cisbp
DIR_OUTPUT=$DIR_PROJ/output/cisbp_known_motif_pwm_similarity

while read -r motif_name
do
	echo -n "tomtom comparing $motif_name ... "
	# tomtom -no-ssc -oc $DIR_OUTPUT/$motif_name -eps -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 1 -bfile $FN_BACKGROUND $DIR_PFM/$motif_name $FN_TARGET_PFM
	tomtom -no-ssc -oc $DIR_OUTPUT/$motif_name -eps -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 1 $DIR_PFM/$motif_name $FN_TARGET_PFM
	sed -n '2,$p' $DIR_OUTPUT/$motif_name/tomtom.txt | awk '{print $2, $5}' > $DIR_OUTPUT/$motif_name.tomtom
	# rm -r $DIR_OUTPUT/$motif_name/
	echo "done"
done < $FN_QUERY_LIST
