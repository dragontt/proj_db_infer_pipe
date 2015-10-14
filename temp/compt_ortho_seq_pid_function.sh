while read rid
do
	FILE=$1/$rid
	FN_OUT=$3/$(basename $FILE)
	clustalo --infile $FILE --outfile ${FN_OUT}.out --seqtype dna --distmat-out ${FN_OUT}.pim --full --percent-id --force
done < $2
