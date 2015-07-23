#/usr/bin/python
import sys
import argparse

""" Compute the average conservation scores of the genes, using the pre-aligned conservation score file (src: http://hgdownload.soe.ucsc.edu/downloads.html#species) """

seq_types = ['genome', 'exon']

def parse_args(argv):
	parser = argparse.ArgumentParser(description='Compute gene conservation scores.')
	parser.add_argument('-t', '--seq_type', dest='seq_type', help='options: %s' % seq_types)
	pasrer.add_argument('-w', '--fn_conserv_scores', dest='fn_conserv_scores', help='format: wigFix')
	parser.add_argument('-f', '--fn_seqs', dest='fn_seqs')
	parser.add_argument('-g', '--fn_gids', dest='fn_gids')
	parser.add_argument('-o', '--fn_output', dest='fn_output')
	parsed = parser.parse_args(argv[1:])
	return parsed

def main(argv):
	parsed = parse_args(argv)

	

if __name__ == "__main__":
    main(sys.argv)