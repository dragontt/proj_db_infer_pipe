#!/usr/bin/python

"""
Compute Spearman's rank order correlation of the NetProphet score rankings and the 
database (CIS-BP) motif-target scan score rankings.
"""

import sys
import argparse
import glob
import os.path
sys.path.append('/home/mblab/ykang/proj_db_infer_pipe/scripts/')
import compt_pw_corr
import time


def parse_args(argv):
    parser = argparse.ArgumentParser(description="Compute rank order correlations of np scores and datrabase motif-target score rankings")
    parser.add_argument('-n', '--dir_np', dest='dir_np', help='NetProphet rankings')
    parser.add_argument('-l', '--fn_tf_list', dest='fn_tf_list')
    parser.add_argument('-d', '--dir_database', dest='dir_database', help="Database motif-target score rankings")
    parser.add_argument('-o', '--dir_output', dest='dir_output')
    parsed = parser.parse_args(argv[1:])
    return parsed


def main(argv):
    parsed = parse_args(argv)

    use_tf_subset = False if parsed.fn_tf_list == None else True
    parsed.dir_np = check_dir(parsed.dir_np)
    parsed.dir_database = check_dir(parsed.dir_database)
    parsed.dir_output = check_dir(parsed.dir_output)

    # parse the list of query motifs and their rankings
    if (use_tf_subset):
        query_tf_fns = []
        lines = open(parsed.fn_tf_list, "r").readlines()
        for line in lines:
            query_tf_fns.append(parsed.dir_np + line.strip())
    else:
        query_tf_fns = glob.glob(parsed.dir_np + "/*")    
    
    [query_tfs, ranks_query, query_tfs_skip] = parse_all_ranks(query_tf_fns)
    print "Valid TF count:", len(query_tfs), "; Skipped TF:", query_tfs_skip

    # parse the list of database motifs and their rankings	
    db_motif_fns = glob.glob(parsed.dir_database + "/*")
    [db_motifs, db_ranks, db_motifs_skip] = parse_all_ranks(db_motif_fns)
    print "CIS-BP motif count:", len(db_motifs)

    # compute ranked list correlations
    for i in range(len(query_tfs)):
    	tic = time.clock()

    	writer = open(parsed.dir_output + query_tfs[i], "w")
        for j in range(len(db_motifs)):
            corr = compt_pw_corr.compute_corr(db_ranks[j], ranks_query[i])
            writer.write("%s\t%.5f\n" % (db_motifs[j], corr))
        writer.close()
        
        toc = time.clock()
    	print query_tfs[i], i+1, "/", len(query_tfs), "Time elapsed:", toc-tic, "sec"


def check_dir(dirname):
    if not dirname.endswith('/'):
        dirname += '/'
    return dirname


def parse_all_ranks(filenames):
    queries_valid = []
    queries_skip = []
    ranks = []
    for i in range(len(filenames)):
        query = os.path.basename(filenames[i])
        if os.stat(filenames[i]).st_size == 0:
            queries_skip.append(query)
        else:
           queries_valid.append(query)
           ranks.append(parse_rank(filenames[i]))
    return [queries_valid, ranks, queries_skip]


def parse_rank(filename):
    rank = {}
    lines = open(filename, 'r').readlines()
    for line in lines:
        linesplit = line.split()
        rank[linesplit[0]] = linesplit[1]
    return rank


if __name__ == "__main__":
    main(sys.argv)
