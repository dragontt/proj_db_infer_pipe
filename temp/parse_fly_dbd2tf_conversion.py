#!/usr/bin/python

lines = open('/home/mblab/ykang/proj_db_infer_pipe/resources/fly_base/fbgn_fbtr_fbpp_fb_2015_01.tsv','r').readlines()
fbpp2fbgn = {}
for i in range(6, len(lines)-1):
	line = lines[i].strip().split()
	if len(line) > 2:
		fbpp2fbgn[line[2]] = line[0]

writer = open('/home/mblab/ykang/proj_db_infer_pipe/resources/fly_aa_seq/pids.dbd2fbgn','w')
lines = open('/home/mblab/ykang/proj_db_infer_pipe/resources/fly_aa_seq/pids.dbd','r').readlines()
for i in range(len(lines)):
	line = lines[i].strip()
	fbgn = fbpp2fbgn[line.split('.')[0]]
	writer.write('%s\t%s\n' % (line, fbgn))
writer.close()

