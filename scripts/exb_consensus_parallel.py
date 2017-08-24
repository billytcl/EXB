#!/usr/bin/env python

import sys;
from optparse import OptionParser
import csv
import re
from collections import Counter
import numpy
from multiprocessing import Pool
from itertools import izip_longest, islice
import time
import math
import string

csv.field_size_limit(sys.maxsize)

def main():
	usage = "usage: exb_consensus_parallel -i input decoder file -o output -n numprocs."
	
	parser = OptionParser(usage=usage)
	parser.add_option("-i", "--input", action="store", type="string", dest="infile")
	parser.add_option("-o", "--output", action="store", type="string", dest="outfile")
	parser.add_option("-n", "--numprocs", action="store", type="int", dest="numprocs")
	
	(option,args) = parser.parse_args()
	
	if len(args) != 0:
		parser.error("Incorrect # of args")
		
	pool = Pool(option.numprocs)
	
	results = []

	#r = csv.reader(open(option.infile,"r"), delimiter='\t')
	
	jobs = []
	
	with open(option.infile,"r") as f:
		grouper_start_time = time.time()
		
		while True:
			print >>sys.stderr, "Chunking from file..."
			grouper_start_time = time.time()
			read_chunk = list(islice(f,1000000))
			
			
			if not read_chunk:
				break
			
			jobs = []
			for chunk in grouper(int(math.ceil(1000000/option.numprocs)),read_chunk):
				if chunk == None:
					continue
				jobs.append(chunk)
			print >>sys.stderr, "Chunking time was %s seconds" % (time.time() - grouper_start_time)
			results.extend(pool.map(proc_chunk,jobs))
	
	w = open(option.outfile,"w")
	
	for result_row in results:
		w.writelines(result_row)
	
	w.close()
	
def grouper(n, iterable):
	#copied from https://coderwall.com/p/sme8mg
	return izip_longest(*[iter(iterable)]*n,fillvalue=None)
	
def proc_chunk(chunk):
	output = []
	start_time = time.time()
	
	for row in chunk:
		quickPerf = False
		if row == None:
			continue #taking care of padding from iterator
		row = row.split('\t')
		
		#could have read in as a dict...
		row[0] = re.sub('[\(\'\)\,\[\]]','',row[0])
		row[1] = re.sub('[\(\'\)\,\[\]]','',row[1])
		seqs = row[1].split()
		seq_f = seqs[0::2]
		seq_r = seqs[1::2]
				
		for i in range(len(seq_f)):
			seq_f[i] = seq_f[i].ljust(151)
		
		for i in range(len(seq_r)):
			seq_r[i] = seq_r[i].ljust(151)
		
		
		consensus_f = [];
		consensus_f_stats = [];
		consensus_r = [];
		consensus_r_stats = [];
		
		#quick check for perfect duplicates
		acc_f = seq_f[0]
		acc_r = seq_r[0]
		
		for x in seq_f:
			if not (sum([ord(a) ^ ord(b) for a,b in zip(acc_f,x)]) == 0):
				quickPerf = False
				break
			if not (sum([ord(a) ^ ord(b) for a,b in zip(acc_r,x)]) == 0):
				quickPerf = False
				break
			quickPerf = True
			
		if quickPerf:
			consensus_f = acc_f
			consensus_f_stats = '0' * len(acc_f)
			consensus_r = acc_r
			consensus_r_stats = '0' * len(acc_r)
			output.append("%d\t%s\t%s\t%s\t%s\n" % (len(seq_f), consensus_f, consensus_f_stats, consensus_r, consensus_r_stats))
			continue
		
		
		for i in range(max(len(x) for x in seq_f)):
			bases_f = [x[i] for x in seq_f]
			bases_f = filter(string.strip, bases_f)
			if bases_f == []:
				continue
			c = Counter(bases_f)
			c_stats = sorted([x[1] for x in c.most_common()],reverse=True)
			
			if len(c_stats) == 1:
				consensus_f.append(c.most_common(1)[0][0])
				consensus_f_stats.append(str(0))
			elif c_stats[0] == c_stats[1]:
				consensus_f.append('N')
				#consensus_f_stats.append(str(2))
				consensus_f_stats.append('X')
			else:
				consensus_f.append(c.most_common(1)[0][0])
				c_len = len(c)
				if 'N' in c and c_len > 4:
					c_len = 4
				#consensus_f_stats.append(str(1))
				consensus_f_stats.append(str(c_len))
			
		consensus_f = ''.join(consensus_f)	
		consensus_f_stats = ''.join(consensus_f_stats)

		
		for i in range(max(len(x) for x in seq_r)):
			bases_r = [x[i] for x in seq_r]
			bases_r = filter(string.strip,bases_r)
			if bases_r == []:
				continue

			c = Counter(bases_r)
			c_stats = sorted([x[1] for x in c.most_common()],reverse=True)
			
			
			if len(c_stats) == 1:
				consensus_r.append(c.most_common(1)[0][0])
				consensus_r_stats.append(str(0))
			elif c_stats[0] == c_stats[1]:
				consensus_r.append('N')
				consensus_r_stats.append('X')
				#consensus_r_stats.append(str(2))
			else:
				consensus_r.append(c.most_common(1)[0][0])
				c_len = len(c)
				if 'N' in c and c_len > 4:
					c_len = 4
				#consensus_r_stats.append(str(1))
				consensus_r_stats.append(str(c_len))
		
		consensus_r = ''.join(consensus_r)	
		consensus_r_stats = ''.join(consensus_r_stats)

		
		output.append("%d\t%s\t%s\t%s\t%s\t%s\n" % (len(seq_f), row[0], consensus_f, consensus_f_stats, consensus_r, consensus_r_stats))
		
	#output[-1] = output[-1].strip() 
	#get rid of last newline
	
	run_time = time.time() - start_time
	print >>sys.stderr, "Processed a chunk of size %d in %d seconds." % (len(chunk), run_time)
	return output
	
if __name__ == "__main__":
	main()

		
