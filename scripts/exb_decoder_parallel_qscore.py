#!/usr/bin/env python

import time
import sys
from collections import Counter
import numpy as np
from sklearn.utils.extmath import cartesian
import string
import gzip
from itertools import islice,izip_longest, izip, groupby
from optparse import OptionParser
import csv
import distance
import dpath
from multiprocessing import Pool, Manager

def grouper(n, iterable):
	#copied from https://coderwall.com/p/sme8mg
	return izip_longest(*[iter(iterable)]*n,fillvalue=None)

def main():
	usage = "usage: exb_decoder_parallel_qscore.py -a read1.fq.gz -b read2.fq.gz -o output -n ncores. Currently only supports gzipped files."
	parser = OptionParser(usage=usage)
	parser.add_option("-a", "--read1fastq", action="store", type="string", dest="read1fastq")
	parser.add_option("-b", "--read2fastq", action="store", type="string", dest="read2fastq")
	parser.add_option("-n", "--numprocs", action="store", type="int", dest="numprocs")
	parser.add_option("-o", "--output", action="store", type="string", dest="outfile")
	(option,args) = parser.parse_args()
	
	if len(args) != 0:
		parser.error("Incorrect # of args")
		
	start_time_master = time.time()
	
	w = csv.writer(open(option.outfile,"w"), delimiter='\t')
	jobs = []
	final_bcs = {}
	results = []
	chunk_size = option.numprocs*100000
	pool = Pool(option.numprocs)

	
	with gzip.open(option.read1fastq) as f1, gzip.open(option.read2fastq) as f2:
		grouper_start_time = time.time()
		
		while True:
			print >>sys.stderr, "Chunking from file..."
			grouper_start_time = time.time()
			r1fq_chunk = list(islice(f1,chunk_size*4))[1::2]
			r2fq_chunk = list(islice(f2,chunk_size*4))[1::2]
			
			#zip with phred information
			z_chunk = zip(list(grouper(2,map(str.strip,r1fq_chunk))), list(grouper(2,map(str.strip,r2fq_chunk))))
			if not z_chunk:
				break
			
			jobs = []
			for chunk in grouper(100000, list(z_chunk)):
				if chunk == None:
					continue
				jobs.append(list(chunk))
				#print len(jobs)
			print >>sys.stderr, "Chunking time was %s seconds" % (time.time() - grouper_start_time)
			results.extend(pool.map(proc_chunk,jobs))
	pool.close()
	#merge tables
	final_bcs = {}
	
	del jobs
	
	for result in results:
		if result == None:
			continue
		for key,val in result.items():
			
			if final_bcs.has_key(key):
				final_bcs[key] += val
			else:
				final_bcs[key] = val
	
	del results
	
	
	#molecular imputation for high error reads
	print >>sys.stderr, "Starting molecular imputator..."
	impute_bc = []
	count_imputes = 0
	num_imputes = len(dpath.util.search(final_bcs,'(*(*)*)'))
	print "Number of ambiguous reads %d" % num_imputes
	for results in dpath.util.search(final_bcs,'(*(*)*)'):
		count_imputes += 1
		if count_imputes % 100000 == 0:
			print "Imputing at %d" % count_imputes
		impute_bc_results = [x for x in results]
		n_impute = 1
		base_bc = []
		accumulator = []
		for x in impute_bc_results:
			if isinstance(x,tuple):
				n_impute *= len(x)
		
		if n_impute >= 64 * 64: #read is kind of messed up; skip
			print >>sys.stderr, "Record has %d possibilities for imputation; skipping..." % n_impute
			del final_bcs[tuple(results)]
			continue
		#copy until reach a list
		for element in impute_bc_results:
			if isinstance(element,int):
				accumulator = []
				if base_bc:
					if isinstance(base_bc[0],list):
						for record in base_bc:
							accumulator.append(record + [element])
				else:
					accumulator.append(base_bc + [element])
				base_bc = accumulator
			if isinstance(element,tuple): #we hit a multielement position
				bc_multiplier = [x for x in element]
				accumulator = []
				for record in bc_multiplier:
					if base_bc:
						if isinstance(base_bc[0],list):
						   for base_bc_record in base_bc:
								accumulator.append(base_bc_record + [record])
						else:
							accumulator.append(base_bc + [record])
					else:
						accumulator.append([record])
				base_bc = accumulator
		#print base_bc
		
		
		#check bctable for each element in base_bc
		found_flag=False
		for record in base_bc:
			#print "Attempting to impute %s with %s" % (tuple(results),tuple(record))
			if tuple(record) in final_bcs and not found_flag: #element exists, append it to that entry
				final_bcs[tuple(record)] += final_bcs[tuple(results)]
				del final_bcs[tuple(results)] #remove the entry 
				found_flag=True
		if not found_flag:
			final_bcs[tuple(base_bc[0])] = final_bcs[tuple(results)] #it doesn't matter which one we pick
			del final_bcs[tuple(results)] #remove it
	
	
	print >>sys.stderr,"Finished imputator, final table with %d entries and %d total molecules" % (len(final_bcs), 0.5 * sum(len(x) for x in final_bcs.values()))
	
	#write to disk
	print >>sys.stderr,"Finished at %d seconds, writing to disk." % (time.time() - start_time_master)
	
	
	
	for key,val in final_bcs.items():
		w.writerow([key,val])
	run_time_master = time.time() - start_time_master
	print >>sys.stderr, "Final dictionary size %d with %d entries in %d seconds" % (len(final_bcs), .5 * sum(len(x) for x in final_bcs.values()), run_time_master)
	
	
#######
	
	
def proc_chunk(chunk):
	start_time = time.time()
	generator = np.matrix(((1,0,0,1,1,1),(0,1,0,1,2,3),(0,0,1,1,3,2)))
	check = np.concatenate((np.mod(-1 * np.matrix(((1,1,1),(1,2,3),(1,3,2))),4), np.matrix(((1,0,0),(0,1,0),(0,0,1)))),axis=1) #check = (-AT | I)

	error_patterns = cartesian(([0,1,2,3],[0,1,2,3],[0,1,2,3],[0,1,2,3],[0,1,2,3],[0,1,2,3]))

	syndromes = {}
	for r in error_patterns:
		syndrome = tuple(np.mod(check*np.matrix(r).transpose(),4).reshape(-1,).tolist()[0])
		if syndrome in syndromes:
			syndromes[syndrome].append(r)
		else:
			syndromes[syndrome] = [r]

	msgs = cartesian(([0,1,2,3],[0,1,2,3],[0,1,2,3]))
	codes = np.mod(msgs * generator,4)
	codes_acgt = []
	for r in codes:
		codes_acgt.append(''.join(str(x) for x in r.tolist()[0]).translate(string.maketrans("0123","ACGT")))
		
	
		
	#start fetching from file
	raws_r1 = []
	raws_r2 = []
	quals_r1 = []
	quals_r2 = []
	syns_r1 = []
	syns_r2 = []
	read_r1 = []
	read_r2 = []
	n_bc = []
	n_r1 = []
	n_r2 = []
	
	num_records = 0
	
	if not chunk:
		return {}
	
	#print "Beginning reading reads at %d seconds." % (time.time() - start_time)
	for record in chunk:
		if record == None:
			continue
		
		(lines_r1_wqual,lines_r2_wqual) = record
		#split into read and qual
		(lines_r1,lines_qual_r1) = lines_r1_wqual
		(lines_r2,lines_qual_r2) = lines_r2_wqual
		#for lines_r1,lines_r2 in chunk:
		if lines_r1 == None or lines_r2 == None:
			break
		num_records += 1
		#print "lines", lines_r1, lines_r2
		if ((distance.levenshtein(lines_r1[12:26],'CGGCGGTGCGCACG') < 7) and (distance.levenshtein(lines_r2[12:26],'CGGCGGTGCGCACG') < 7)):
			#print "unit 1 ok"
			if ((distance.levenshtein(lines_r1[32:46],'CGCCGCCGCAGCGC') < 7) and (distance.levenshtein(lines_r2[32:46],'CGCCGCCGCAGCGC') < 7)):
				#print "unit 2 ok"
				if ((distance.levenshtein(lines_r1[52:71],'AGATGTGTATAAGAGACAG') < 11) and (distance.levenshtein(lines_r2[52:71],'AGATGTGTATAAGAGACAG') < 11)):
					#print "unit 3 ok"
					code_r1_1 = lines_r1[6:12]
					code_r1_2 = lines_r1[26:32]
					code_r1_3 = lines_r1[46:52]
					code_r2_1 = lines_r2[6:12]
					code_r2_2 = lines_r2[26:32]
					code_r2_3 = lines_r2[46:52]
					
					#extract barcodes with N for later
					if ('N' in code_r1_1 + code_r1_2 + code_r1_3 + code_r2_1 + code_r2_2 + code_r2_3):
						n_bc.append([code_r1_1,code_r1_2,code_r1_3,code_r2_1,code_r2_2,code_r2_3])
						n_r1.append(lines_r1)
						n_r2.append(lines_r2)
						continue
					
					read_r1.append(lines_r1)
					read_r2.append(lines_r2)
										
					raws_r1.append([code_r1_1,code_r1_2,code_r1_3])
					quals_r1.append([lines_qual_r1[6:12],lines_qual_r1[26:32],lines_qual_r1[46:52]])
					
					raws_r2.append([code_r2_1,code_r2_2,code_r2_3])
					quals_r2.append([lines_qual_r2[6:12],lines_qual_r2[26:32],lines_qual_r2[46:52]])
					
					syns_r1.append([tuple(np.mod(check*np.matrix([int(x) for x in code_r1_1.translate(string.maketrans("ACGT","0123"))]).transpose(),4).reshape(-1,).tolist()[0]),tuple(np.mod(check*np.matrix([int(x) for x in code_r1_2.translate(string.maketrans("ACGT","0123"))]).transpose(),4).reshape(-1,).tolist()[0]),tuple(np.mod(check*np.matrix([int(x) for x in code_r1_3.translate(string.maketrans("ACGT","0123"))]).transpose(),4).reshape(-1,).tolist()[0])])
					syns_r2.append([tuple(np.mod(check*np.matrix([int(x) for x in code_r2_1.translate(string.maketrans("ACGT","0123"))]).transpose(),4).reshape(-1,).tolist()[0]),tuple(np.mod(check*np.matrix([int(x) for x in code_r2_2.translate(string.maketrans("ACGT","0123"))]).transpose(),4).reshape(-1,).tolist()[0]),tuple(np.mod(check*np.matrix([int(x) for x in code_r2_3.translate(string.maketrans("ACGT","0123"))]).transpose(),4).reshape(-1,).tolist()[0])])
					
					
					
	#print "Building syndromes at %d seconds." % (time.time() - start_time)
	#create the barcode table and process reads and syndromes
	bctable = {}
	for syn_r1,raw_r1,record_r1,syn_r2,raw_r2,record_r2,qual_r1,qual_r2 in zip(syns_r1,raws_r1,read_r1,syns_r2,raws_r2,read_r2,quals_r1,quals_r2):
		#print count
		bctemp = []
		
		for unit,unit_raw,unit_qual in zip(syn_r1,raw_r1,qual_r1):
			zeros = [6-np.count_nonzero(x) for x in syndromes[unit]]
			zeros_count = Counter(zeros)
			#print zeros_count
			if (max(zeros_count.keys()) == 6):
				bctemp.append(codes_acgt.index(unit_raw))
				continue

			corrected = []
			probs = []
			for error_pattern in syndromes[unit]:
				y = [int(x) for x in unit_raw.translate(string.maketrans("ACGT","0123"))]
				x = error_pattern.tolist()
				if (6-np.count_nonzero(x) == max(zeros_count.keys())):
					temp = np.mod([a - b for a,b in zip(y,x)],4)
					corrected.append(codes_acgt.index(''.join([str(x).translate(string.maketrans("0123","ACGT")) for x in temp])))
					p_temp = 1
					for base,basequal in zip(error_pattern,unit_qual):
						if base == 0: #no change
							p_temp = p_temp * (1 - (10 ** (-0.1 * (ord(basequal)-33))))
						else: #base change
							p_temp = p_temp * (10 ** (-0.1 * (ord(basequal)-33)))
					probs.append(p_temp)
			if len(corrected) == 1:
				bctemp.append(corrected[0])
			else:
				#calculate maximums to impute
				max_prob = max(probs)
				max_prob_pos = [i for i,j in enumerate(probs) if j == max_prob]
				if len(max_prob_pos) == 1:
					bctemp.append(corrected[max_prob_pos[0]])
				else:
					bctemp.append(tuple([corrected[i] for i in max_prob_pos]))
				#bctemp.append(tuple(corrected))

		for unit,unit_raw,unit_qual in zip(syn_r2,raw_r2,qual_r2):
			zeros = [6-np.count_nonzero(x) for x in syndromes[unit]]
			zeros_count = Counter(zeros)
			#print zeros_count
			if (max(zeros_count.keys()) == 6):
				bctemp.append(codes_acgt.index(unit_raw))
				continue

			corrected = []
			probs = []
			
			for error_pattern in syndromes[unit]:
				y = [int(x) for x in unit_raw.translate(string.maketrans("ACGT","0123"))]
				x = error_pattern.tolist()
				if (6-np.count_nonzero(x) == max(zeros_count.keys())):
					temp = np.mod([a - b for a,b in zip(y,x)],4)
					corrected.append(codes_acgt.index(''.join([str(x).translate(string.maketrans("0123","ACGT")) for x in temp])))
					p_temp = 1
					for base,basequal in zip(error_pattern,unit_qual):
						if base == 0: #no change
							p_temp = p_temp * (1 - (10 **(-0.1 * (ord(basequal)-33))))
						else: #base change
							p_temp = p_temp * (10 **(-0.1 * (ord(basequal)-33)))
					probs.append(p_temp)
			if len(corrected) == 1:
				bctemp.append(corrected[0])
			else:
				#calculate maximums to impute
				max_prob = max(probs)
				max_prob_pos = [i for i,j in enumerate(probs) if j == max_prob]
				if len(max_prob_pos) == 1:
					bctemp.append(corrected[max_prob_pos[0]])
				else:
					bctemp.append(tuple([corrected[i] for i in max_prob_pos]))
				#bctemp.append(tuple(corrected))
		if bctable.has_key(tuple(bctemp)):
			bctable[tuple(bctemp)] += [record_r1.strip(), record_r2.strip()]
		else:
			bctable[tuple(bctemp)] = [record_r1.strip(), record_r2.strip()]
	
	
	
	#start the N imputator
	#print >>sys.stderr, "Starting N imputator at %d seconds" % (time.time() - start_time)
	
	for bc_record,record_r1,record_r2 in zip(n_bc,n_r1,n_r2):
		bctemp = []
		for unit in bc_record:
			accumulator = []
			dists = [int(distance.levenshtein(unit,x)) for x in codes_acgt]
			dists_count = Counter(dists)
			for record in codes_acgt:
				if distance.levenshtein(unit,record) == min(dists_count.keys()):
					accumulator.append(codes_acgt.index(record))
			if len(accumulator) == 1:
				bctemp.append(accumulator[0])
			else:
				bctemp.append(tuple(accumulator))
		if bctable.has_key(tuple(bctemp)):
			bctable[tuple(bctemp)] += [record_r1.strip(), record_r2.strip()]
		else:
			bctable[tuple(bctemp)] = [record_r1.strip(), record_r2.strip()]
	
	run_time = time.time() - start_time
	
	print >>sys.stderr, "Processed a chunk with %d records in %d seconds; dict size %d with %d values" % (num_records, run_time, len(bctable), 0.5 * sum(len(x) for x in bctable.values()))

	return bctable
	
if __name__ == "__main__":
	main()

