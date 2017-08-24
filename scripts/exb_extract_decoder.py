#!/usr/bin/env python

import sys
from optparse import OptionParser
import csv
import string

def main():
	usage = "usage: infinity_extract -i consensus_file -o outputfile.fa."
	parser = OptionParser(usage=usage)
	parser.add_option("-i", "--input", action="store", type="string", dest="infile")
	
	(option,args) = parser.parse_args()
	
	if len(args) != 0:
		parser.error("Incorrect # of args")
		
	r = csv.reader(open(option.infile,"r"), delimiter='\t')
	
	
	for row in r:
		#print read1
		print ">B:%s-D:%d-R:%s/1\n%s" % (row[1].translate(string.maketrans("+/ ","-__")), int(row[0]), row[2][0:6] + row[4][0:6], row[2][71::])
		#print read2
		print ">B:%s-D:%d-R:%s/2\n%s" % (row[1].translate(string.maketrans("+/ ","-__")), int(row[0]), row[2][0:6] + row[4][0:6], row[4][71::])


if __name__ == "__main__":
	main()

		
