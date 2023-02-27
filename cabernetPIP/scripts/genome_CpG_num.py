#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import pysam
import re

def find_dinucleotide (seqin, di_n):
    '''
    find first pos of dinucleotides (such as CpG or GpC) in the seqin.
    di_n = "CG" or di_n = "CG" and so on.
    di_pos = [3, 13]
    '''
    sta = 0
    di_pos = []
    while sta < len(seqin):
        aa = seqin.find(di_n,sta)
        if aa == -1:
            return di_pos
            break
        di_pos.append(aa)
        sta=aa+1
    return di_pos


	

if __name__ == "__main__":
	genome_reference = sys.argv[1]
	dinucleotide = sys.argv[2]
	file_out = sys.argv[3]

	fastafile = pysam.Fastafile(genome_reference)	

	di_pos = []
	file_o=open(file_out,'w')

	references = fastafile.references
	for cc in references:
	#	if re.search(r'random', cc) or re.search(r'chrUn', cc) or re.search(r'alt', cc) or re.search(r'chrM', cc):
	#		continue
		chrom_seq = fastafile.fetch(cc).upper()
		sta = 0
		while sta < len(chrom_seq):
			aa = chrom_seq.find(dinucleotide, sta)
			if aa == -1:
				break
			di_pos.append(aa)
			sta=aa+1
			file_o.write(cc+"\t"+str(aa)+"\n")
	
	file_o.close()
	print "genome genome_reference  get %s sites %s" % (dinucleotide, len(di_pos))	
	fastafile.close()

	##find_dinucleotide(chrom_seq, dinucleotide)


	#file_o=open(file_out,'w')

	#sta = 0
	#di_pos = []
	#while sta < len(chrom_seq):
	#	aa = chrom_seq.find(dinucleotide, sta)
	#	if aa == -1:
	#		break
	#	di_pos.append(aa)
	#        sta=aa+1
	#	file_o.write(chrom+"\t"+str(aa)+"\n")
	#file_o.close()

	#print "genome genome_reference chrom %s get %s %s sites." % (chrom, dinucleotide, len(di_pos))





