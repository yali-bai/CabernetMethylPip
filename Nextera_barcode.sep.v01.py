import os
import sys
import pandas as pd
import gzip
import re
import gzip
import argparse
from collections import defaultdict
from collections import Counter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import Levenshtein
import difflib

import math
import numpy as np
import time

from Bio import SeqIO
import gzip
from time import clock
import gc
import time
import subprocess
from time import clock
import shutil
import multiprocess
import signal


'''
Mark Meta tags from read1 & read2.
"@readid:metaA-metaB "

The fragment is: meta16_Tn5_fragment_Tn5_meta16
Write by Bailey in 20190810.
'''

parser = argparse.ArgumentParser(description="The code is used to seperate reads by different cell barcode. Update in 2019.08.10 by Bailey. \nCopyright (c) Bailey")

parser.add_argument("-R1", "--Read1_fastq", dest = 'Read1_fastq',
                    help="Read1 fastq file path.")
parser.add_argument("-R2", "--Read2_fastq",  dest = 'Read2_fastq',
                    help="Read2 fastq file path.")

parser.add_argument("-T5","--T5_tage",  dest = 'T5_tag_file',
                    help="T5_tag_file .")

parser.add_argument("-T7","--T7_tage",  dest = 'T7_tag_file',
                    help="T7_tag_file .")

parser.add_argument("-o", "--outdir", dest = 'outdir',
                    help="The outdir of output files")

parser.add_argument("--min_dis",
                    type = int,
                    default=2,
                    help="min dist of (index+barcode seq) and seq (Default 2).")
parser.add_argument("--Tn5",
                    type = str,
                    default='AGATGTGTATAAGAGACAG',
                    help="The Tn5 seq after barcode(Default AGATGTGTATAAGAGACAG)..")

parser.add_argument("--prefix",
                    type = str,
		    help="The out name ")


parser.add_argument("--seprate",
                    action="store_true",
                    help="Seprate reads by barcode (Default False).")
parser.add_argument("--threads",
                    type = int,
                    default=2,
                    help="multiprocess (Default 2).")

parser.add_argument("--bysize",
		    action="store_false",
		    help="split sequences by size (Default True)")

parser.add_argument("--bypart",
		    dest = 'bypart',
		    action="store_true",
                    help="split sequences by parts (Default False)")

parser.add_argument("--split_size",
		type = int,
		default=1000000,
		help="split sequences into multi parts with N sequences")

parser.add_argument("--split_part",
                type = int,
                default=10,
		help="split sequences into N parts (Default 10)")


args = parser.parse_args()

Read1_fastq = args.Read1_fastq
Read2_fastq = args.Read2_fastq
T5_tag_file = args.T5_tag_file
T7_tag_file = args.T7_tag_file
outdir = args.outdir
bypart = args.bypart
bysize = args.bysize
split_size = args.split_size
split_part = args.split_part

if bypart:
	split_a = "p"
	split_s = split_part
else:
	split_a = "s"
	split_s = split_size


prefix = args.prefix
#if not args.prefix:
#	basename =re.split(r'[._]R1', os.path.basename(Read1_fastq))[0]
min_dis = args.min_dis
Tn5 = args.Tn5
seprate = args.seprate
threads = args.threads

print "seprate: %s " % seprate

def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def get_meta_barcode_seq(file_in):
    '''barcode file:
        "bar_1S  GATATG"
    or meta tags file:
        "META16-1        GGCACCGAAAA"
    '''
    barcode_seq = {}
    bac = pd.read_table(file_in, header=None)
    barcode_seq = bac.to_dict()[0]
    #barcode_seq = bac.set_index(0).to_dict()[1]
    return barcode_seq

def choose_barcode(seq_record, meta_tags):
    read_seq = seq_record.seq
    
    meta_dd = ''
    tmp_meta_score = 9999
    tmp_meta_seq = ''
    for meta in meta_tags.keys():
        metatag = meta_tags[meta]
        query_cut_seq = str(read_seq[:len(metatag)])
        #hamming_dis = Levenshtein.hamming(query_cut_seq, metatag)
	#seq_dis = Levenshtein.distance(query_cut_seq, metatag)
	seq_dis = Levenshtein.hamming(query_cut_seq, metatag)
        if seq_dis < tmp_meta_score:
            tmp_meta_score = seq_dis
            meta_dd = meta
            tmp_meta_seq = metatag
    
    tn5_sta = 0
    tn5_dd = ''
    if read_seq[len(tmp_meta_seq):len(tmp_meta_seq)+len(Tn5)] < 6:
        tn5_sta = len(tmp_meta_seq)
    else:
        tn5_mindis = 9999
        for ss in range(3, len(read_seq)-19-10):
            cut_seq = read_seq[ss:ss+len(Tn5)]
            tn5_dis = Levenshtein.hamming(str(cut_seq), Tn5)
	    #tn5_dis = Levenshtein.distance(str(cut_seq), Tn5)
	   
            if tn5_dis < tn5_mindis:
                tn5_mindis = tn5_dis
                tn5_sta = ss

    if tmp_meta_score > min_dis+5:
        meta_dd = "Puzzled"
    
    new_seq = SeqRecord(read_seq[tn5_sta+len(Tn5):], id=seq_record.id, description =seq_record.description, name=seq_record.name)
    new_seq.letter_annotations["phred_quality"] = seq_record.letter_annotations['phred_quality'][tn5_sta+len(Tn5):]
    
    return meta_dd, new_seq

def choose_barcode2(seq_record, meta_tags):
    read_seq = seq_record.seq

    meta_dd = ''
    tmp_meta_score = 9999
    tmp_meta_seq = ''

    tn5_sta = 0
    tn5_dd = ''
    tn5_mindis = 9999
    for ss in range(3, len(Tn5)):
	cut_seq = read_seq[ss:ss+len(Tn5)]
	tn5_dis = Levenshtein.hamming(str(cut_seq), Tn5)
	if tn5_dis < tn5_mindis:
		tn5_mindis = tn5_dis
		tn5_sta = ss

    query_cut_seq = str(read_seq[:ss])
    for meta in meta_tags.keys():
	metatag = meta_tags[meta]
	#if (len(metatag)!=tn5_sta):
	#	continue
	#seq_dis = Levenshtein.hamming(query_cut_seq, metatag)
	seq_dis = Levenshtein.distance(query_cut_seq, metatag)
        if seq_dis < tmp_meta_score:
            tmp_meta_score = seq_dis
            meta_dd = meta
            tmp_meta_seq = metatag

    if tn5_mindis > min_dis+3:
	meta_dd = "Puzzled"
    #if tmp_meta_score > min_dis:
    #    meta_dd = "Puzzled"

    new_seq = SeqRecord(read_seq[tn5_sta+len(Tn5):], id=seq_record.id, description =seq_record.description, name=seq_record.name)
    new_seq.letter_annotations["phred_quality"] = seq_record.letter_annotations['phred_quality'][tn5_sta+len(Tn5):]

    return meta_dd, new_seq


def main(Read1_fastq, Read2_fastq, outdir, out_basename):
    global lock
    lock = multiprocess.Lock()

    out_meta_com = []
    read1_meta = defaultdict(dict)
    handle1 = gzip.open(Read1_fastq, "r")
    if not seprate:
	    print "Not Sep"
	    out1 = gzip.open('%s/%s.R1.fastq.gz' % (outdir, out_basename), "wt")
   	    out2 = gzip.open('%s/%s.R2.fastq.gz' % (outdir, out_basename), "wt") 

    handle2 = gzip.open(Read2_fastq, "r")
    records2 = SeqIO.parse(handle2, "fastq")
    

    meta_all = set()
    for seq_record1 in SeqIO.parse(handle1, "fastq"):
        seq_record2 = records2.next()
        old_id = seq_record1.id
        T5_tag, new_seq1 = choose_barcode(seq_record1, T5_tags)
        T7_tag, new_seq2 = choose_barcode(seq_record2, T7_tags)
	T5_tag = str(T5_tag)
	T7_tag = str(T7_tag)

	#print T5_tag,T7_tag
	if re.search(r'Puzz', T5_tag) or re.search(r'Puzz', T7_tag):
                continue

	merged_id = "%s:T5_%s.T7_%s" % (seq_record1.id, T5_tag, T7_tag)
        new_seq1.id = merged_id
        new_seq2.id = merged_id
	meta_all.add("%s_%s" % (T5_tag, T7_tag))

	#if len(new_seq1.seq) < 20:
	#	continue

	#if re.search(r'Puzz', merged_id):
	#	continue
        
        seq1_quality = ''.join([chr(cc+33) for cc in new_seq1.letter_annotations['phred_quality']])
        seq2_quality = ''.join([chr(cc+33) for cc in new_seq2.letter_annotations['phred_quality']])

	if seprate:
		out1 = gzip.open('%s/%s.%s_%s.R1.fastq.gz' % (outdir, out_basename, T5_tag, T7_tag), "a")
		out2 = gzip.open('%s/%s.%s_%s.R2.fastq.gz' % (outdir, out_basename, T5_tag, T7_tag), "a")
		out_meta_com.append('%s/%s.%s_%s' % (outdir, out_basename, T5_tag, T7_tag))
        
	lock.acquire()
        out1.write('@%s\n'%merged_id)        
        out1.write(str(new_seq1.seq)+"\n")
        out1.write('+\n')
        out1.write(seq1_quality+"\n")
	lock.release()
	        

	lock.acquire()
        out2.write('@%s\n'%merged_id)
        out2.write(str(new_seq2.seq)+"\n")
        out2.write('+\n')
        out2.write(seq2_quality+"\n")
	lock.release()

	if seprate:
		out1.close()
		out2.close()

    
    handle1.close()
    handle2.close()

    if not seprate:
    	out1.close()
    	out2.close()
	
    return meta_all
	

    
    
if __name__ == '__main__':
    t0=clock()
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if not os.path.exists("%s/%s_split" % (outdir, prefix)):
        os.makedirs("%s/%s_split" % (outdir, prefix))
    if not os.path.exists("%s/%s_sep" % (outdir, prefix)):
        os.makedirs("%s/%s_sep" % (outdir, prefix))    
    
    r1_basename = os.path.basename(Read1_fastq).split(r'.fq.gz')[0].split(r'.fastq.gz')[0]
    r2_basename = os.path.basename(Read2_fastq).split(r'.fq.gz')[0].split(r'.fastq.gz')[0]
    raw_basename = os.path.basename(Read1_fastq).split('.R1_val')[0].split('.R1.fastq.gz')[0]
    
    T5_tags = get_meta_barcode_seq(T5_tag_file)
    T7_tags = get_meta_barcode_seq(T7_tag_file)
    
    #Args_m = ['/share/home/baiyl/basic/tools/seqkit split2 -%s %s -1 %s -2 %s -O %s -f' % 
    #      (split_a, split_s, Read1_fastq, Read2_fastq, "%s/%s_split" % (outdir,prefix))]
    #subprocess.check_call(Args_m, shell=True)
    
    lock = multiprocess.Lock()
    pool = multiprocess.Pool(threads,init_worker)
    results_q = [] 
    
    
    R1_files = []
    R2_files = []
    
    files = os.listdir("%s/%s_split" % (outdir,prefix))
    files_n = len(files) / 2
    for part in range(1,files_n+1):
        read1 = "%s/%s_split/%s.part_%s.fq.gz" % (outdir,prefix, r1_basename, "{:0>3}".format(part))
        read2 = "%s/%s_split/%s.part_%s.fq.gz" % (outdir,prefix, r2_basename, "{:0>3}".format(part))
        
        out_basename = '%s.part_%s'%(prefix, "{:0>3}".format(part))
        print part, read1, read2
        print out_basename
        
        #main(read1, read2, "%s/%s_sep" % (outdir, prefix), '%s.part_%s'%(raw_basename, "{:0>3}".format(part)))
        result = pool.apply_async(main, (read1, read2, "%s/%s_sep" % (outdir, prefix), out_basename, ))
	#main(read1, read2, "%s/%s_sep" % (outdir, prefix), out_basename)
	#main(read1, read2, "%s/%s_sep" % (outdir, prefix), out_basename)
        results_q.append(result)
    
        r1_tmp = '%s/%s_sep/%s.R1.f*q.gz' % (outdir, prefix, out_basename)
        r2_tmp = '%s/%s_sep/%s.R2.f*q.gz' % (outdir, prefix, out_basename)
	if not seprate:
        	R1_files.append(r1_tmp)
        	R2_files.append(r2_tmp)

    try:
        print "Waiting 1 seconds"
        time.sleep(1)
    except KeyboardInterrupt:
        print "Caught KeyboardInterrupt, terminating workers"
        pool.terminate()
        pool.join()
    else:
        print "Quitting normally"
        pool.close()
        pool.join()
        print 'Pool finished'
    metas = set()
    for result in results_q:
        print(result.get())
        out_meta = result.get()
        metas = metas | out_meta

    if seprate:
	for mm in metas:
		Args_m = ['zcat %s/%s_sep/%s.part*.%s.R1.f*q.gz | gzip > %s/%s.%s.R1.fastq.gz' % (outdir, prefix, prefix, mm, outdir, prefix, mm)]
		subprocess.check_call(Args_m, shell=True)
		Args_m = ['zcat %s/%s_sep/%s.part*.%s.R2.f*q.gz | gzip > %s/%s.%s.R2.fastq.gz' % (outdir, prefix, prefix, mm, outdir, prefix, mm)]
		subprocess.check_call(Args_m, shell=True)
		 
    else:
	    Args_m = ['zcat %s |gzip > %s/%s.R1.fastq.gz' % (' '.join(R1_files), outdir, prefix)]
	    subprocess.check_call(Args_m, shell=True)
	
	    Args_m = ['zcat %s |gzip > %s/%s.R2.fastq.gz' % (' '.join(R2_files), outdir, prefix)]
	    subprocess.check_call(Args_m, shell=True)

    outdir1 = "%s/%s_split" % (outdir, prefix)
    outdir2 = "%s/%s_sep" % (outdir, prefix)
    if os.path.exists(outdir1):
        shutil.rmtree(outdir1)
    
    if os.path.exists(outdir2):
        shutil.rmtree(outdir2)
    
    t1 = clock()
    print 'Time spend all: %s' % str(t1-t0)
    
    
