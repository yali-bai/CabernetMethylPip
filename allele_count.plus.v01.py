import os
import sys
import re
import pandas as pd
import numpy as np
from collections import defaultdict
from collections import Counter
import matplotlib
from matplotlib.ticker import FuncFormatter

from time import clock
import multiprocess
import time
import subprocess

import pysam

import argparse
import signal
import gzip


args = sys.argv[1:]
parser = argparse.ArgumentParser(
    description="The code to call the modifications in the bam/sam reads. Update in 2019.05.10 by Bailey.",
    formatter_class = argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-i", "--input_bam",
                    required=True,
                    help="input bam format file containing sequencing reads.")

parser.add_argument("-o", "--output_bam",
                    type=str,
                    required=True,
                    help="The output BAM file.")

parser.add_argument("-a", "--output_allele",
                    type=str,
                    required=True,
                    help="The output allele file.")

parser.add_argument("-m", "--minPos_dis",
                    type=int,
                    default=1,
                    help="The min position distance value to define the same allele.")

parser.add_argument("-s", "--seq_type",
		    type=str,
		    default = "meta",
		    help="The sequence amplification type. (nextera | meta). meta means with barcode.")

parser.add_argument("--MapQ_flag",
		    action='store_true', 
		    default=True,
	   	    help="Mapping Quality flag, if it's Ture, low mapping quality reads will not consider in allele count with M start.")

parser.add_argument("--minMQ",
                    type=int, 
                    default=40,
                    help="Min Mapping Quality sum of read1 and read2, default 40.")

parser.add_argument("--thread",
                        type=int,
                        default = 10,
                        help="thread, default 1")
args = parser.parse_args(args)


bam_in = args.input_bam
outbam = args.output_bam
outall = args.output_allele
threads = args.thread
seq_type = args.seq_type
minPos_dis = 1
MapQ_flag=args.MapQ_flag
minMQ = args.minMQ

def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):

        qname = read.query_name
        #print qname
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]
            
def Count_allele(bam_in,chrom, outbam, outallele):
    global lock
    lock = multiprocess.Lock()


    samfile = pysam.AlignmentFile(bam_in, "rb")
    samfile_b = pysam.AlignmentFile(bam_in, "rb")

    header = samfile.header
    outf = pysam.AlignmentFile(outbam, "wb", header=header)
    
    Allele_index = 1
    
    
    read_dict = defaultdict()
    new_read = defaultdict()
    allele_count = []
    read_n = 0
    for read in samfile.fetch(chrom):
        readID = read.query_name
        query_name_cut =  re.split(':',read.query_name)
            
        if readID in read_dict:
            if not "flag" +str(read.flag) in read_dict[readID]:
                read.tags += [('AE',read_dict[readID][0])]
                read_dict[readID].append("flag" +str(read.flag))
                read.query_name = new_read[readID]

                lock.acquire()
                outf.write(read)
                lock.release()
            continue

        read_n +=1 
        #trim_t = query_name_cut[len(query_name_cut)-1]
        
        if MapQ_flag:
            r1MQ,r2MQ = re.split(",",read.get_tag('MQ'))
           
            if int(r1MQ)+int(r2MQ) < minMQ:
        #        trim_t = query_name_cut[len(query_name_cut)-1]
                
                #read.query_name = 'M'+str(Allele_index) +':'+chrom+':r'+str(read_n)+':'+trim_t+":"+read.get_tag('XG')+":"+read.get_tag('MQ')
		read.query_name = 'M'+str(Allele_index) +':'+chrom+':r'+str(read_n)+':'+read.get_tag('XG')+":"+read.get_tag('MQ')
                read.tags += [('AE',Allele_index)]
                read_dict[readID] = [Allele_index, "flag" +str(read.flag)]

                new_read[readID] = read.query_name

                lock.acquire()
                outf.write(read)
                lock.release()

		Allele_index += 1			
		continue


	if seq_type != "nextera":	
	        b1, b2 = readID.split(r':')[-2].split(r'_')
		meta_tag = [[b1, b2]]
	        Tags = [b1, b2]

        frag_sta, frag_end = map(int, read.get_tag('FG').split(','))
        chrom = read.reference_name
        read.tags += [('AE',Allele_index)]
        read_dict[readID] = [Allele_index, "flag" +str(read.flag)]
        
        #trim_t = query_name_cut[len(query_name_cut)-1]
        if MapQ_flag:
            #read.query_name = 'A'+str(Allele_index) +':'+chrom+':r'+str(read_n)+':'+trim_t+":"+read.get_tag('XG') + ":"+read.get_tag('MQ')
		read.query_name = 'A'+str(Allele_index) +':'+chrom+':r'+str(read_n)+':'+read.get_tag('XG') + ":"+read.get_tag('MQ')
        else:
            #read.query_name = 'A'+str(Allele_index) +':'+chrom+':r'+str(read_n)+':'+trim_t+":"+read.get_tag('XG')
		read.query_name = 'A'+str(Allele_index) +':'+chrom+':r'+str(read_n)+':'+read.get_tag('XG')

        new_read[readID] = read.query_name
        lock.acquire()
        outf.write(read)
        lock.release()
        
        half_tag = 0
        same_pos = 1
        #samfile_b = pysam.AlignmentFile(bam_in, "rb")
        for rr in samfile_b.fetch(chrom, frag_sta, frag_end):
            qname = rr.query_name
            if qname in read_dict:
                if not "flag" +str(rr.flag) in read_dict[qname]:
                    rr.tags += [('AE',read_dict[qname][0])]
                    read_dict[qname].append("flag" +str(rr.flag))
                    rr.query_name = new_read[qname]

                    lock.acquire()
                    outf.write(rr)
                    lock.release()
                continue

            query_name_cut =  re.split(':',rr.query_name)
            if MapQ_flag:
                r1MQ,r2MQ = re.split(",",rr.get_tag('MQ'))
           
                if int(r1MQ)+int(r2MQ) < minMQ:
                    #trim_t = query_name_cut[len(query_name_cut)-1]
                    #rr.query_name = 'M'+str(Allele_index) +':'+chrom+':r'+str(read_n)+':'+trim_t+":"+rr.get_tag('XG')+":"+rr.get_tag('MQ')
		    rr.query_name = 'M'+str(Allele_index) +':'+chrom+':r'+str(read_n)+':'+rr.get_tag('XG')+":"+rr.get_tag('MQ')
                    rr.tags += [('AE',Allele_index)]
                    read_dict[qname] = [Allele_index, "flag" +str(rr.flag)]

                    new_read[qname] = rr.query_name
                    lock.acquire()
                    outf.write(rr)
                    lock.release()
                    continue

            ss, ee = map(int, rr.get_tag('FG').split(','))
            dis_s = abs(frag_sta - ss)
            dis_e = abs(frag_end - ee)
            if dis_s <= minPos_dis or dis_e <= minPos_dis:
                read_dict[qname] = [Allele_index, "flag" +str(rr.flag)]
                rr.tags += [('AE',Allele_index)]
                same_pos += 1
                
                read_n +=1

                #trim_t = query_name_cut[len(query_name_cut)-1]
                if MapQ_flag:
                #    rr.query_name = 'A'+str(Allele_index) +':'+chrom+':r'+str(read_n)+':'+trim_t+":"+rr.get_tag('XG') + ":"+rr.get_tag('MQ')
			rr.query_name = 'A'+str(Allele_index) +':'+chrom+':r'+str(read_n)+':'+rr.get_tag('XG') + ":"+rr.get_tag('MQ')

                else:
                    #rr.query_name = 'A'+str(Allele_index) +':'+chrom+':r'+str(read_n)+':'+trim_t+":"+rr.get_tag('XG')
			rr.query_name = 'A'+str(Allele_index) +':'+chrom+':r'+str(read_n)+':'+rr.get_tag('XG')
                new_read[qname] = rr.query_name
                outf.write(rr)
                
		if seq_type != "nextera":
			b21, b22 = qname.split(r':')[-2].split(r'_')
                	meta_tag.append([b21, b22])
                	Tags.extend([b21, b22])      

	if seq_type != "nextera":
        	most_tags = Counter(Tags).most_common(2)
        	most_tag1 = most_tags[0][0]
        	
        	same_tag = 0 
        	if len(most_tags) == 1:
        	    allele_count.append([chrom, str(frag_sta), str(frag_end), str(same_pos), str(same_tag), str(half_tag), format(1.0, '.2f'), format(1.0, '.2f'), most_tag1+','+most_tag1])
        	    Allele_index += 1 
        	    continue


        	most_tag2 = most_tags[1][0]
        	for b1, b2 in meta_tag:
        	    if (most_tag1 == b1 and most_tag2 == b2) or (most_tag1 == b2 and most_tag2 == b1):
        	        same_tag += 1
        	    elif most_tag1 == b1 or most_tag2 == b2 or most_tag1 == b2 or most_tag2 == b1:
        	        half_tag += 1


        	allele_count.append([chrom, str(frag_sta), str(frag_end), str(same_pos), str(same_tag), str(half_tag), format(1.0*same_tag/same_pos, '.2f'), format(1.0*(same_tag+half_tag)/same_pos, '.2f'), most_tag1+','+most_tag2])
		allele_count = pd.DataFrame(allele_count)
	        allele_count.to_csv(outallele, header = None, index=None,sep='\t')


        Allele_index += 1 

    
    outf.close()
    samfile.close()
    samfile_b.close()
    return


def run_all():

    t0=clock()

    samfile = pysam.AlignmentFile(bam_in, "rb")
    assert samfile.check_index(),  "ErrorType: %s file does not have index file." % bam_in
    print 'Start Analysis threads: %s' % threads


    header = pd.DataFrame(samfile.header['SQ'])
    Chrs = header['SN']
    samfile.close()

    pool = multiprocess.Pool(threads,init_worker)
    results_q = []
    tmp_ReadOut = []
    for chrom in Chrs:
        if re.search(r'_',chrom):
            continue
        if chrom == "chrM":
            continue

        print chrom
        file_out_tmp = outbam + '.' + chrom + ".bam"
        allele_out_t = outall + '.' + chrom + ".allele.txt"
        
        result = pool.apply_async(Count_allele, (bam_in, chrom, file_out_tmp, allele_out_t ))
	#Count_allele(bam_in, chrom, file_out_tmp, allele_out_t)
	#results_q.append(result)
        tmp_ReadOut.extend([file_out_tmp, allele_out_t])
        
    try:
        print "Waiting 1 seconds"
        time.sleep(5)
    except KeyboardInterrupt:
        print "Caught KeyboardInterrupt, terminating workers"
        pool.terminate()
        pool.join()
    else:
        print "Quitting normally"
        pool.close()
        pool.join()
        print 'Pool finished'

    for result in results_q:
        print(result.get())
        
    t1=clock()
    print "Bam file count alleles spend %s" % str(t1-t0)

    Args_m = ['samtools merge -f -@ %s %s %s.chr*.bam ' % (threads, outbam, outbam)]
    print Args_m
    subprocess.check_call(Args_m, shell=True)

    if seq_type != "nextera":
    	Args_m = ['cat %s.chr*.allele.txt > %s.allele.txt ' % (outall, outall)]
    	print Args_m
    	subprocess.check_call(Args_m, shell=True)

    for my_file in tmp_ReadOut:
        if os.path.exists(my_file):
            os.remove(my_file)
    t2=clock()

    print "Final count allele to %s in %s" % (outbam, str(t2-t1))


if __name__ == "__main__":

    run_all()



