import pysam
import sys
import os
import re
from collections import defaultdict
import pandas as pd
import argparse
import subprocess
from time import clock
import pandas as pd
import re
import signal
import gzip
import multiprocess
import time
import subprocess


'''
New flag: 
read old sta & end: ("SE", '%s,%s'%(old_sta, old_end)),
frag old sta & end: ("FG", '%s,%s'%(fragment_5_end, fragment_3_end)),
update by Bailey; 20200606

'''


def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

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
    
parser.add_argument("-r", "--reference",
                    dest = 'reference',
                    type=str,
                    required=True,
                    help="Reference DNA sequence in FASTA format, such as genome.fa")
    
parser.add_argument("-g","--Gap_len",
                    type=int,
                    default = 9,
                    help="The Tn5 gap length. default 9")

parser.add_argument("--CH_low",
                    default = 0.3,
                    help="The lower mCH ratio not bigger than CH_high. default 0.3")
parser.add_argument("--CH_high",
                        type=int,
                        default = 0.3,
                        help="The higher mCH ratio not lower than CH_high. default 0.3")
parser.add_argument("--min_CH_num",
                        type=int,
                        default = 4,
                        help="The min total CH number of both ends. default 4")
parser.add_argument("--min_mCH_num",
                        type=int,
                        default = 2,
                        help="The min mCH number if another end without CHs, Default is 2")
parser.add_argument("--min_BaseQuality",
                        type=int,
                        default = 25,
                        help="Min Base quality. default 25")
parser.add_argument("--min_MappingQuality",
                        type=int,
                        default = 30,
                        help="Min Map quality. default 30")
parser.add_argument("--Read_type",
                        type=str,
                        default = 'PE',
                        help="PE: pair end read, SE: single end read, default PE")
parser.add_argument("--thread",
                        type=int,
                        default = 10,
                        help="thread, default 1")
args = parser.parse_args(args)


input_bam = args.input_bam
output_bam = args.output_bam
reference = args.reference

Gap_len = args.Gap_len
CH_low = args.CH_low
CH_high = args.CH_high
min_CH_num = args.min_CH_num
min_mCH_num = args.min_mCH_num
minBQ = args.min_BaseQuality
minMQ = args.min_MappingQuality
 
Read_type = args.Read_type
threads = args.thread
        
fastafile = pysam.Fastafile(reference)


def get_MD_list (MDstring):
    '''
    Get the MD tag list. MD_tag: MD:Z:12G29G24
    MDstring: 12G29G24
    MD list: [[12,'G'], [29,'G']]
    '''
 #   print 'MDstring',MDstring
    pat = "[0-9]+[GCATN^]+"
    MD_s_list = re.findall(pat,MDstring)
    MD_list = [re.split(r'(\d+)',aa)[1:] for aa in MD_s_list]
    if re.search(r'[0-9]+$',MDstring):
        
        MD_final = re.split(r'([A-Z^]+)',re.findall(r'[GCATN^]*[0-9]+$',MDstring)[-1])[-1]
        MD_list.append([MD_final,''])
        
 #   print MD_list
    return MD_list

def formatSeqByCigar_I(seq, cigar, MD_list, pre_cut, off_cut):
    ''' Get cigarString.
    seq = query seq
    cigarString = '52M1I13M'
    cigar = [(52,M),(1,I),(13,M)]
    MD_list = [[12,'G'], [29,'G']]
    formatSeq = 'ATTTCCC' <which remove insert, add deletion, get mapped seq>
    '''
    cigar_index = {'M':0, 'I':1, 'D':2}
    formatSeq = ''
    pointer = 0; qstart = 0; qend = -1; origin_seq_len = 0
            
    new_cigar = cigar

    pre_cut_1 = pre_cut
    while pre_cut > 0:
        pair = new_cigar[0]
        operation = pair[0]
        cigar_len = pair[1]
        del new_cigar[0]
        if operation == 0 or operation == 1:
            if cigar_len >= pre_cut:
                cigar_len = cigar_len - pre_cut
                pre_cut = 0
                new_cigar.insert(0,(operation,cigar_len))   
            else:
                pre_cut = pre_cut - cigar_len
    
    off_cut_1 = off_cut  
    while off_cut > 0:
        pair = new_cigar[-1]
        operation = pair[0]
        cigar_len = int(pair[1])
        new_cigar.pop(-1)
        if operation == 0 or operation == 1:
            if cigar_len >= off_cut:
                cigar_len = cigar_len - off_cut
                off_cut = 0
                new_cigar.append((operation,cigar_len))
            else:
                off_cut = off_cut - cigar_len

    MD_list_new = MD_list   
  #  print MD_list_new
    while pre_cut_1 > 0:
        pair = MD_list_new[0]
        pre_base = int(pair[0])
        mis_base = pair[1]
        del MD_list_new[0]
        if pre_base >= pre_cut_1:
            pre_base = pre_base - pre_cut_1
            pre_cut_1 = 0
            MD_list_new.insert(0,[str(pre_base),mis_base])
        elif mis_base != '':
            pre_cut_1 = pre_cut_1 - pre_base
            if len(mis_base) >= pre_cut_1:
                cut_n = len(mis_base) - pre_cut_1
                mis_base = mis_base[-(len(mis_base) - pre_cut_1):]
                MD_list_new.insert(0,[str(0),mis_base])
                
                pre_cut_1 = 0
            else:
                pre_cut_1 = pre_cut_1 - len(mis_base)
            
 #   print MD_list_new
    while off_cut_1 > 0:
        pair = MD_list_new[-1]
        pre_base = int(pair[0])
        mis_base = pair[1]
        MD_list_new.pop(-1)
  #      print pair, off_cut_1

        if mis_base != '':
            if len(mis_base) >= off_cut_1:
                cut_n = len(mis_base) - off_cut_1
                mis_base = mis_base[:(len(mis_base) - off_cut_1)]
                MD_list_new.append([str(pre_base),mis_base])

                off_cut_1 = 0
            else:
                off_cut_1 = off_cut_1 - len(mis_base)
        if off_cut_1 > 0:
            if pre_base >= off_cut_1:
                pre_base = pre_base - off_cut_1
                off_cut_1 = 0
                MD_list_new.append([str(pre_base), ''])
            else:
                off_cut_1 = off_cut_1 - pre_base
  #      print MD_list_new
            
    MD_string = "".join(''.join(ss) for ss in MD_list_new)
    #print new_cigar, MD_string
    return new_cigar, MD_string

def formatSeqByCigar_I(seq, cigar, MD_list, pre_cut, off_cut):
    #
    cigar_index = {'M':0, 'I':1, 'D':2}
    formatSeq = ''
    pointer = 0; qstart = 0; qend = -1; origin_seq_len = 0

    new_cigar = cigar

    pre_cut_1 = pre_cut
    inerstion_pre = 0
    while pre_cut > 0:
        pair = new_cigar[0]
        operation = pair[0]
        cigar_len = pair[1]
        del new_cigar[0]
        if operation == 0:# or operation == 1:
            if cigar_len >= pre_cut:
                cigar_len = cigar_len - pre_cut
                pre_cut = 0
                new_cigar.insert(0,(operation,cigar_len))
            else:
                pre_cut = pre_cut - cigar_len
	elif operation == 1:
	    inerstion_pre += cigar_len
	    	
                
    off_cut_1 = off_cut
    inerstion_off = 0
    while off_cut > 0:
        pair = new_cigar[-1]
        operation = pair[0]
        cigar_len = int(pair[1])
        new_cigar.pop(-1)
        if operation == 0:# or operation == 1:
            if cigar_len >= off_cut:
                cigar_len = cigar_len - off_cut
                off_cut = 0
                new_cigar.append((operation,cigar_len))
            else:
                off_cut = off_cut - cigar_len
	elif operation == 1:
            inerstion_off += cigar_len

    MD_list_new = MD_list
  #  print MD_list_new
    while pre_cut_1 > 0:
        pair = MD_list_new[0]
        pre_base = int(pair[0])
        mis_base = pair[1]
        del MD_list_new[0]
        if pre_base >= pre_cut_1:
            pre_base = pre_base - pre_cut_1
            pre_cut_1 = 0
            MD_list_new.insert(0,[str(pre_base),mis_base])
        elif mis_base != '':
            pre_cut_1 = pre_cut_1 - pre_base
            if len(mis_base) >= pre_cut_1:
                cut_n = len(mis_base) - pre_cut_1
                mis_base = mis_base[-(len(mis_base) - pre_cut_1):]
                MD_list_new.insert(0,[str(0),mis_base])

                pre_cut_1 = 0
            else:
                pre_cut_1 = pre_cut_1 - len(mis_base)

    while off_cut_1 > 0:
        pair = MD_list_new[-1]
        pre_base = int(pair[0])
        mis_base = pair[1]
        MD_list_new.pop(-1)
  #      print pair, off_cut_1

        if mis_base != '':
            if len(mis_base) >= off_cut_1:
                cut_n = len(mis_base) - off_cut_1
                mis_base = mis_base[:(len(mis_base) - off_cut_1)]
                MD_list_new.append([str(pre_base),mis_base])

                off_cut_1 = 0
            else:
                off_cut_1 = off_cut_1 - len(mis_base)
        if off_cut_1 > 0:
            if pre_base >= off_cut_1:
                pre_base = pre_base - off_cut_1
                off_cut_1 = 0
                MD_list_new.append([str(pre_base), ''])
            else:
                off_cut_1 = off_cut_1 - pre_base
  #      print MD_list_new
    aa = 0
    for cc in new_cigar:   
        cigar_len = int(cc[1])
     #   print cigar_len,cc
        if cigar_len == 0:
            del new_cigar[aa]
        aa +=1

    MD_string = "".join(''.join(ss) for ss in MD_list_new)
    #print new_cigar, MD_string   
    #print new_cigar, MD_string
    return new_cigar, MD_string,inerstion_pre, inerstion_off               


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

def trim_read (read, chrom, clip_5_end, clip_3_end, trim_t, fragment_5_end, fragment_3_end, MapQ_flag):
        r1_sta = read.reference_start
        r1_end = read.reference_end
        r1_seq = read.get_reference_sequence
        r1_map_len = read.reference_length
        r1_query_len = read.query_length
        read_id = read.query_name
        old_sta = r1_sta
        old_end = r1_end

        if clip_5_end < 0:
            clip_5_end = 0
        if clip_3_end < 0:
            clip_3_end = 0
        
        pre_cut = 0
        off_cut = 0

        
        r1_end_gap = r1_query_len - read.query_alignment_end

        if read.query_alignment_start < clip_5_end:
            pre_cut = clip_5_end - read.query_alignment_start
            r1_sta = r1_sta + pre_cut
        if r1_end_gap < clip_3_end:
            off_cut = clip_3_end - r1_end_gap
            r1_end = r1_end - off_cut
       

        if r1_end < r1_sta:
            #print 'r1_end < r1_sta', r1_sta, r1_end
            return

#        print read_id,clip_5_end, clip_3_end, pre_cut, off_cut, r1_sta, r1_end, read.reference_start, read.reference_end
        
        read_MDstring = read.get_tag('MD')
        read_MD_list = get_MD_list(read_MDstring)
        cigarstring = read.cigar
        New_cigar, New_MD, inerstion_pre, inerstion_off = formatSeqByCigar_I(read.seq, cigarstring, read_MD_list, pre_cut, off_cut)

	pre_cut = pre_cut + inerstion_pre
	off_cut = off_cut + inerstion_off

        a = pysam.AlignedSegment()
        a.query_name = read.query_name + ":" + trim_t
        a.flag = read.flag
        a.reference_id = read.reference_id
        a.reference_start = r1_sta
        a.mapping_quality = read.mapping_quality
        a.cigar = (New_cigar)
        a.next_reference_id = read.next_reference_id
        if pre_cut == 0:
            if off_cut == 0:
                a.query_sequence = read.seq
                a.query_qualities = read.query_qualities
                XM_tag = read.get_tag('XM')
            else:
                a.query_sequence = read.seq[:-off_cut]
                a.query_qualities = read.query_qualities[:-off_cut]
                XM_tag = read.get_tag('XM')[:-off_cut]
        else:
            if off_cut == 0:
                a.query_sequence = read.seq[pre_cut:]
                a.query_qualities = read.query_qualities[pre_cut:]
                XM_tag = read.get_tag('XM')[pre_cut:]
            else:
                a.query_sequence = read.seq[pre_cut:-off_cut]
                a.query_qualities = read.query_qualities[pre_cut:-off_cut]
                XM_tag = read.get_tag('XM')[pre_cut:-off_cut]

        a.tags += (("MD", New_MD),
                  ('XG', read.get_tag('XG')),
                  ("NM", read.get_tag('NM')),
                  ("XM", XM_tag),
                  ("SE", '%s,%s'%(old_sta, old_end)),
                  ("FG", '%s,%s'%(fragment_5_end, fragment_3_end)),
                  ("MQ", MapQ_flag),
                  ("XR", read.get_tag('XR')))
        #a.set_tags(("MD", New_MD),('PG', read.get_tag('PG')),('XG', read.get_tag('XG')),("NM", read.get_tag('NM')),("XM", XM_tag),("XR", read.get_tag('XR')))

#        print "new v:",a.query_name,a

        return a

    

def PE_process(samfile,chrom, outbam):
    global lock
    lock = multiprocess.Lock()
    
    samfile = pysam.AlignmentFile(samfile, "rb")
    header = samfile.header
    outf = pysam.AlignmentFile(outbam, "wb", header=header)
    
    
   # aa = 0
    
    for read1, read2 in read_pair_generator(samfile,chrom):
        
        
        read_id = read1.query_name
       # print read_id
        
        r1_sta = read1.reference_start
        r2_sta = read2.reference_start
        r1_end = read1.reference_end
        r2_end = read2.reference_end
        r1_MQ = read1.mapping_quality
        r2_MQ = read2.mapping_quality

        
        r1_bq = read1.query_alignment_qualities
        r2_bq = read2.query_alignment_qualities
        
        read1_XM_tag = read1.get_tag('XM')
        read2_XM_tag = read2.get_tag('XM')
        read1_XG_tag = read1.get_tag('XG')
        read2_XG_tag = read2.get_tag('XG')
        
        if read1.is_reverse:
            fragment_5_end = r2_sta
            fragment_3_end = r1_end
     
            pre_MQ = r2_MQ
            off_MQ = r1_MQ
        else:
            fragment_5_end = r1_sta
            fragment_3_end = r2_end
            
            pre_MQ = r1_MQ
            off_MQ = r2_MQ

                
        MapQ_flag = str(pre_MQ) + "," + str(off_MQ)
        if fragment_3_end - fragment_5_end > Gap_len + Gap_len:  
            
            Gap_pre_len = Gap_len
            Gap_off_len = Gap_len
            
            if read1_XG_tag == 'CT' and read2_XG_tag == 'CT':
                trim_t = 'trim_off'
                Gap_pre_len = 0
                Gap_off_len = Gap_len
            elif read1_XG_tag == 'GA' and read2_XG_tag == 'GA':
                trim_t = 'trim_pre'
                Gap_pre_len = Gap_len
                Gap_off_len = 0
            else:
                trim_t = 'trim_both'
                
             
            
            r1_pre_tim = Gap_pre_len - (r1_sta - fragment_5_end)
            r1_off_tim = Gap_off_len - (fragment_3_end - r1_end)
            r2_pre_tim = Gap_pre_len - (r2_sta - fragment_5_end)
            r2_off_tim = Gap_off_len - (fragment_3_end - r2_end)
            
    #        print read_id, r1_pre_tim, r1_off_tim, r2_pre_tim, r2_off_tim
       #     print read_id, r1_sta, r1_end, r2_sta, r2_end
            
            r1_new = trim_read (read1, chrom, r1_pre_tim, r1_off_tim, trim_t, fragment_5_end, fragment_3_end, MapQ_flag)
            r2_new = trim_read (read2, chrom, r2_pre_tim, r2_off_tim, trim_t, fragment_5_end, fragment_3_end, MapQ_flag)

            r1_new.next_reference_start = r2_new.reference_start
            r2_new.next_reference_start = r1_new.reference_start

            
            lock.acquire()
            outf.write(r1_new)
            outf.write(r2_new)
            lock.release()
            
        else:     
          #  print 'Fragment length not enough, drop out'
            continue
    
    outf.close()
    samfile.close()
    return

def SE_process(samfile,chrom, outbam):
    global lock
    lock = multiprocess.Lock()
    
    samfile = pysam.AlignmentFile(samfile, "rb")
    header = samfile.header
    outf = pysam.AlignmentFile(outbam, "wb", header=header)
    
    #aa=0

    for read1 in samfile.fetch(chrom):
        
        read_id = read1.query_name

        r1_sta = read1.reference_start
        r1_end = read1.reference_end   

        r1_MQ = read1.mapping_quality
        r1_bq = read1.query_alignment_qualities        
        
        read1_XM_tag = read1.get_tag('XM')
        read1_XG_tag = read1.get_tag('XG')

        fragment_5_end = r1_sta
        fragment_3_end = r1_end
        
        pre_bq = r1_bq[:Gap_len]
        off_bq = r1_bq[-Gap_len:]
    
        MapQ_flag = 1
        if fragment_3_end - fragment_5_end > Gap_len + Gap_len:
            Gap_pre_len = Gap_len
            Gap_off_len = Gap_len
                
            if read1_XG_tag == 'CT':
                trim_t = 'trim_off'
                Gap_pre_len = 0
                Gap_off_len = Gap_len
            elif read1_XG_tag == 'GA':
                trim_t = 'trim_pre'
                Gap_pre_len = Gap_len
                Gap_off_len = 0
            else:
                trim_t = 'trim_both'
            
            r1_pre_tim = Gap_pre_len - (r1_sta - fragment_5_end)
            r1_off_tim = Gap_off_len - (fragment_3_end - r1_end)
             
      #      print r1_pre_tim, r1_off_tim
        
            r1_new = trim_read(read1, chrom, r1_pre_tim, r1_off_tim, trim_t, fragment_5_end, fragment_3_end, MapQ_flag)

            lock.acquire()
            outf.write(r1_new)
            lock.release()
            
        else:
     #       print 'Fragment length not enough, drop out'
            continue
        
    
    outf.close()
    samfile.close()
    return

        
def run_all():

    t0=clock()

    samfile = pysam.AlignmentFile(input_bam, "rb")
    assert samfile.check_index(),  "ErrorType: %s file does not have index file." % input_bam
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
        file_out_tmp = output_bam + '.' + chrom + ".bam"
        if Read_type == 'PE':
            
            result = pool.apply_async(PE_process, (input_bam, chrom, file_out_tmp, ))
        #    PE_process(input_bam, chrom, file_out_tmp)
        else:
            result = pool.apply_async(SE_process, (input_bam, chrom, file_out_tmp, ))
        
        results_q.append(result)
        tmp_ReadOut.append(file_out_tmp)
        
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
    print "Change Bam file by CHs spend %s" % str(t1-t0)
        
    Args_m = ['samtools merge -f -@ %s %s %s' % (threads, output_bam, " ".join(tmp_ReadOut))]
    subprocess.check_call(Args_m, shell=True)

    for my_file in tmp_ReadOut:
        if os.path.exists(my_file):
            os.remove(my_file)
    t2=clock()

    
    print "Final merged Bam file %s in %s" % (output_bam, str(t2-t1))
    
        
if __name__ == "__main__":

    run_all()
    
    
        
        


        
        

