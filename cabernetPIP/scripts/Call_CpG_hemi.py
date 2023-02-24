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
from numpy import  *
import argparse
import signal
import gzip
import shutil
import operator
from functools import reduce
import timeit

def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

#```
#input :
#chr1 3001629 0 A4 CT 24,24
#chr1 3001629 0 A4 GA 42,42
#chr1 3003226 0 A5 CT 24,24
#chr1 3003339 0 A6 CT 42,42
#chr1 3003339 1 A6 GA 42,42
#chr1 3003379 0 A6 CT 42,42
#chr1 3003379 1 A6 GA 42,42
#```


CpG_context_f = '/share/home/baiyl/merlot/TAPS-merlot/workflow/Analysis/allele_dropout/metaEM_ado/SRC_bin/data/bismark_mm10/CpG_context_TE_2.allele_count.XG.3.sortn.chr1.txt'

CpG_context_f = sys.argv[1]


args = sys.argv[1:]
parser = argparse.ArgumentParser(
    description="The code to call the modifications in the bam/sam reads. Update in 2019.05.10 by Bailey.",
    formatter_class = argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-i", "--input",
                    required=True,
                    help="input CpGs context file of allele methylation reads.")

parser.add_argument("-o", "--outdir",
                    type=str,
                    required=True,
                    help="The output dir .")

parser.add_argument("--sampleID",
                    type=str,
                    required=True,
                    help="The sampleID to output.")

parser.add_argument("--species",
                    type=str,
                    default = 'mm10',
                    help="The species of genome.such as mm10, hg38")

parser.add_argument("--MapQ_flag",
                    action='store_true', 
                    default=True,
                    help="Mapping Quality flag, if it's Ture, low mapping quality reads will not consider in allele count with M start.")

parser.add_argument("--threads",
                        type=int,
                        default = 2,
                        help="thread, default 2")

args = parser.parse_args(args)
CpG_context_f = args.input
outdir = args.outdir
sampleID = args.sampleID
species = args.species
threads = args.threads
MapQ_flag = args.MapQ_flag

bedtools='/share/newdata4/baiyl/basic/tools/bedtools2/bin/bedtools'


if species == 'hg38':
    chrs = range(1,23)
elif species == 'mm10':
    chrs = range(1,20)
    chrs.extend(['X','Y'])
    
def POS_Allele_process(chrom, pos, pos_all, outf, outf_sameAllele, outf_posMerge):
    Pos_Methyl = []
    pos_all[0] = pos_all[0].astype('int')
    Alleles = [mm for mm in set(pos_all[1])]
    strand_news = {'CT':"+",'GA':"-"}
    Alleles_methyl = []
    depth=0
    for allele in Alleles:
        strands = list(set(pos_all.loc[pos_all[1]==allele, 2]))
        methyls_dict = defaultdict()

        allele_m = []
        for strand in strands:
            strand_new = strand_news[strand]
            
            same_all = pos_all.loc[(pos_all[1]==allele)&(pos_all[2]==strand),]
            depth_all = same_all.shape[0]
            mean_methyl = mean(same_all.iloc[:,0])
            mean_methyl_01 = int(bool(mean_methyl-0.5>0))
            
            support_ratio = 1.0*same_all.loc[same_all[0]==mean_methyl_01,].shape[0]/depth_all
            MapQ_all = mean(map(int,reduce(operator.add,same_all.loc[:,3:4].values.tolist())))

            Allele_results = '\t'.join([chrom, str(pos), str(len(Alleles)), allele, strand_new, str(mean_methyl), str(mean_methyl_01), str(support_ratio),str(depth_all), str(MapQ_all)])+'\n'
            outf.write(Allele_results)
            
            methyls_dict[strand_new] =(mean_methyl_01)
            allele_m.append(mean_methyl_01)
	    depth += 1
            

        Pos_Methyl.append(mean(allele_m))
        if len(methyls_dict) == 2:
            Allele_results = '\t'.join([chrom, str(pos), str(len(Alleles)), allele, str(mean(methyls_dict.values())), str(methyls_dict['+']), str(methyls_dict['-'])])+'\n'
            outf_sameAllele.write(Allele_results)
    
    Pos_results = '\t'.join([chrom, str(pos), str(len(Alleles)), str(mean(Pos_Methyl)), ",".join(map(str,Pos_Methyl)),str(depth)])+'\n'
    outf_posMerge.write(Pos_results)
        
#%timeit POS_Allele_process
    

def call_Methl_values(CpG_context_f, chrom, outfile):
    outf = open(outfile, 'w')
    outf_sameAllele = open(outfile+'.sameAllele.txt','w')
    outf_posMerge = open(outfile+'.posMerge.txt','w')
    global lock
    lock = multiprocess.Lock()
    
#input :
#    CpG_methyl = pd.read_table(CpG_context_f,header=0,index_col=None,sep=' ')
    CpG_methyl = pd.read_table(CpG_context_f,header=None,index_col=None,sep='\t')
    CpG_methyl.columns = ['chrom','Pos','Pos2','Reads']
    #aa = 0
    t1 = clock()
    for mm in CpG_methyl.index:
     #   aa +=1
     #   if aa>10:
     #       break
        pos=CpG_methyl.loc[mm,'Pos']
        Pos_alleles = pd.DataFrame([re.split(',',cpgs) for cpgs in re.split(':',CpG_methyl.loc[mm,'Reads'])])
        if Pos_alleles.shape[0]==0:
            continue
        POS_Allele_process(chrom, pos, Pos_alleles, outf, outf_sameAllele, outf_posMerge)

    #    print "time: %s" % str(clock()-t1)
        t1 = clock()
        CpG_methyl.drop(index=CpG_methyl.index)

        
    outf_posMerge.close()
    outf_sameAllele.close()     
    outf.close()
    
    
if __name__ == '__main__':
    t0=clock()
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if not os.path.exists("%s/%s_split" % (outdir, sampleID)):
        os.makedirs("%s/%s_split" % (outdir, sampleID))
    print time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))

    pool = multiprocess.Pool(threads,init_worker)
    results_q = []
    tmp_ReadOut = []
    t1=clock()

    lock = multiprocess.Lock()

    print 'Split CpG file'
    print time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))
    for cc in chrs:
        chrom = 'chr' + str(cc) 
	lock.acquire()
        Args_m = ["grep -w %s %s |awk 'BEGIN{pos=\"Pos\";reads=\"Reads\"}{if($2==pos){reads=reads\":\"$3\",\"$4\",\"$5\",\"$6;}else{if(pos!=\"Pos\"){print $1,pos,pos,reads};pos=$2;reads=$3\",\"$4\",\"$5\",\"$6;}}'|sed 's/ /\t/g'> %s/%s_split/%s.%s.CpG_context.tmp" % (chrom, CpG_context_f, outdir, sampleID, sampleID, chrom)]
        subprocess.check_call(Args_m, shell=True)
	lock.release()
        CpG_context_tmp = '%s/%s_split/%s.%s.CpG_context.tmp' % (outdir, sampleID, sampleID, chrom)
        
        chr_size = os.path.getsize(CpG_context_tmp)
        print chrom, chr_size
        split_num = int(chr_size / 10000000+1)

	lock.acquire()
        Args_m = ['%s split -n %s -i %s -p %s' %
              (bedtools, split_num, CpG_context_tmp, "%s/%s_split/%s.%s" % (outdir,sampleID, sampleID, chrom))]
        subprocess.check_call(Args_m, shell=True)
	lock.release()

    for cc in chrs:
	chrom = 'chr' + str(cc)
	CpG_context_tmp = '%s/%s_split/%s.%s.CpG_context.tmp' % (outdir, sampleID, sampleID, chrom)
	chr_size = os.path.getsize(CpG_context_tmp)
        split_num = int(chr_size / 10000000+1)

        for part in range(1,split_num+1):
            CpG_context_tmp_2 = '%s/%s_split/%s.%s.%s.bed'%(outdir,sampleID, sampleID, chrom, "{:0>5}".format(part))
            out_basename = '%s/%s_split/%s.%s.%s.CpG_allele_m.tmp' % (outdir,sampleID, sampleID, chrom, "{:0>5}".format(part))

            result = pool.apply_async(call_Methl_values, (CpG_context_tmp_2, chrom, out_basename ))
            results_q.append(result)
	    #call_Methl_values(CpG_context_tmp_2, chrom, out_basename)
        
        
    print time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))

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

    #for result in results_q:
    #    print(result.get())

    print time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))
    
    
    Args_m = ['cat %s/%s_split/%s.chr*.CpG_allele_m.tmp |gzip > %s/%s.CpG_allele_m.txt.gz' % (outdir, sampleID, sampleID, outdir, sampleID)]
    subprocess.check_call(Args_m, shell=True)
    Args_m = ['cat %s/%s_split/%s.chr*.CpG_allele_m.tmp.sameAllele.txt |gzip > %s/%s.CpG_allele_m.sameAllele.txt.gz' % (outdir, sampleID, sampleID, outdir, sampleID)]
    subprocess.check_call(Args_m, shell=True)
    Args_m = ['cat %s/%s_split/%s.chr*.CpG_allele_m.tmp.posMerge.txt |gzip > %s/%s.CpG_allele_m.posMerge.txt.gz' % (outdir, sampleID, sampleID, outdir, sampleID)]
    subprocess.check_call(Args_m, shell=True)

#    outdir_tmp = "%s/%s_split" % (outdir, sampleID)
#    if os.path.exists(outdir_tmp):
#        shutil.rmtree(outdir_tmp)
#        
    print time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))    
    print "Allele cpgs call finished!!!"
    

