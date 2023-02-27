import sys
import os
import re
from collections import defaultdict
import pandas as pd
import operator
import multiprocess
import subprocess
import argparse
import signal
from time import clock
import time
import shutil

def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

args = sys.argv[1:]
parser = argparse.ArgumentParser(
    description="The code to call the modifications in the bam/sam reads. Update in 2019.05.10 by Bailey.",
    formatter_class = argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--phased_read_file",
                    required=True,
                    help="input phased read file.")

parser.add_argument("--cpg_context_file",
                    required=True,
                    help="input CpG context file.")

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

parser.add_argument("--threads",
                        type=int,
                        default = 1,
                        help="thread, default 1")

args = parser.parse_args(args)
phased_read_file = args.phased_read_file
cpg_context_file = args.cpg_context_file
outdir = args.outdir
sampleID = args.sampleID
species = args.species
threads = args.threads

    
if species == 'hg38':
    chrs = range(1,23)
elif species == 'mm10':
    chrs = range(1,20)
    chrs.extend(['X','Y'])

def sep_parternal (read_phased_f, CpG_context_f, CpG_sep_out):
    read_phased = pd.read_table(read_phased_f,header =None,index_col=None)
    CpG_context = pd.read_table(CpG_context_f,header=0,index_col =0,sep=' ')
    
    read_phased = read_phased[0].tolist()
    CpG_context_read = CpG_context.index.tolist()
    CpG_readID = set(read_phased).intersection(set(CpG_context_read))


    CpG_context_cut = CpG_context.loc[list(CpG_readID),]
    CpG_context_cut.to_csv(CpG_sep_out, header=True,index=True,sep='\t')
    return
    
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

    print 'Split CpG file'
    print time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))
    for cc in chrs:
        chrom = 'chr' + str(cc)   
        
        Args_m = ["grep -w %s %s|sort|awk 'BEGIN{readid=\"readID\";alleles=\"Alleles\"}{if($1==readid){alleles=alleles\":\"$3\",\"$4\",\"$5}else{print readid,alleles;readid=$1;alleles=$3\",\"$4\",\"$5}}' > %s/%s_split/%s.%s.CpG_context.tmp" % (chrom, cpg_context_file, outdir, sampleID, sampleID, chrom)]
        subprocess.check_call(Args_m, shell=True)
        CpG_context_f = '%s/%s_split/%s.%s.CpG_context.tmp' % (outdir, sampleID, sampleID, chrom)
        
        Args_m = ["less %s|grep :%s:|awk '$2~\"M\"'> %s/%s_split/%s.%s.read_phased_M.tmp" % (phased_read_file, chrom, outdir, sampleID, sampleID, chrom)]
        subprocess.check_call(Args_m, shell=True)
        read_Maternal = '%s/%s_split/%s.%s.read_phased_M.tmp' % (outdir, sampleID, sampleID, chrom)
        
        Args_m = ["less %s|grep :%s:|awk '$2~\"P\"'> %s/%s_split/%s.%s.read_phased_P.tmp" % (phased_read_file, chrom, outdir, sampleID, sampleID, chrom)]
        subprocess.check_call(Args_m, shell=True)
        read_Paternal = '%s/%s_split/%s.%s.read_phased_P.tmp' % (outdir, sampleID, sampleID, chrom)
        
        out_Maternal = '%s/%s_split/%s.%s.cpgs_phased_M.tmp' % (outdir, sampleID, sampleID, chrom)
        out_Paternal = '%s/%s_split/%s.%s.cpgs_phased_P.tmp' % (outdir, sampleID, sampleID, chrom)
        
        result = pool.apply_async(sep_parternal, (read_Maternal, CpG_context_f, out_Maternal ))
        results_q.append(result)
        result = pool.apply_async(sep_parternal, (read_Paternal, CpG_context_f, out_Paternal ))
        results_q.append(result)

#        sep_parternal(read_Maternal, CpG_context_f, out_Maternal )
    

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

    Args_m = ['cat %s/%s_split/%s.chr*.cpgs_phased_M.tmp > %s/%s.cpgs_phased_M.txt' % (outdir, sampleID, sampleID, outdir, sampleID)]
    subprocess.check_call(Args_m, shell=True)
    Args_m = ['cat %s/%s_split/%s.chr*.cpgs_phased_P.tmp > %s/%s.cpgs_phased_P.txt' % (outdir, sampleID, sampleID, outdir, sampleID)]
    subprocess.check_call(Args_m, shell=True)
    

    outdir_tmp = "%s/%s_split" % (outdir, sampleID)
    if os.path.exists(outdir_tmp):
        shutil.rmtree(outdir_tmp)

    print time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))
    print "Reads Allele separate finished!!!"
        
        
        


