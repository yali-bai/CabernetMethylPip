import pysam
from collections import defaultdict
import pandas as pd
import sys
import os
import re
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


BamFile = "../../../../Data_run/yuan_20200922_nano_scEM_5mC_5hmC/Output_s1/bismark_mm10/D200914_hmC_L1C_5A_S1.mm10.bam"
SnpFile = "data/C57BL_6NJ.DBA_1J.snp.bed"

#Bam_File = argv[1]
#Snp_File = argv[2]

def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

    
args = sys.argv[1:]
parser = argparse.ArgumentParser(
    description="The code to call the modifications in the bam/sam reads. Update in 2019.05.10 by Bailey.",
    formatter_class = argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-b", "--inbam",
                    required=True,
                    help="input Bam file.")

parser.add_argument("-s", "--snp",
                    required=True,
                    help="input snp bed file of alleles.")

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
Bam_File = args.inbam
outdir = args.outdir
sampleID = args.sampleID
species = args.species
threads = args.threads
Snp_File = args.snp


if species == 'hg38':
    chrs = range(1,23)
elif species == 'mm10':
    chrs = range(1,20)
    chrs.extend(['X','Y'])


def Get_alleles(SnpFile, BamFile, outfile):
    snpPos2 = pd.read_table(SnpFile,header=None,index_col =None)
    samfile = pysam.Samfile(BamFile, "rb")
    outf = open(outfile, 'w')
    
    for mm in snpPos2.index:
        QueryChrom, QueryPos, tmp, MatAllele, PatAllele, OverlappingGene = snpPos2.loc[mm,]
        QueryPos = int(QueryPos)-1

        for pileupcolumn in samfile.pileup(QueryChrom, QueryPos, QueryPos+1, truncate = True, stepper = 'all'):
            if pileupcolumn.pos == QueryPos:
                for pileupread in pileupcolumn.pileups:
                    Alignment = pileupread.alignment
                    ReadName = Alignment.query_name

                    if pileupread.query_position == None:
                        continue
                    BaseQualAsciiCode = ord(Alignment.qual[pileupread.query_position]) - 33   # Ascii code

                  #  MappingQual = Alignment.mapping_quality

                    Allele = Alignment.query_sequence[pileupread.query_position]

                   # print ReadName, QueryChrom, QueryPos, BaseQualAsciiCode, MappingQual,MatAllele, PatAllele,Allele
                    Allele_results = '\t'.join([ReadName,QueryChrom, str(QueryPos), str(BaseQualAsciiCode), MatAllele,PatAllele,Allele])+'\n'
                    outf.write(Allele_results)
    outf.close()
    samfile.close()
    
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
        Args_m = ["grep -w %s %s > %s/%s_split/%s.%s.Snps.tmp" % (chrom, Snp_File, outdir, sampleID, sampleID, chrom)]
        subprocess.check_call(Args_m, shell=True)
	Snp_tmp = '%s/%s_split/%s.%s.Snps.tmp' % (outdir, sampleID, sampleID, chrom)
        out_basename = '%s/%s_split/%s.%s.read_alleles.tmp' % (outdir, sampleID, sampleID, chrom)

        result = pool.apply_async(Get_alleles, (Snp_tmp, Bam_File, out_basename ))
        results_q.append(result)
        #Get_alleles(Snp_tmp, Bam_File, out_basename)

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
        
    for result in results_q:
        print(result.get())

    print time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))

    Args_m = ['cat %s/%s_split/%s.chr*.read_alleles.tmp > %s/%s.read_alleles.txt' % (outdir, sampleID, sampleID, outdir, sampleID)]
    subprocess.check_call(Args_m, shell=True)
    

    outdir_tmp = "%s/%s_split" % (outdir, sampleID)
    if os.path.exists(outdir_tmp):
        shutil.rmtree(outdir_tmp)

    print time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))
    print "Reads Allele separate finished!!!"
        
        
        

