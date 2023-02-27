import sys
import os
from collections import defaultdict
import pandas as pd
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

def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)


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

parser.add_argument("--minAlleleNum",
                    type=int,
                    default = 1,
                    help="The min Allele Number in a read. default: 1")

parser.add_argument("--minDV",
                    type=float,
                    default = 0.2,
                    help="The min difference value(D-Values) of paternal and maternal. default:0.2")

parser.add_argument("--minBQ",
                    type=int,
                    default = 30,
                    help="The min base quality of allele in the reads.")

parser.add_argument("--threads",
                        type=int,
                        default = 1,
                        help="thread, default 1")

args = parser.parse_args(args)
allele_context_f = args.input
outdir = args.outdir
sampleID = args.sampleID
species = args.species
minAlleleNum = args.minAlleleNum
minDV = args.minDV
minBQ = args.minBQ
threads = args.threads

if species == 'hg38':
    chrs = range(1,23)
elif species == 'mm10':
    chrs = range(1,20)
    chrs.extend(['X','Y'])


def Read_Sep(Strand, reads_Pos_alleles):
    Filter_Pos_alleles = reads_Pos_alleles.loc[reads_Pos_alleles[1] > minBQ,]
    if Filter_Pos_alleles.shape < 1:
        Filter_Pos_alleles = reads_Pos_alleles
    
    if Strand == 'CT':
        Change_alleles = Filter_Pos_alleles.loc[(Filter_Pos_alleles[2] == "C")|(Filter_Pos_alleles[3] == "C"),]
        nonChange_alleles = Filter_Pos_alleles.loc[(Filter_Pos_alleles[2] != "C")&(Filter_Pos_alleles[3] != "C"),]
        ref_allele = 'C'
        ref_allele2 = 'T'
        ref_pair = ['C','T']
    else:
        Change_alleles = Filter_Pos_alleles.loc[(Filter_Pos_alleles[2] == "G")|(Filter_Pos_alleles[3] == "G"),]
        nonChange_alleles = Filter_Pos_alleles.loc[(Filter_Pos_alleles[2] != "G")&(Filter_Pos_alleles[3] != "G"),]
        ref_allele = 'G'
        ref_allele2 = 'A'
        ref_pair = ['A','G']

    Change_alleles = Change_alleles.loc[[mm for mm in Change_alleles.index if sorted(Change_alleles.loc[mm,2:3].values)!= ref_pair],]
    Total_alleles = Change_alleles.shape[0] + nonChange_alleles.shape[0]
    Maternal_alleles = 0
    Paternal_alleles = 0
    maAs = 0
    paAs = 0
    
    if Total_alleles == 0:
        return Total_alleles, Maternal_alleles, Paternal_alleles

    if Change_alleles.shape[0] >0:
        maA = Change_alleles.loc[Change_alleles[2] == ref_allele,]
        maAs = maA.loc[(maA[4]==ref_allele)|(maA[4]==ref_allele2),].shape[0]
	#paAs += maA.loc[(maA[4]!=ref_allele)&(maA[4]!=ref_allele2),].shape[0]
	paAs += maA.loc[maA[4]==maA[3],].shape[0]
        paA = Change_alleles.loc[Change_alleles[3] == ref_allele,]
        paAs += paA.loc[(paA[4]==ref_allele)|(paA[4]==ref_allele2),].shape[0]
	#maAs += paA.loc[(paA[4]!=ref_allele)&(paA[4]!=ref_allele2),].shape[0]
	maAs += paA.loc[paA[4] == paA[2],].shape[0]

    if nonChange_alleles.shape[0] >0:
        Maternal_alleles = nonChange_alleles.loc[nonChange_alleles[2] == nonChange_alleles[4],].shape[0]
        Paternal_alleles = nonChange_alleles.loc[nonChange_alleles[3] == nonChange_alleles[4],].shape[0]
    
    Maternal_alleles = maAs + Maternal_alleles
    Paternal_alleles = paAs + Paternal_alleles
    
    return Total_alleles, Maternal_alleles, Paternal_alleles


def Sep_parental_reads(read_allels_f, outfile):
	#read_allels_f="/share/home/baiyl/merlot/TAPS-merlot/workflow/Analysis/allele_dropout/SNP_ado/SRC/test_data/D200914_hmC_L1C_5A_S1.2.read_allels.tmp"
#readID Alleles
#A100:chr1:r155:trim_off:CT:42,42 3012372,41,T,C,T:3012372,41,T,C,T
#A100:chr1:r156:trim_pre:GA:42,42 3012372,41,T,C,T:3012372,41,T,C,T
#A104:chr1:r165:trim_pre:GA:42,42 3012869,41,T,G,A:3012869,41,T,G,A
#A105:chr1:r166:trim_pre:GA:42,42 3012869,37,T,G,T:3012869,41,T,G,T
#A106:chr1:r167:trim_pre:GA:42,42 3012869,41,T,G,T
#A111:chr1:r174:trim_off:CT:40,40 3014189,27,G,T,G
#A111:chr1:r175:trim_off:CT:42,42 3014189,37,G,T,G

    outf = open(outfile, 'w')
    read_allels = pd.read_table(read_allels_f, header = 0, index_col = None,sep =' ')
    #read_allels = pd.read_table(read_allels_f, header = None, index_col = None,sep =' ')

#    aa = 0
    for mm in read_allels.index:
#        aa +=1
#        if aa > 20:
#            break
        readID = read_allels.loc[mm,'readID']
        reads_Pos_alleles = pd.DataFrame([re.split(',',allele) for allele in re.split(':',read_allels.loc[mm,'Alleles'])])

        AlleleN, chrom, readN, strand, MapQ = re.split(":",readID)

        Total_alleles, Maternal_alleles, Paternal_alleles = Read_Sep(strand, reads_Pos_alleles)
        Maternal_aR = Paternal_aR = 0
        if Total_alleles >= minAlleleNum:
            Maternal_aR = Maternal_alleles / Total_alleles
            Paternal_aR = Paternal_alleles / Total_alleles

            if Maternal_aR - Paternal_aR > minDV:
                Parental = "M"
            elif Paternal_aR - Maternal_aR > minDV:
                Parental = "P"
            else:
                Parental = "N"
        else:
            Parental = "N"
     #   print readID, Parental, Total_alleles, Maternal_aR, Paternal_aR
        Sep_out = '\t'.join([readID, Parental, str(Total_alleles), str(Maternal_alleles), str(Paternal_alleles)])+'\n'
        outf.write(Sep_out)
        read_allels.drop(index=mm)

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

    print time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))
    for cc in chrs:
        chrom = 'chr' + str(cc)
        
        #"test_data/D200914_hmC_L1C_5A_S1.2_split/D200914_hmC_L1C_5A_S1.2.chr10.read_allels.tmp |awk 'BEGIN{ID="readID";alleles="Alleles"}{if($1==ID){alleles=alleles":"$3","$4","$5","$6","$7;}else{print ID,alleles;ID=$1;alleles=$3","$4","$5","$6","$7}}'"
        Args_m = ["grep -w %s %s|sort |awk 'BEGIN{ID=\"readID\";alleles=\"Alleles\"}{if($1==ID){alleles=alleles\":\"$3\",\"$4\",\"$5\",\"$6\",\"$7;}else{print ID,alleles;ID=$1;alleles=$3\",\"$4\",\"$5\",\"$6\",\"$7}}' > %s/%s_split/%s.%s.allele_context.tmp" % (chrom,allele_context_f, outdir, sampleID, sampleID, chrom) ]
        subprocess.check_call(Args_m, shell=True)
        
        allele_context_tmp = '%s/%s_split/%s.%s.allele_context.tmp' % (outdir, sampleID, sampleID, chrom)
        out_basename = '%s/%s_split/%s.%s.read_phased.tmp' % (outdir, sampleID, sampleID, chrom)

        result = pool.apply_async(Sep_parental_reads, (allele_context_tmp, out_basename ))
        results_q.append(result)
       # Sep_parental_reads(allele_context_tmp, out_basename )
        
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
    Args_m = ['cat %s/%s_split/%s.chr*.read_phased.tmp |gzip > %s/%s.Read_Phased.txt.gz' % (outdir, sampleID, sampleID, outdir, sampleID)]
    subprocess.check_call(Args_m, shell=True)

    outdir_tmp = "%s/%s_split" % (outdir, sampleID)
    #if os.path.exists(outdir_tmp):
    #    shutil.rmtree(outdir_tmp)

    print time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))
    print "Separate Reads to Parental by alleles finished!!!"




