# -*- coding: utf-8 -*-
import gzip
import pandas as pd
import numpy as np
import sys
import os
import re
from collections import defaultdict
from collections import Counter
from time import clock
from numpy import mean
import time
import signal
import multiprocess
import subprocess

#CpGs_f = '/share/home/baiyl/merlot/TAPS-merlot/workflow/Analysis/CNV/huang_200522_Xten_scEM-mouseTE_ICM_cultureM/Output_s1/Allele/ICM_1.CpG_read.bed'
outdir = sys.argv[1]
sampleId = sys.argv[2]
species = sys.argv[3]
threads = int(sys.argv[4])

print 'Input: %s' % CpGs_f
print 'Start process input'
outdir = os.path.split(os.path.realpath(allele_bed_f))[0]

def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def get_Methyls(CpGs_sub):
    Methyls = []
    for mm in set(CpGs_sub[3]):
        allele_sub = CpGs_sub.loc[(CpGs_sub[3]== mm),]
        methyl_value = mean(allele_sub[1])
       # methyl_value = mean(CpGs_sub.loc[(CpGs_sub[12]== mm),'Methyl'])
        #print methyl_value
        
        Himi_type = 'Full'
        if methyl_value != 0 and methyl_value !=1:
            Himi_type = 'Himi'
            methyl_3end = round(mean(allele_sub.loc[allele_sub[2]== 'off',1]))
            methyl_5end = round(mean(allele_sub.loc[allele_sub[2]== 'pre',1]))
            methyl_both = round(mean(allele_sub.loc[allele_sub[2]== 'both',1]))
            
            #methyl_3end = round(mean(CpGs_sub.loc[(CpGs_sub[12]== mm)&(CpGs_sub[8]== 'trim_off'),'Methyl']))
            #methyl_5end = round(mean(CpGs_sub.loc[(CpGs_sub[12]== mm)&(CpGs_sub[8]== 'trim_pre'),'Methyl']))
            #methyl_both = round(mean(CpGs_sub.loc[(CpGs_sub[12]== mm)&(CpGs_sub[8]== 'trim_both'),'Methyl']))
           # print methyl_3end, methyl_5end, methyl_both
            if methyl_3end!='nan' and not methyl_5end !='nan' and abs(methyl_3end-methyl_5end)>0.5:
                methyl_value = 0.5
            elif methyl_3end!='nan' and methyl_both!='nan' and abs(methyl_3end-methyl_both)>0.5:
                methyl_5end = methyl_both
                methyl_value = 0.5
            elif methyl_5end!='nan' and methyl_both!='nan' and abs(methyl_5end-methyl_both)>0.5:
                methyl_3end = methyl_both
                methyl_value = 0.5
            else:
                Himi_type = 'Puzz'
                if methyl_value < 0.3:
                    methyl_value = 0
                elif methyl_value > 0.7:
                    methyl_value = 1
                else:
                    methyl_value = 0.5

        Methyls.append([Himi_type, methyl_value,len(allele_sub)])
        
    return Methyls

def main(CpGs_file, CpGs_out):
    print 'Start get CpG allele methylation'
    t0=clock()
    size_mark = 100000
    aa = 0
    bb = 1
    
    CpG_s = ''
    Methyl_data = []
    CpGs_o2 = open(CpGs_out, 'w')
#    with gzip.open(CpGs_f,'rb') as ff:
    with open(CpGs_file,'rb') as ff:
        for line in ff:
            aa += 1
            if aa == size_mark * bb:
                print 'Processed %s cpgs in time %s min' % (bb*100000, (clock()-t0)/60)
                bb += 1
                
            cpg, methyl,trim,allele = line.split(' ')
            if CpG_s == '' or CpG_s == cpg:
                CpG_s = cpg
                Methyl_data.append([cpg, float(methyl),trim,allele])
                
            else:
                CpGs_sub = pd.DataFrame(Methyl_data)
                Methyls = get_Methyls(CpGs_sub)
		#allele_n = len(set(CpGs_sub[3]))
                allele_m = ','.join([':'.join(map(str,mm)) for mm in Methyls])
                mean_v = mean([mm[1] for mm in Methyls])
#                aaa = '\t'.join([CpG_s, str(len(Methyls)), allele_m, str(mean_v)])
		aaa = '\t'.join([CpG_s, str(allele_n), allele_m, str(mean_v)])

                CpGs_o2.write(aaa+"\n")
                
                CpG_s = cpg
                Methyl_data = [[cpg, float(methyl),trim,allele]]
           
    print 'Total spend  %s min' % str((clock()-t0)/60)
    CpGs_o2.close()
    
    
if __name__ == "__main__":
	if species == 'hg38':
		chrs = range(1,23)
	elif species == 'mm10':
		chrs = range(1,20)
	chrs.extend(['X','Y'])

	pool = multiprocess.Pool(threads,init_worker)
	results_q = []
	tmp_ReadOut = []
	t1=clock()


	for cc in chrs:
		print cc
		chrom = 'chr' + str(cc)

#		Args_m = ['grep -w %s %s > %s.%s.tmp' % (chrom, CpGs_f, CpGs_f, chrom)]
#		subprocess.check_call(Args_m, shell=True)
		
		CpGs_f_chr = '%s/%s.%s.cpgs_read.bed' % (outdir, sampleId, chrom)
		CpGs_o_chr = '%s/%s.%s.allele_methyl.txt' % (outdir, sampleId, chrom)

#		CpGs_f_chr = '%s.%s.tmp' % (CpGs_f, chrom)
#		CpGs_o_chr = '%s.%s.tmp' % (CpGs_o, chrom)
		if not os.path.exists(CpGs_f_chr):
			continue
		elif not os.path.getsize(CpGs_f_chr):
                        continue

		main(CpGs_f_chr, CpGs_o_chr)
#		result = pool.apply_async(main, (CpGs_f_chr, CpGs_o_chr ))
		tmp_ReadOut.extend([CpGs_f_chr, CpGs_o_chr])
   
    
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
    
	Args_m = ['cat %s/%s.chr*.allele_methyl.txt > %s/%s.allele_methyl.txt' % (outdir, sampleId, outdir, sampleId)]
	subprocess.check_call(Args_m, shell=True)
#	for my_file in tmp_ReadOut:
#            if os.path.exists(my_file):
#                os.remove(my_file)
	t2=clock()

	print "Final get CpG allele %s in %s" % (sampleId, str(t2 - t1))

 
