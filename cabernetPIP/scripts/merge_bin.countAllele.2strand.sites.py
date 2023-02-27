import pandas as pd
import numpy as np
import os
import sys
import collections
from numpy import mean,median,ptp,var,std
from decimal import Decimal

bed_file = sys.argv[1]
bed_dict = collections.defaultdict(list)

#bed_f_tmp = pd.read_table(bed_file,header=None, index_col = None)
#bed_f_tmp['index'] = ["%s_%s" % (bed_f_tmp.loc[mm,0], bed_f_tmp.loc[mm,1]) for mm in bed_f_tmp.index]
#
#
#for ii,kk in enumerate(bed_f_tmp.index):
"""
chr1    3003720 3003721 .       -1      -1      .       -1      0
chr1    3003884 3003885 chr1    3003885 3003885 2       1.0     0
chr1    3003884 3003885 chr1    3003885 3003885 2       0.0     0
"""


ff = open(bed_file, 'r')
while 1:
    rline = ff.readline()
    if not rline:
	break

    line = rline.split()
    kk = line[0]+"_"+line[1]
    #print line
    if line[5]== "-1":
        tmp = []
        tmp.extend(line[:3])
        tmp.extend([0,"Nan",[]])
	#tmp.extend([0,"Nan"])
        bed_dict[kk] = tmp
        #print kk, bed_dict[kk]
    else:
        if bed_dict[kk] == []:
            bed_dict[kk].extend(line[:3])
            #bed_dict[kk].extend([line[6],0,[round(Decimal(line[7]),2)]])
	    bed_dict[kk].extend([line[6],0,[]])

	if line[7] == "0.5":
	    bed_dict[kk][4] = 1
	else:
	    bed_dict[kk][4] = 0

        bed_dict[kk][5].append(round(Decimal(line[7]),2))
	

#columns = ['feature', 'chrom', 'sta', 'end', 'ref_CpGs', 'detected_CpGs', 'detected_mCpGs', 'detected_umCpGs','Methylation_info', 'mean', 'median', 'ptp', 'var', 'std']
#columns = ['feature', 'chrom', 'sta', 'end', 'ref_CpGs', 'detected_CpGs', 'hemi_num', 'CpG_methyl_values', 'mean', 'median', 'var', 'std']
#print '\t'.join(columns)

for kk,vv in bed_dict.iteritems():
    #vv = map(str, vv)
    #print '\t'.join(vv)
    #print '\t'.(vv)
    values_v = vv[5]
    values_m = "Nan"
    allele_s = "Nan"

##    vv[6] = ','.join(map(str,vv[6]))
##    vv = map(str, vv)
#
#    new_v = [vv[0]+"_"+vv[1]]
#    new_v.extend(vv[:7])
    if values_v != []:
    	#mean_v = str(mean(values_v))
	mean_v = mean(values_v)
	#print len(values_v), vv[3]
	#print int(vv[3])== 2
		
	if len(values_v) == 2 and int(vv[3]) == 2:
		if mean_v == 1:
			values_m = 1
		elif mean_v == 0:
			values_m = 0
		elif values_v[0] != values_v[1]:
			values_m = 3
		else:
			values_m = 2
		allele_s = 1
	elif int(vv[3]) > 2:
		values_m = 4
		allele_s = 2
	else:
		if mean_v == 1:
			values_m = 1
		elif mean_v == 0:
			values_m = 0
		else:
			values_m = 2
		allele_s = 0

    vv.extend([values_m, allele_s])
    vv = map(str, vv)
    print '\t'.join(vv)
#    else:
#	new_v.extend(['', '', '', ''])
#
#    #print '\t'.join(new_v)
   



