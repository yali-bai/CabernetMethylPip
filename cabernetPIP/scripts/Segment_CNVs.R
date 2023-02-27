
#'''
#library(DNAcopy)
#data(coriell)
#CNA.object <- CNA(cbind(coriell$Coriell.13330),
#                  coriell$Chromosome,coriell$Position,
#                  data.type="logratio",sampleid="c13330")
#smoothed.CNA.object <- smooth.CNA(CNA.object)
#segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)
#'''

library(DNAcopy)
library(data.table)
#library(tidyverse)
library(dplyr)

#cnv_f = '/Users/Bai/bai_project/XIE_Group/TAPS-merlot/1_Analysis/ESCC_smallBulk/CNVs/ESCC.1Mb.nor.txt'
#cnv_f = "/Users/aliyabai/Project/TAPS-merlot/10\ Mouseq_Embryo/CNVs/Merged_CNVs/D20.1Mb.nor.txt"
#outdir = "/Users/aliyabai/Project/TAPS-merlot/10\ Mouseq_Embryo/CNVs/Merged_CNVs/"


args = commandArgs(T)
cnv_f = args[1]
outdir = args[2]
CNVs_data <- fread(cnv_f)



Samples = colnames(CNVs_data)[5:length(CNVs_data)]
for(ss in Samples){
  #for(ss in colnames(CNVs_data)[5:length(CNVs_data)]){
  new_id = strsplit(ss, ".correctCH.rmdup")[[1]]
  #new_id = ss
  ss_data <- CNA(cbind(CNVs_data %>% select(all_of(ss))),
                 CNVs_data$chromosome,
                 CNVs_data$start,
                 data.type = 'logratio', 
                 sampleid = new_id
  )
  ss_data.seg <- segment(ss_data, verbose = 1)
  ss_data.seg.sum <- segments.summary(ss_data.seg)
  
  write.csv(ss_data.seg.sum, 
            paste0(outdir,"/seg.",new_id,'.tsv'),
            row.names = TRUE)
  
}




#plot(ESCC.CNA.seg)
#zoomIntoRegion(ESCC.CNA.seg, chrom=1, sampleid= "ESCC_119N")



