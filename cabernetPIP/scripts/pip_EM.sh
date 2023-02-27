#!/usr/bin/bash
## langurage
python='/share/newdata4/baiyl/basic/tools/anaconda2/bin/python'
Rscript=/share/newdata4/baiyl/basic/tools/anaconda2/lib/R/bin/Rscript

## soft Path
bismark="/share/newdata4/baiyl/basic/tools/Bismark_v0.19.0/bismark"
bowtie2="/share/newdata4/baiyl/basic/tools/anaconda2/bin"
picard=~/./basic/tools/picard-tools-2.4.1/picard.jar
bismark_methylation_extractor="/share/newdata4/baiyl/basic/tools/Bismark_v0.19.0/bismark_methylation_extractor"
samtools='/share/newdata4/baiyl/basic/tools/anaconda2/bin/samtools'
bedtools='/share/newdata4/baiyl/basic/tools/bedtools2/bin/bedtools'
sambamba="/share/newdata4/baiyl/basic/tools/sambamba_v0.6.6"

## home made code
SRC_DIR=`pwd`/../
SRC_DIR=/share/home/baiyl/merlot/TAPS-merlot/workflow/RUN_src/src_2.0

Trim_script=${SRC_DIR}/scripts/trim_fq.sh
Remove_mC_filledGap=$SRC_DIR/scripts/Remove_mC_filledGap.py
Allele_count_py=${SRC_DIR}/scripts/allele_count.plus.v01.py
Call_CpG_Hemi=${SRC_DIR}/scripts/Call_CpG_hemi.py
Phase_SNP_read=${SRC_DIR}/scripts/Phase_SNP_read.v01.py
Sep_parental_reads=${SRC_DIR}/scripts/Sep_parental_reads.py
Sep_CpGs=${SRC_DIR}/scripts/Sep_CpGs.py
countAllele2strand=${SRC_DIR}/scripts/merge_bin.countAllele.2strand.sites.py

## reference
#----------------------- Bismark Ref --------------------
hg38_bis=/share/home/baiyl/basic/database/human/hg38/Bismark_hg38
mm10_bis=/share/home/baiyl/basic/database/mouse/mm10/Index/Bismark_mm10
full_puc19_bis=/share/home/baiyl/basic/database/puC19/puc19_bismark
Ecoli_bis=/share/home/baiyl/basic/database/cut_ecoli/ecoli_bismark
lambda_bis=/share/home/baiyl/basic/database/NEB_lambda/lambda_Bismark
msp1_bis=/share/home/baiyl/basic/database/msp1/msp1_bwaindex
hmC5_bis=/share/home/baiyl/database/5hmC_zymo/hmC5_bismark
clai_bis=/share/home/baiyl/database/clai/bismark_clai

#----------------------- Bwa Ref ----------------------
hg38_ref_fa=/share/newdata4/baiyl/basic/database/human/hg38/Index/bwa/hg38.fa
mm10_ref_fa=/share/home/baiyl/basic/database/Illumina_iGENOMES/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa
lambda_ref_fa=/share/home/baiyl/database/NEB_lambda/BWAIndex/genome.fa
hmC5_ref_fa=/share/home/baiyl/database/5hmC_zymo/BWAIndex/5hmC_zymo.fa
full_puc19=/share/home/baiyl/basic/database/puC19/BWAIndex/genome.fa
clai_fa=/share/home/baiyl/database/clai/clai.fa

#----------------------- Blacklist --------------------
hg38_blacklist=/share/home/baiyl/basic/database/human/blacklist/hg38.blacklist.bed
mm10_blacklist=/share/home/baiyl/database/mouse/mm10/Blacklist/mm10.blacklist.bed

#----------------------- ref_snp --------------------
ref_snp=/share/home/baiyl/merlot/TAPS-merlot/workflow/Analysis/allele_dropout/SNP_ado/SRC/data/C57BL_6NJ.DBA_1J.snp.bed


## variable
declare -A Ref_bis
Ref_bis["hg38"]=$hg38_bis
Ref_bis["mm10"]=$mm10_bis
Ref_bis["hg38_mm10"]=$hg38_mm10_bis
Ref_bis["hg38_ecoli_puc19"]=$hg38_ecoli_puc19_bis
Ref_bis["mm10_ecoli_puc19"]=$mm10_ecoli_puc19_bis
Ref_bis["fullpuc19"]=$full_puc19_bis
Ref_bis["puc19"]=$song_puc19_bis
Ref_bis["ecoli"]=$Ecoli_bis
Ref_bis["lambda"]=$lambda_bis
Ref_bis["msp1"]=$msp1_bis
Ref_bis["hmC5"]=$hmC5_bis
Ref_bis["clai"]=$clai_bis

declare -A Ref_fa
Ref_fa["hg38"]=$hg38_ref_fa
Ref_fa["mm10"]=$mm10_ref_fa
Ref_fa["puc19"]=$cutpuc19_ref_fa
Ref_fa["lambda"]=$lambda_ref_fa
Ref_fa["ecoli"]=$cutecoli_ref_fa
Ref_fa["hg38_mm10"]=$hg38_mm10_ref
Ref_fa["hg38_mm10_ecoli_puc19"]=$hg38_mm10_ecoli_puc19
Ref_fa["hg38_ecoli_puc19"]=$hg38_ecoli_puc19
Ref_fa["mm10_ecoli_puc19"]=$mm10_ecoli_puc19
Ref_fa["fullpuc19"]=$full_puc19
Ref_fa["hmC5"]=$hmC5_ref_fa
Ref_fa["clai"]=$clai_fa

declare -A Blacklist
Blacklist["hg38"]=$hg38_blacklist
Blacklist["mm10"]=$mm10_blacklist


m_threshold=0.5



## Args
temp=`getopt -o i:s:p:m:o:t: --long indir:,sep_type:,prefix:,meta_type:,tag_file:,outdir:,trim_type:,read1:,read2:,read_len:,clip_len:,thread:,species:,seq_type:,partition:  -- "$@"`
if [ $? != 0 ] ; then echo "terminating..." >&2 ; exit 1 ; fi
eval set -- "$temp"

while true ; do
        case "$1" in
		-i|--indir) echo "Input dir is $2"; indir=$2; shift 2;;
                -s|--sep_type) sep_type=$2; echo "Separate the reads with barcodes. Default: meta_mark" ; shift 2;;
                -p|--prefix) prefix=$2; echo "Reads prefix is $prefix" ; shift 2;;
                -m|--meta_type) meta_type=$2; echo "meta type is $meta_type"; shift 2 ;;
                --tag_file) tag_file=$2; echo "Tag file is $tag_file"; shift 2 ;;

                -o|--outdir) echo "Output dir is $2" ;outdir=$2; shift 2;;
                --trim_type) trim_type=$2;echo "option b, argument \`$2'" ; shift 2;;
                --read1)
                        case "$2" in
                                "") echo "option read1, no argument"; shift 2 ;;
                                *) Read1=$2; echo "Read1 is $Read1" ; shift 2 ;;
                        esac ;;
                --read2)
                        case "$2" in
                                "") echo "option read2, no argument"; shift 2 ;;
                                *) Read2=$2; echo "Read2 is $2" ; shift 2 ;;
                        esac ;;
                --read_len) read_len=$2; echo "minReadlen is $2" ; shift 2 ;;
                --clip_len) clip_l=$2; echo "minCliplen is $2" ; shift 2 ;;
		--species) species=$2; echo "Species is $2"; shift 2 ;;
		--partition) partition=$2; echo "partition is $2"; shift 2 ;;
		--seq_type) seq_type=$2; echo "seq_type is $2"; shift 2 ;;
		--CallHemi) seq_type=$2; echo "CallHemi is $2"; shift 2 ;;
		--Parental_Phase) seq_type=$2; echo "Parental_Phase is $2"; shift 2 ;;
		-t|--thread) thread=$2; echo "multi-thread is $thread"; shift 2 ;;

                --) shift ; break ;;
                *) echo "internal error!" ; exit 1 ;;
        esac
done


### -------------------------------------------------------------------
###### judge the args
###
if [ ! -s $Read1 ];then
        echo $Read1 "Read1 not exist"
        exit
fi
if [ ! -s $Read2 ];then
        echo $Read2 "Read2 not exist"
        exit
fi

if [[ ! -n $read_len ]];then read_len=20;fi
if [[ ! -n $trim_type ]];then trim_type=basic_trim;fi
if [[ ! -n $outdir ]];then  outdir=".";fi
if [[ ! -n $species ]];then echo "Please choose a species";exit;fi
if [[ ! -n $seq_type ]];then seq_type=nextera;fi


cores=`expr $SLURM_CPUS_PER_TASK / 5`
if [[ $cores == 0 ]];then
	 cores=1
fi

ref_fa=${Ref_fa[$species]}
ref_bis=${Ref_bis[$species]}
Blacklist=${Blacklist[$species]}

sampleid=`echo $Read1|sed 's/.*\///g'|sed 's/[._]R1.f.*//g'|sed 's/.R1_val_1.fq.gz//g'`
echo $sampleid
id=$sampleid

if [[ ! -f $SRC_DIR/data/genome.length.txt ]];then
        genome_len=`less $ref_fa |grep -v ">"|awk 'BEGIN{sum=0}{sum+=length($1)}END{print sum}'|xargs`
        echo -e "${species}\t$genome_len" >> $SRC_DIR/data/genome.length.txt
else
        genome_len=`grep -w ${species} $SRC_DIR/data/genome.length.txt|awk '{print $2}'|uniq|xargs`

        if [[ $genome_len == "" ]];then
                genome_len=`less $ref_fa |grep -v ">"|awk 'BEGIN{sum=0}{sum+=length($1)}END{print sum}'|xargs`
                echo -e "${species}\t$genome_len" >> $SRC_DIR/data/genome.length.txt

        fi
fi
CpG_num=`grep -w ${species} $SRC_DIR/data/genome.CpG_num.txt|awk '{print $2}'|uniq|xargs`

###
### --------------- functions -----------------
###
function trim_array {
mkdir -p $outdir/trim
t1=`ls $outdir/trim/${sampleid}[._]R1_val_1.fq.gz`
t2=`ls $outdir/trim/${sampleid}[._]R2_val_2.fq.gz`
echo "Reads need to be trimed: $Read1, $Read2"

if [[ -f $t1 ]] && [[ -f $t2 ]];then
        echo "Trim already finished"
        return
else
       rm -fr $outdir/trim/${sampleid}*
fi



if [[ ! -s $clip_len ]];then
bash $Trim_script \
        -o $outdir/trim \
        -t $trim_type \
        --read1 $Read1 \
        --read2 $Read2 \
        --read_len $read_len
else
bash $Trim_script \
        -o $outdir/trim \
        -t $trim_type \
        --read1 $Read1 \
        --read2 $Read2 \
        --read_len $read_len \
	--clip_len $clip_len
fi

}

function bismark_call {
trim_Read1=$outdir/trim/${sampleid}.R1_val_1.fq.gz
trim_Read2=$outdir/trim/${sampleid}.R2_val_2.fq.gz
mkdir -p $outdir/bismark_${species}

$bismark --multicore $cores \
        --fastq --non_directional --unmapped \
        --nucleotide_coverage \
        --path_to_bowtie $bowtie2 --bowtie2 \
        --genome_folder $ref_bis \
        --output_dir $outdir/bismark_${species} \
        --temp_dir $outdir/bismark_${species}/temp_bismark_${sampleid} \
	-1 $trim_Read1 -2 $trim_Read2


outbam=`ls $outdir/bismark_${species}/${sampleid}*_bismark_bt2_pe.bam|xargs`


if [[ $species == "hg38" ]] || [[ $species == "mm10" ]]
then
$samtools view -bh -q 1 -F 4 $outbam |$samtools sort -@ $SLURM_CPUS_PER_TASK |$bedtools intersect -a - -b $Blacklist -v > $outdir/bismark_${species}/${sampleid}.${species}.bam
else
$samtools view -bh -q 1 -F 4 $outbam |$samtools sort -@ $SLURM_CPUS_PER_TASK > $outdir/bismark_${species}/${sampleid}.${species}.bam
fi


$samtools index $outdir/bismark_${species}/${sampleid}.${species}.bam
$samtools stats $outdir/bismark_${species}/${sampleid}.${species}.bam > $outdir/bismark_${species}/${sampleid}.${species}.bam.stat
$samtools depth $outdir/bismark_${species}/${sampleid}.${species}.bam > $outdir/bismark_${species}/${sampleid}.${species}.depth

}



function Remove_mC_filledGap {

echo "Remove the mC filled Gap region caused by Tn5"
Read_type=$1  ## PE or SE
inbam=$outdir/bismark_${species}/${sampleid}.${species}.bam

date
if [[ ! -s $outdir/bismark_${species}/${sampleid}.correctCH.bam ]];then
$python $Remove_mC_filledGap \
        -i $inbam \
        -o $outdir/bismark_${species}/${sampleid}.correctCH.bam \
        -r $ref_fa \
        -g 9 \
        --thread $SLURM_CPUS_PER_TASK \
        --Read_type $Read_type
date
fi

$samtools sort -@ $SLURM_CPUS_PER_TASK $outdir/bismark_${species}/${sampleid}.correctCH.bam > $outdir/bismark_${species}/${sampleid}.correctCH.sorted.bam
$samtools index $outdir/bismark_${species}/${sampleid}.correctCH.sorted.bam
$samtools stats $outdir/bismark_${species}/${sampleid}.correctCH.sorted.bam > $outdir/bismark_${species}/${sampleid}.correctCH.sorted.stat
$samtools depth $outdir/bismark_${species}/${sampleid}.correctCH.sorted.bam > $outdir/bismark_${species}/${sampleid}.correctCH.sorted.depth

date
## remove duplication reads
java -jar $picard MarkDuplicates \
       I=${outdir}/bismark_${species}/${id}.correctCH.sorted.bam \
       O=${outdir}/bismark_${species}/${id}.correctCH.rmdup.bam \
       M=${outdir}/bismark_${species}/${id}.correctCH.rmdup.txt \
       READ_NAME_REGEX=null \
       REMOVE_DUPLICATES=true
$samtools index ${outdir}/bismark_${species}/${id}.correctCH.rmdup.bam
date
}

function call_methyl {
echo "Starting call the CpG methylation..."
date

$samtools sort -n -@ $SLURM_CPUS_PER_TASK ${outdir}/bismark_${species}/${id}.correctCH.rmdup.bam > ${outdir}/bismark_${species}/${id}.correctCH.rmdup.sortn.bam
$samtools index ${outdir}/bismark_${species}/${id}.correctCH.rmdup.sortn.bam
date

$bismark_methylation_extractor -p \
   --multicore $cores \
   --comprehensive \
   --no_overlap \
   --bedGraph \
   --counts \
   --buffer_size 20G \
   --report \
   --cytosine_report \
   --genome_folder $ref_bis \
   -o ${outdir}/bismark_${species} \
   ${outdir}/bismark_${species}/${id}.correctCH.rmdup.sortn.bam

date

mkdir -p $outdir/sites


CpG_f=`ls ${outdir}/bismark_${species}/${id}.*CpG_report.txt.CpG_report.txt|grep correctCH|grep rmdup |xargs`

awk 'BEGIN{chr="Chr";pos="Pos";methyl="Methyl";cov="UMethyl"}{OFS="\t";if($3=="-"){$2=$2-1};if(pos==$2){methyl+=$4;cov+=$5}else{sum=methyl+cov;if(sum==0){rr=""}else{rr=methyl/sum};print chr,pos,methyl,cov,rr;chr=$1;pos=$2;methyl=$4;cov=$5}}' $CpG_f  > $outdir/sites/${sampleid}.${species}.CpG_methy.txt
grep -v Pos $outdir/sites/${sampleid}.${species}.CpG_methy.txt|grep -v "_" |sort -V > $outdir/sites/${sampleid}.${species}.CpG_methy.tmp
cat <(echo -e "${sampleid}.cov\t${sampleid}.ratio") <(cut -f 5 $outdir/sites/${sampleid}.${species}.CpG_methy.tmp) > $outdir/sites/${sampleid}.${species}.CpG_methy.csv

cat <(echo -e "chrom\tPos\t") <(cut -f 1,2 $outdir/sites/${sampleid}.${species}.CpG_methy.tmp) > $outdir/sites/index.csv

}

function count_allele {
## identify the reads from the same allele

mkdir -p ${outdir}/Count_allele

date
aa=`$samtools view ${outdir}/Count_allele/${id}.correctCH.rmdup.allele.bam|head -n 5|wc -l`
if [[ ! -s ${outdir}/Count_allele/${id}.correctCH.rmdup.allele.bam ]] || [[ $aa == 0 ]];then
$python $Allele_count_py \
        --thread $SLURM_CPUS_PER_TASK \
        --seq_type $seq_type \
	-i ${outdir}/bismark_${species}/${id}.correctCH.rmdup.bam \
	-o ${outdir}/Count_allele/${id}.correctCH.rmdup.allele.bam \
	-a ${outdir}/Count_allele/${id}.correctCH.rmdup.allele.txt \
	--MapQ_flag \
	-m 2
date
fi

bb=`$samtools view ${outdir}/Count_allele/${id}.correctCH.rmdup.allele.sortn.bam|head -n 5|wc -l`
if [[ ! -s ${outdir}/Count_allele/${id}.correctCH.rmdup.allele.sortn.bam ]] || [[ $bb == 0 ]];then
$samtools sort -n -@ $SLURM_CPUS_PER_TASK ${outdir}/Count_allele/${id}.correctCH.rmdup.allele.bam > ${outdir}/Count_allele/${id}.correctCH.rmdup.allele.sortn.bam
$samtools index ${outdir}/Count_allele/${id}.correctCH.rmdup.allele.sortn.bam
$samtools sort -@ $SLURM_CPUS_PER_TASK ${outdir}/Count_allele/${id}.correctCH.rmdup.allele.bam > ${outdir}/Count_allele/${id}.correctCH.rmdup.allele.sorted.bam
$samtools index ${outdir}/Count_allele/${id}.correctCH.rmdup.allele.sorted.bam

$samtools view -h ${outdir}/Count_allele/${id}.correctCH.rmdup.allele.bam|awk '$1~"^A"||$1~"^@"'|$samtools sort > ${outdir}/Count_allele/${id}.correctCH.rmdup.alleleA.bam
$samtools depth ${outdir}/Count_allele/${id}.correctCH.rmdup.alleleA.bam > ${outdir}/Count_allele/${id}.correctCH.rmdup.alleleA.depth

date
fi

if [[ ! -s ${outdir}/Count_allele/CpG_context_${id}.correctCH.rmdup.allele.sortn.txt ]];then
$bismark_methylation_extractor -p \
   --multicore $cores \
   --comprehensive \
   --no_overlap \
   --bedGraph \
   --counts \
   --buffer_size 20G \
   --report \
   --cytosine_report \
   --genome_folder $ref_bis \
   -o ${outdir}/Count_allele \
   ${outdir}/Count_allele/${id}.correctCH.rmdup.allele.sortn.bam
fi
mkdir -p $outdir/sites

CpG_f=`ls ${outdir}/Count_allele/${id}.*CpG_report.txt.CpG_report.txt|grep correctCH|grep rmdup |grep allele|xargs`
awk 'BEGIN{chr="Chr";pos="Pos";methyl="Methyl";cov="UMethyl"}{OFS="\t";if($3=="-"){$2=$2-1};if(pos==$2){methyl+=$4;cov+=$5}else{sum=methyl+cov;if(sum==0){rr=""}else{rr=methyl/sum};print chr,pos,methyl,cov,rr;chr=$1;pos=$2;methyl=$4;cov=$5}}' $CpG_f  > $outdir/sites/${sampleid}.${species}.CpG_methy.txt
grep -v Pos $outdir/sites/${sampleid}.${species}.CpG_methy.txt|grep -v "_" |sort -V > $outdir/sites/${sampleid}.${species}.CpG_methy.tmp
cat <(echo -e "${sampleid}.cov\t${sampleid}.ratio") <(cut -f 5 $outdir/sites/${sampleid}.${species}.CpG_methy.tmp) > $outdir/sites/${sampleid}.${species}.CpG_methy.csv
if [[ ! -f $outdir/sites/index.csv ]];then
cat <(echo -e "chrom\tPos\t") <(cut -f 1,2 $outdir/sites/${sampleid}.${species}.CpG_methy.tmp) > $outdir/sites/index.csv
fi

date

awk '$1!~"^M"{split($1,aa,":");if(aa[4]=="GA"){$4=$4-1};print $3,$4,$5,aa[1],aa[4],aa[5]}' ${outdir}/Count_allele/CpG_context_${id}.correctCH.rmdup.allele.sortn.txt |sed 's/Z/1/g'|sed 's/z/0/g'|sort -V > ${outdir}/Count_allele/CpG_context_${id}.allele.clean.txt
set -x
if [[ -s ${outdir}/Count_allele/${id}.allele.CpG_allele_m.posMerge.txt.gz ]];then
aa=`less ${outdir}/Count_allele/${id}.allele.CpG_allele_m.posMerge.txt.gz|cut -f 1|uniq|wc -l|xargs`
label=no
if [[ $aa != 21 ]] && [[ $species == "mm10" ]];then label=yes;fi
if [[ $aa != 24 ]] && [[ $species == "hg38" ]];then label=yes;fi
else
label=yes
fi
date
echo "Get the CpG Hemi values"
if [[ $label == "yes" ]];then
$python $Call_CpG_Hemi \
        -i ${outdir}/Count_allele/CpG_context_${id}.allele.clean.txt \
        -o ${outdir}/Count_allele \
        --sampleID ${id}.allele \
        --species $species \
        --MapQ_flag \
        --threads $SLURM_CPUS_PER_TASK
fi
date

if [[ ! -f ${outdir}/stat/chrom.allele.stat.txt ]];then
echo -e "SampleID Chrom chr_mean_allele_methyl allele2_2strand_CpGs CpGs_2strand_hemi CpGs_2strand_hemi_ratio CpGs_2strand_fullM CpGs_2strand_fullM_ratio CpGs_2strand_fullUM CpGs_2strand_fullUM_ratio" > ${outdir}/stat/chrom.allele.stat.txt
fi

chrs=(`seq 1 19` "X" "Y")
for cc in ${chrs[@]}
do
       cc_mean_allele_methyl=`awk -v cc=chr$cc '$1==cc{methyl+=$3;sum+=1}END{print 100*methyl/sum"%"}' ${outdir}/Count_allele/CpG_context_${id}.allele.clean.txt`
       allele2_2strand=`less ${outdir}/Count_allele/${id}.allele.CpG_allele_m.sameAllele.txt.gz |awk -v cc=chr${cc} '$3<3&&$1==cc {if($5==0.5){Hemi+=1}else if($5==1){full_m+=1}else{full_u+=1}}END{all=Hemi+full_m+full_u;print all,Hemi,Hemi/all,full_m,full_m/all,full_u,full_u/all}'`
	echo -e "${id} $cc $cc_mean_allele_methyl $allele2_2strand" >> ${outdir}/stat/chrom.allele.stat.txt
done
date

set -x
Alleledepth=$(wc -l ${outdir}/Count_allele/${id}.correctCH.rmdup.alleleA.depth|awk '{print $1}')
AlleleReadCov=`echo "scale=4;100*$Alleledepth/$genome_len"|bc|xargs`

CpGs_Allele_type=`less ${outdir}/Count_allele/${id}.allele.CpG_allele_m.posMerge.txt.gz|awk -v CpG_num=$CpG_num 'BEGIN{a1=0;a2=0;a3=0}{if($3==1){a1+=1}else if($3==2){a2+=1}else{a3+=1}}END{sum=a1+a2+a3;print sum,sum/CpG_num,a1/CpG_num,a2/CpG_num,a3/CpG_num}'`
Mean_Allele_Methyl=`awk '{sum+=$3}END{print 100*sum/NR"%"}' Output_s1/Count_allele/CpG_context_${id}.allele.clean.txt`

CpGs_2allele_2strand=`less ${outdir}/Count_allele/${id}.allele.CpG_allele_m.sameAllele.txt.gz |awk '$3<3 {print $1,$2}'|sort|uniq|wc -l`
CpGs_2allele_2strand_Cov=`echo "scale=4;100*$CpGs_2allele_2strand/$CpG_num"|bc|xargs`
allele2_strand2_type=`less ${outdir}/Count_allele/${id}.allele.CpG_allele_m.sameAllele.txt.gz |awk '$3<3{if($5==0.5){Hemi+=1}else if($5==1){full_m+=1}else{full_u+=1}}END{all=Hemi+full_m+full_u;print all,Hemi,Hemi/all,full_m/all,full_u/all}'`

if [[ ! -f ${outdir}/stat/Allele_methyl.stat.txt ]];then
echo "SampleID species AlleleReadCov Mean_Allele_Methyl CpGs_with_2strand_Cov Alleles_2strand Alleles_2strand_Hemi HemiAllelesRatio FullMAllelesRatio FullUMAllelesRatio" > ${outdir}/stat/Allele_methyl.stat.txt
fi

echo -e "${id} ${species} " \
  "${AlleleReadCov}% " \
  "${Mean_Allele_Methyl} " \
  "${CpGs_2allele_2strand_Cov}% " \
  "${allele2_strand2_type}" >> ${outdir}/stat/Allele_methyl.stat.txt

set +x
}

function Parental_Phase {
## Assign the reads to parental

mkdir -p ${outdir}/Parental_Phase
inbam=${outdir}/Count_allele/${id}.correctCH.rmdup.allele.sorted.bam

date
if [[ ! -s ${outdir}/Parental_Phase/${id}.read_alleles.txt ]];then
$python $Phase_SNP_read \
       -b $inbam \
       -s $ref_snp \
       -o ${outdir}/Parental_Phase \
       --sampleID ${id} \
       --species $species \
       --threads $SLURM_CPUS_PER_TASK
fi
date
aa=`less ${outdir}/Parental_Phase/${id}.Read_Phased.txt.gz|head -n 5|wc -l|xargs`
#if [[ ! -s ${outdir}/Parental_Phase/${id}.Read_Phased.txt.gz ]];then
if [[ $aa == 0 ]];then
$python $Sep_parental_reads \
       -i ${outdir}/Parental_Phase/${id}.read_alleles.txt \
       -o ${outdir}/Parental_Phase \
       --sampleID ${id} \
        --species $species \
       --minAlleleNum 1 \
        --minDV 0.2 \
        --minBQ 10 \
        --threads $SLURM_CPUS_PER_TASK
fi

date
if [[ ! -s ${outdir}/Parental_Phase/${id}.cpgs_phased_M.txt ]];then
$python $Sep_CpGs \
       --phased_read_file ${outdir}/Parental_Phase/${id}.Read_Phased.txt.gz \
       --cpg_context_file ${outdir}/Count_allele/CpG_context_${id}.correctCH.rmdup.allele.sortn.txt \
        --outdir ${outdir}/Parental_Phase \
        --sampleID ${id} \
        --species $species \
        --threads $SLURM_CPUS_PER_TASK
date


awk '$1!~"^M"{split($1,aa,":");split($2,bb,":");for(cc in bb){split(bb[cc],dd,",");if(aa[4]=="GA"){dd[2]=dd[2]-1};print dd[1],dd[2],dd[3],aa[1],aa[4],aa[5]}}' ${outdir}/Parental_Phase/${id}.cpgs_phased_M.txt|grep -v Alleles|sed 's/z/0/g'|sed 's/Z/1/g'|sort -V > ${outdir}/Parental_Phase/${id}.cpgs_phased_M.CpG_context.txt
awk '$1!~"^M"{split($1,aa,":");split($2,bb,":");for(cc in bb){split(bb[cc],dd,",");if(aa[4]=="GA"){dd[2]=dd[2]-1};print dd[1],dd[2],dd[3],aa[1],aa[4],aa[5]}}' ${outdir}/Parental_Phase/${id}.cpgs_phased_P.txt|grep -v Alleles|sed 's/z/0/g'|sed 's/Z/1/g'|sort -V > ${outdir}/Parental_Phase/${id}.cpgs_phased_P.CpG_context.txt

fi

date
if [[ -s ${outdir}/Parental_Phase/${id}.Maternal.CpG_allele_m.posMerge.txt.gz ]];then
label=no
aa=`less ${outdir}/Parental_Phase/${id}.Maternal.CpG_allele_m.posMerge.txt.gz |cut -f 1|uniq|wc -l|xargs`
if [[ $aa != 21 ]] && [[ $species == "mm10" ]];then label=yes;fi
if [[ $aa != 24 ]] && [[ $species == "hg38" ]];then label=yes;fi
else
label=yes
fi

date
echo "Get the Matennal CpG Hemi values"
if [[ $label == "yes" ]];then
$python $Call_CpG_Hemi \
        -i ${outdir}/Parental_Phase/${id}.cpgs_phased_M.CpG_context.txt \
        -o ${outdir}/Parental_Phase \
        --sampleID ${id}.Maternal \
        --species $species \
        --MapQ_flag \
        --threads $SLURM_CPUS_PER_TASK
date
fi
if [[ -s ${outdir}/Parental_Phase/${id}.Paternal.CpG_allele_m.posMerge.txt.gz ]];then
label=no
aa=`less ${outdir}/Parental_Phase/${id}.Paternal.CpG_allele_m.posMerge.txt.gz |cut -f 1|uniq|wc -l|xargs`
if [[ $aa != 21 ]] && [[ $species == "mm10" ]];then label=yes;fi
if [[ $aa != 24 ]] && [[ $species == "hg38" ]];then label=yes;fi
else
label=yes
fi

date
echo "Get the Patennal CpG Hemi values"
if [[ $label == "yes" ]];then
$python $Call_CpG_Hemi \
        -i ${outdir}/Parental_Phase/${id}.cpgs_phased_P.CpG_context.txt \
        -o ${outdir}/Parental_Phase \
        --sampleID ${id}.Paternal \
        --species $species \
        --MapQ_flag \
        --threads $SLURM_CPUS_PER_TASK
date
fi

$bedtools intersect -wao -a <(less ${outdir}/Parental_Phase/${id}.Maternal.CpG_allele_m.posMerge.txt.gz|awk '{$2=$2" "$2;print $0}'|sed 's/ /\t/g') -b <(less ${outdir}/Parental_Phase/${id}.Paternal.CpG_allele_m.posMerge.txt.gz|awk '{$2=$2" "$2;print $0}'|sed 's/ /\t/g') |awk '$9!~"-1"' > ${outdir}/Parental_Phase/${id}.MP.CpG_allele_m.posMerge.txt
$bedtools intersect -wao -a <(less ${outdir}/Parental_Phase/${id}.Maternal.CpG_allele_m.sameAllele.txt.gz |awk '$3<3 {$2=$2" "$2;print $0}'|sed 's/ /\t/g') -b <(less ${outdir}/Parental_Phase/${id}.Paternal.CpG_allele_m.sameAllele.txt.gz |awk '$3<3 {$2=$2" "$2;print $0}'|sed 's/ /\t/g')  |awk '$10!~"-1"' > ${outdir}/Parental_Phase/${id}.MP.CpG_allele_m.sameAllele.txt
date

Read_Phased=`less ${outdir}/Parental_Phase/${id}.Read_Phased.txt.gz|awk 'BEGIN{N=0;M=0;P=0}{if($2=="M"){M+=1}else if($2=="P"){P+=1}else{N+=1}}END{sum=M+P+N;print sum,M,M/sum,P,P/sum,N,N/sum}'`

M_Mean_Allele_Methyl=`awk '{sum+=$3}END{print 100*sum/NR"%"}' ${outdir}/Parental_Phase/${id}.cpgs_phased_M.CpG_context.txt`
P_Mean_Allele_Methyl=`awk '{sum+=$3}END{print 100*sum/NR"%"}' ${outdir}/Parental_Phase/${id}.cpgs_phased_P.CpG_context.txt`

M_CpGs_2allele_2strand=`less ${outdir}/Parental_Phase/${id}.Maternal.CpG_allele_m.sameAllele.txt.gz |awk '$3<3 {print $1,$2}'|sort|uniq|wc -l`
P_CpGs_2allele_2strand=`less ${outdir}/Parental_Phase/${id}.Paternal.CpG_allele_m.sameAllele.txt.gz |awk '$3<3 {print $1,$2}'|sort|uniq|wc -l`
M_CpGs_2allele_2strand_Cov=`echo "scale=4;100*$M_CpGs_2allele_2strand/$CpG_num"|bc|xargs`
P_CpGs_2allele_2strand_Cov=`echo "scale=4;100*$P_CpGs_2allele_2strand/$CpG_num"|bc|xargs`
M_2allele_2strand_type=`less ${outdir}/Parental_Phase/${id}.Maternal.CpG_allele_m.sameAllele.txt.gz |awk '$3<3{if($5==0.5){Hemi+=1}else if($5==1){full_m+=1}else{full_u+=1}}END{all=Hemi+full_m+full_u;print all,Hemi,Hemi/all,full_m,full_m/all,full_u,full_u/all}'`
P_2allele_2strand_type=`less ${outdir}/Parental_Phase/${id}.Paternal.CpG_allele_m.sameAllele.txt.gz |awk '$3<3{if($5==0.5){Hemi+=1}else if($5==1){full_m+=1}else{full_u+=1}}END{all=Hemi+full_m+full_u;print all,Hemi,Hemi/all,full_m,full_m/all,full_u,full_u/all}'`

Overlap_CpGs=`wc -l ${outdir}/Parental_Phase/${id}.MP.CpG_allele_m.posMerge.txt|awk '{print $1}'|xargs`
Overlap_CpGs_Cov=`echo "scale=4;100*$Overlap_CpGs/$CpG_num"|bc|xargs`
Overlap_2strand=`wc -l ${outdir}/Parental_Phase/${id}.MP.CpG_allele_m.sameAllele.txt|awk '{print $1}'|xargs`
Overlap_2strand_Cov=`echo "scale=4;100*$Overlap_2strand/$CpG_num"|bc|xargs`

if [[ ! -f ${outdir}/stat/Parental_Phase.stat.txt ]];then
echo -e "SampleID\tspecies\tAlignedReads	MphasedReads	MReadsRatio	PphasedReads	PReadsRatio	NonPhased	NonPhasedRatio	M_Mean_Allele_Methyl	P_Mean_Allele_Methyl	Overlap_CpGs	Overlap_CpGs_cov	Overlap_CpGs_2starnd	Overlap_CpGs_2strand_cov	MCpGs_2strand_cov	PCpGs_2strand_cov	MCpGs_2strand	MCpGs_2strand_Hemi	MCpGs_HemiRatio	MCpGs_2strand_FullUM	MCpGs_2strand_FullM_Ratio	PCpGs_2strand_FullM	PCpGs_2strand_FullUM_Ratio	PCpGs_2strand	PCpGs_2strand_Hemi	PCpGs_HemiRatio	PCpGs_2strand_FullM	PCpGs_2strand_FullM_Ratio	PCpGs_2strand_FullUM	PCpGs_2strand_FullUM_Ratio" |xargs > ${outdir}/stat/Parental_Phase.stat.txt
fi

echo -e "$sampleid\t${species}\t" \
  "${Read_Phased}\t" \
  "${M_Mean_Allele_Methyl}\t" \
  "${P_Mean_Allele_Methyl}\t" \
  "${Overlap_CpGs}\t${Overlap_CpGs_Cov}%\t" \
  "${Overlap_2strand}\t${Overlap_2strand_Cov}%\t" \
  "${M_CpGs_2allele_2strand_Cov}%\t${P_CpGs_2allele_2strand_Cov}%\t" \
  "${M_2allele_2strand_type}\t${P_2allele_2strand_type}" >> ${outdir}/stat/Parental_Phase.stat.txt

chrs=(`seq 1 19` "X" "Y")
for cc in ${chrs[@]}
do
       ccM_Mean_Allele_Methyl=`awk -v cc=chr$cc '$1==cc{methyl+=$3;sum+=1}END{print 100*methyl/sum"%"}' ${outdir}/Parental_Phase/${id}.cpgs_phased_M.CpG_context.txt`
       ccP_Mean_Allele_Methyl=`awk -v cc=chr$cc '$1==cc{methyl+=$3;sum+=1}END{print 100*methyl/sum"%"}' ${outdir}/Parental_Phase/${id}.cpgs_phased_P.CpG_context.txt`
       allele2_2strand_P=`less ${outdir}/Parental_Phase/${id}.Paternal.CpG_allele_m.sameAllele.txt.gz |awk -v cc=chr${cc} '$3<3&&$1==cc {if($5==0.5){Hemi+=1}else if($5==1){full_m+=1}else{full_u+=1}}END{all=Hemi+full_m+full_u;print all,Hemi,Hemi/all,full_m,full_m/all,full_u,full_u/all}'`
       allele2_2strand_M=`less ${outdir}/Parental_Phase/${id}.Maternal.CpG_allele_m.sameAllele.txt.gz |awk -v cc=chr${cc} '$3<3&&$1==cc {if($5==0.5){Hemi+=1}else if($5==1){full_m+=1}else{full_u+=1}}END{all=Hemi+full_m+full_u;print all,Hemi,Hemi/all,full_m,full_m/all,full_u,full_u/all}'`
	echo -e "${id} $cc $ccM_Mean_Allele_Methyl $ccP_Mean_Allele_Methyl $allele2_2strand_M $allele2_2strand_P" >> ${outdir}/stat/Phased.chrom.allele.stat.txt
done

}


function allele_sites {


mkdir -p $outdir/cov_bin

### CpG with 2 strand
less ${outdir}/Count_allele/${id}.allele.CpG_allele_m.sameAllele.txt.gz |awk 'BEGIN{OFS="\t"}{print $1,$2,$2,$3,$5}' > ${outdir}/Count_allele/${id}.sameAllele.bed

$bedtools intersect -a ${SRC_DIR}/data/${species}.CpGs.bed -b ${outdir}/Count_allele/${id}.sameAllele.bed -wao > $outdir/cov_bin/${id}.sameAllele.sites.tmp

echo "Calculate the hemi abundance"
$python $countAllele2strand $outdir/cov_bin/${id}.sameAllele.sites.tmp > $outdir/cov_bin/${id}.sameAllele.sites.bed


### CpG sites
less $outdir/Count_allele/${id}.allele.CpG_allele_m.posMerge.txt.gz|awk 'BEGIN{OFS="\t"}{print $1,$2,$2,$4}' > $outdir/Count_allele/${id}.allele.CpG_allele_m.posMerge.bed

$bedtools intersect -a ${SRC_DIR}/data/${species}.CpGs.bed -b $outdir/Count_allele/${id}.allele.CpG_allele_m.posMerge.bed -wao > $outdir/Count_allele/${id}.allele.CpG.posMerge.tmp

awk 'BEGIN{while((getline t < ARGV[1]) > 0)last++;close(ARGV[1]);pos="";methyl=""}{if(pos!=$2){if(pos!=""){print $1,$2,methyl};pos=$2;if($2!=$5-1){methyl="Nan"}else{methyl=$7}}else{if($2==$5-1){methyl=$7}};if(last==FNR){print $1,$2,methyl}}' $outdir/Count_allele/${id}.allele.CpG.posMerge.tmp > $outdir/Count_allele/${id}.allele.CpG.posMerge.MeanMethyl.tmp


}


function fragment_len {
echo "Calculate the fragment length"
java -jar $picard CollectInsertSizeMetrics \
        I=$outdir/bismark_${species}/${sampleid}.${species}.bam \
        O=$outdir/bismark_${species}/${sampleid}.${species}.insert_size_metrics.txt \
        H=$outdir/bismark_${species}/${sampleid}.${species}.insert_size_histogram.pdf \
        M=0.5
}

function Align_stat {
mkdir -p $outdir/stat

Total_read_pair=`grep "Total reads processed" $outdir/trim/${sampleid}[._]R1*_trimming_report.txt | sed 's/\s\s*/ /g' | cut -d " " -f 4 | sed 's/,//g'|xargs`
if [[ $Total_read_pair == "" ]];then
Total_read_pair=`grep "Total reads processed" $outdir/trim/${sampleid}[._]R2*_trimming_report.txt | sed 's/\s\s*/ /g' | cut -d " " -f 4 | sed 's/,//g'|xargs`
fi

Total_reads=$(($Total_read_pair * 2 ))
Total_bases=`echo "$Total_reads * 150"|bc`

if [[ ! -s $outdir/trim/${sampleid}.R1_val_1_fastqc/fastqc_data.txt ]];then
unzip -n $outdir/trim/${sampleid}[._]R*_val_1_fastqc.zip -d $outdir/trim
unzip -n $outdir/trim/${sampleid}[._]R*_val_2_fastqc.zip -d $outdir/trim
fi

Clean_reads=$(($(grep "Total Sequences" $outdir/trim/${sampleid}.R1_val_1_fastqc/fastqc_data.txt |sed 's/Total Sequences\s//g'|xargs) + $(grep "Total Sequences" $outdir/trim/${sampleid}.R2_val_2_fastqc/fastqc_data.txt |sed 's/Total Sequences\s//g'|xargs)))
Clean_read_base=$(($(grep "Total written (filtered)" $outdir/trim/${sampleid}[._]R1*_trimming_report.txt |sed 's/.*:\s//g'|awk '{print $1}'|sed 's/,//g'|xargs) + $(grep "Total written (filtered)" $outdir/trim/${sampleid}[._]R2*_trimming_report.txt |sed 's/.*:\s//g'|awk '{print $1}'|sed 's/,//g'|xargs)))

Filter_reads_r=`echo "scale=2;100*(1-$Clean_reads/$Total_reads)"|bc`
Filter_base_r=`echo "scale=2;100*(1-$Clean_read_base/$Total_bases)"|bc`

reads_aligned=`grep "reads mapped:" $outdir/bismark_${species}/${sampleid}.${species}.bam.stat|cut -f 3|xargs`
base_aligned=`grep "bases mapped (cigar):" $outdir/bismark_${species}/${sampleid}.${species}.bam.stat|cut -f 3|xargs`
Align_rate=`echo "scale=2;100*$reads_aligned/$Clean_reads"|bc|xargs`

Mean_Insertion_size=`grep -A 1 MEAN_INSERT_SIZE $outdir/bismark_${species}/${sampleid}.${species}.insert_size_metrics.txt|cut -f 5|tail -n 1|xargs`

if [[ ! -s $outdir/bismark_${species}/${sampleid}.${species}.depth ]];then
$samtools depth $outdir/bismark_${species}/${sampleid}.${species}.bam > $outdir/bismark_${species}/${sampleid}.${species}.depth
fi

COVERED_REGION=$(wc -l $outdir/bismark_${species}/${sampleid}.${species}.depth|awk '{print $1}')
COVERAGE=`echo "scale=4;100*$COVERED_REGION/$genome_len"|bc|xargs`

if [[ ! -s $outdir/bismark_${species}/${sampleid}.correctCH.sorted.depth ]];then
$samtools depth $outdir/bismark_${species}/${sampleid}.correctCH.sorted.bam > $outdir/bismark_${species}/${sampleid}.correctCH.sorted.depth
fi
Correct_Covs=`wc -l $outdir/bismark_${species}/${sampleid}.correctCH.sorted.depth|awk '{print $1}'|xargs`
CHcorrect_Cov=`echo "scale=4;100*$Correct_Covs/$genome_len"|bc|xargs`
CpG_detected=`awk '$5!=""' $outdir/sites/${sampleid}.${species}.CpG_methy.txt|wc -l`
CpG_cov=`echo "scale=2;100*$CpG_detected/$CpG_num"|bc|xargs`

chrM_readsN=`$samtools idxstats $outdir/bismark_${species}/${sampleid}.${species}.bam|awk '{print $3,$1}'|sort|grep -v "_"|grep chr|awk 'BEGIN{sum=0;chrM=0}{if($2=="chrM"){chrM=$1};sum+=$1}END{print chrM/sum}'`

PERCENT_DUP=$(awk -F "\t" '{print $9}' $outdir/bismark_${species}/${id}.correctCH.rmdup.txt | grep -A 1 "PERCENT_DUPLICATION" | sed -n 2p)
Dup_rate=`echo "scale=4;$PERCENT_DUP*100"|bc|xargs`

raw_depth=`echo "scale=6;$Total_bases / $genome_len "|bc|xargs`
clean_depth=`echo "scale=6;$Clean_read_base / $genome_len "|bc|xargs`
Depth_after_align=`echo "scale=6;$base_aligned / $genome_len "|bc|xargs`
GC_1=`less $outdir/trim/${sampleid}[._]R*_val_1_fastqc/fastqc_data.txt |grep %GC|cut -f 2|xargs`
GC_2=`less $outdir/trim/${sampleid}[._]R*_val_2_fastqc/fastqc_data.txt |grep %GC|cut -f 2|xargs`
GC_content=`echo "scale=1;($GC_1+$GC_2)/2"|bc|xargs`

#if [ ! -d $outdir/Count_allele ];then
if [ ! -s $outdir/Count_allele/${id}.correctCH.rmdup.allele.sortn_splitting_report.txt ];then
bis_stat=$outdir/bismark_${species}/${id}.correctCH.rmdup.sortn_splitting_report.txt
else
bis_stat=$outdir/Count_allele/${id}.correctCH.rmdup.allele.sortn_splitting_report.txt
fi
bis_stat=$outdir/bismark_${species}/${id}.correctCH.rmdup.sortn_splitting_report.txt
#Mean_Allele_Methyl=`awk '{sum+=$3}END{print 100*sum/NR"%"}' Output_s1/Count_allele/CpG_context_${id}.allele.clean.txt`


if [[ $species == "clai" ]];then
#Mean_Allele_Methyl="Nan"
CGmethy_stat="Nan"
CHmethy_stat=`cat $outdir/bismark_${species}/C*context_${sampleid}.correctCH.rmdup.sortn.txt|cut -f 3-|awk '$2>18 && $2<68'|awk 'BEGIN{sum=0;umethy=0}{sum+=1;if($3~"X|H|Z"){methy+=1}}END{print sum,methy/sum}'|sed 's/ /\t/g'`

else
CGmethy_stat=`grep "CpG" $bis_stat |sed 's/.*:\t//g'|awk -v ORS="\t" '{print}'|awk '{CpG=$1+$2; cpg_r=100*$1/CpG;print CpG,cpg_r"%"}'|sed 's/ /\t/g'`
CHmethy_stat=`grep "CH[HG]" $bis_stat |sed 's/.*:\t//g'|awk -v ORS="\t" '{print}'|awk '{CHG=$1+$3;CHH=$2+$4;chg_r=100*$1/CHG;chh_r=100*$2/CHH;print chg_r"%",chh_r"%"}'|sed 's/ /\t/g'`

fi


if [[ ! -f ${outdir}/stat/seq_stat.${species}.txt ]];then
	echo -e "SampleID	Species	Total_reads	Clean_Reads	Reads_Clean_rate	Base_Clean_rate	reads_aligned	Align_rate	Genome_Cov	CHcorrected_GenomeCov	CpG_detected	CpG_coverage	chrM_readsN	MeanFragment_len	DupRate	Raw_depth	Clean_depth	Align_depth	GC_content	Seqed_CpGs	Bismark_mCpG_ConversionRate	mCHG_ratio	mCHH_ratio" > ${outdir}/stat/seq_stat.${species}.txt
fi

echo -e "$sampleid\t${species}\t" \
  "${Total_reads}\t" \
  "${Clean_reads}\t" \
  "${Filter_reads_r}%\t" \
  "${Filter_base_r}%\t" \
  "${reads_aligned}\t" \
  "${Align_rate}%\t" \
  "${COVERAGE}%\t" \
  "${CHcorrect_Cov}%\t" \
  "${CpG_detected}\t" \
  "${CpG_cov}%\t" \
  "${chrM_readsN}\t" \
  "${Mean_Insertion_size}\t" \
  "${Dup_rate}%\t" \
  "${raw_depth}\t" \
  "${clean_depth}\t" \
  "${Depth_after_align}\t" \
  "${GC_content}%\t" \
  "${CGmethy_stat}\t" \
  "${CHmethy_stat}\t" \
  >> $outdir/stat/seq_stat.${species}.txt

}


set -x

trim_array

bismark_call

fragment_len

Remove_mC_filledGap PE

if [[ ! $CallHemi ]];then
    call_methyl
else
    count_allele

    allele_sites

    if [[ ! $Parental_Phase ]];then
        Parental_Phase
    fi
fi

Align_stat


echo "All pipeline Finished"
date

set +x
