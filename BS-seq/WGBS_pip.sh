#!/usr/bin/bash
#SBATCH -D ./
#SBATCH -p compute
#SBATCH --mem=10G
##SBATCH --exclude=node07
#SBATCH --cpus-per-task=2
#SBATCH -o %J.%j.pip.log
#SBATCH -e %J.%j.pip.log

SRC_DIR=/share/home/baiyl/merlot/TAPS-merlot/workflow/RUN_src/src_1.8

trim_galore=/share/newdata4/baiyl/basic/tools/TrimGalore-0.4.5/trim_galore
python='/share/newdata4/baiyl/basic/tools/anaconda2/bin/python'

bismark="/share/newdata4/baiyl/basic/tools/Bismark_v0.19.0/bismark"
bowtie2="/share/newdata4/baiyl/basic/tools/anaconda2/bin"
picard=/datum/data2/baiyl/bin/src/picard-tools-2.4.1/picard.jar
bismark_methylation_extractor="/share/newdata4/baiyl/basic/tools/Bismark_v0.19.0/bismark_methylation_extractor"
samtools='/share/newdata4/baiyl/basic/tools/anaconda2/bin/samtools'

CpG_merge=/share/home/baiyl/merlot/TAPS-merlot/workflow/RUN_src/src_1.8/scripts/bismark_CpGsites_merge.i.pl
genome_CpG_py=/share/home/baiyl/merlot/TAPS-merlot/workflow/RUN_src/src_1.8/call_mod/genome_CpG_num.py

hg38_bis=/share/home/baiyl/basic/database/human/hg38/Bismark_hg38
mm10_bis=/share/home/baiyl/basic/database/mouse/mm10/Index/Bismark_mm10
lambda_bis=/share/home/baiyl/basic/database/NEB_lambda/lambda_Bismark

hg38_ref_fa=/share/newdata4/baiyl/basic/database/human/hg38/Index/bwa/hg38.fa
mm10_ref_fa=/share/home/baiyl/basic/database/Illumina_iGENOMES/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa
lambda_ref_fa=/share/home/baiyl/database/NEB_lambda/BWAIndex/genome.fa

declare -A Ref_fa
Ref_fa["hg38"]=$hg38_ref_fa
Ref_fa["mm10"]=$mm10_ref_fa
Ref_fa["lambda"]=$lambda_ref_fa

declare -A Ref_bis
Ref_bis["hg38"]=$hg38_bis
Ref_bis["mm10"]=$mm10_bis
Ref_bis["lambda"]=$lambda_bis

set -x
#----------------------------- ARGS -----------------------------

ddate=`date -R|awk '{print $2$3}'`

outdir=$1
species=$2
Read1=$3
Read2=$4
read_len=$5
run_SE=PE

methy_threshold=0.5

cores=`expr $SLURM_CPUS_PER_TASK / 5`


ref_fa=${Ref_fa[$species]}
ref_bis=${Ref_bis[$species]}

if [ ! -s $Read1 ];then
        echo $Read1 "Read not exist"
        exit
fi
if [ ! -s $Read2 ];then
        echo $Read2 "Read not exist"
fi

id=`echo $Read1|sed 's/.*\///g'|sed 's/.R1[._].*//g'|xargs`

echo "Start processing $id with $Read1, $Read2 in $species"

##---------------------------- FUNCTIONS --------------------------

mkdir -p $outdir/stat
if [[ ! -f $SRC_DIR/data/genome.CpG_num.txt ]];then
        genome_CpG=`$python $genome_CpG_py $ref_fa CG $SRC_DIR/data/${species}.CpGs.txt|sed 's/.*sites//g'|xargs`
        echo -e "${species}\t$genome_CpG" >> $SRC_DIR/data/genome.CpG_num.txt
else
        genome_CpG=`grep -w ${species} $SRC_DIR/data/genome.CpG_num.txt|awk '{print $2}'|uniq|xargs`
        if [[ $genome_CpG == "" ]];then
                genome_CpG=`$python $genome_CpG_py $ref_fa CG $SRC_DIR/data/${species}.CpGs.txt|sed 's/.*sites//g'|xargs`
                echo -e "${species}\t$genome_CpG" >> $SRC_DIR/data/genome.CpG_num.txt

        fi
fi

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



function :q
:{
mkdir -p $outdir/trim
if [[ -f $outdir/trim/${id}.R1_val_1.fq.gz ]] && [[ -f $outdir/trim/${id}.R2_val_2.fq.gz ]];then
        echo "Trim already finished"
        return
else
        rm -fr $outdir/trim/${id}*
fi

$trim_galore --fastqc --paired \
	--phred33 \
        --retain_unpaired \
	--clip_R1 9 --clip_R2 9 \
	--output_dir $outdir/trim \
	$Read1 $Read2
}

function CpG_stat {
CpG_file=$1
CpG_out=$2

PE_CpG_stat=(`perl $CpG_merge $CpG_file $CpG_out $methy_threshold|xargs`)

Total_CpG_seqed=${PE_CpG_stat[0]}
Total_mCpG_seqed=${PE_CpG_stat[1]}
CpG_methy_depth=${PE_CpG_stat[2]}
Total_CpG_sites=${PE_CpG_stat[3]}
Total_mCpG_sites=${PE_CpG_stat[4]}
CpG_methy_sites=${PE_CpG_stat[5]}
CpG_depth=`echo "scale=3;$Total_CpG_seqed/$genome_CpG"|bc|xargs`
CpG_cover=`echo "scale=3;100*$Total_CpG_sites/$genome_CpG"|bc|xargs`

arr_out=($Total_CpG_seqed $CpG_depth $Total_CpG_sites ${CpG_cover}% $Total_mCpG_seqed $CpG_methy_depth $Total_mCpG_sites $CpG_methy_sites)
outff=`echo $CpG_out|sed 's/txt/stat.txt/g'`
echo ${arr_out[*]} > $outff
}


function se_bismark {
unmap_read=$1
read_t=$2
$bismark --multicore $cores \
        --fastq --non_directional --unmapped \
        --nucleotide_coverage \
        --path_to_bowtie $bowtie2 --bowtie2 \
        --genome_folder $ref_bis \
        --output_dir $outdir/bismark \
        --temp_dir $outdir/bismark/temp_unmap_${id} \
        $unmap_read

rm -fr $outdir/bismark/temp_unmap_${id}
#

if [[ $read_t == "r1" ]];then
new_id=${id}.${species}.unmap_r1
else
new_id=${id}.${species}.unmap_r2
fi
outbam=`ls $outdir/bismark/${new_id}_bismark_bt2.bam|xargs`


$samtools sort -@ $SLURM_CPUS_PER_TASK $outbam > $outdir/bismark/${new_id}.sort.bam
$samtools index $outdir/bismark/${new_id}.sort.bam
$samtools stats $outdir/bismark/${new_id}.sort.bam > $outdir/bismark/${new_id}.sort.bam.stat

java -jar $picard MarkDuplicates \
        I=$outdir/bismark/${new_id}.sort.bam \
        O=$outdir/bismark/${new_id}.rmdup.bam \
        M=$outdir/bismark/${new_id}.rmdup.txt \
        REMOVE_DUPLICATES=true
$samtools index $outdir/bismark/${new_id}.rmdup.bam

$samtools sort -n -@ $SLURM_CPUS_PER_TASK $outdir/bismark/${new_id}.rmdup.bam > $outdir/bismark/${new_id}.sorted.bam
$bismark_methylation_extractor -s \
        --multicore $cores \
        --comprehensive \
        --no_overlap \
        --bedGraph \
        --counts \
        --buffer_size 20G \
        --report \
        --cytosine_report \
        --genome_folder $ref_bis \
        -o $outdir/bismark \
	$outdir/bismark/${new_id}.sorted.bam

rm -fr $outdir/bismark/temp_unmap_${id}

CpG_file=$outdir/bismark/${new_id}.sorted.CpG_report.txt.CpG_report.txt

awk '$4+$5>0' $outdir/bismark/${new_id}.sorted.CpG_report.txt.CpG_report.txt > $outdir/bismark/${new_id}.cut.CpG_report.txt
CpG_file=$outdir/bismark/${new_id}.cut.CpG_report.txt

CpG_stat $CpG_file $outdir/bismark/${new_id}.CpG.txt
}

function bismark_call {
trim_Read1=$outdir/trim/${id}.R1_val_1.fq.gz
trim_Read2=$outdir/trim/${id}.R2_val_2.fq.gz
mkdir -p $outdir/bismark_${species}

$bismark --multicore $cores \
        --fastq --non_directional --unmapped \
        --nucleotide_coverage \
        --path_to_bowtie $bowtie2 --bowtie2 \
        --genome_folder $ref_bis \
        --output_dir $outdir/bismark_${species} \
        --temp_dir $outdir/bismark_${species}/temp_bismark_${id} \
	-1 $trim_Read1 -2 $trim_Read2

outbam=`ls $outdir/bismark_${species}/${id}*_bismark_bt2_pe.bam|xargs`

$samtools view -bh -q 1 -F 4 $outbam |$samtools sort -@ $SLURM_CPUS_PER_TASK > $outdir/bismark_${species}/${id}.${species}.bam

$samtools index $outdir/bismark_${species}/${id}.${species}.bam
$samtools stats $outdir/bismark_${species}/${id}.${species}.bam > $outdir/bismark_${species}/${id}.${species}.bam.stat
$samtools depth $outdir/bismark_${species}/${id}.${species}.bam > $outdir/bismark_${species}/${id}.${species}.depth


java -jar $picard MarkDuplicates \
        I=$outdir/bismark_${species}/${id}.${species}.bam \
        O=$outdir/bismark_${species}/${id}.${species}.rmdup.bam \
        M=$outdir/bismark_${species}/${id}.${species}.rmdup.txt \
        REMOVE_DUPLICATES=true
$samtools index $outdir/bismark_${species}/${id}.${species}.rmdup.bam

$samtools sort -n -@ $SLURM_CPUS_PER_TASK $outdir/bismark_${species}/${id}.${species}.rmdup.bam > $outdir/bismark_${species}/${id}.${species}.sorted.bam

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
   -o $outdir/bismark_${species} \
   $outdir/bismark_${species}/${id}.${species}.sorted.bam

rm -fr $outdir/bismark_${species}/temp_bismark_${id}


awk '$4+$5>0' $outdir/bismark_${species}/${id}.${species}.sort*CpG_report.txt.CpG_report.txt > $outdir/bismark_${species}/${id}.${species}.cut.CpG_report.txt
CpG_file=$outdir/bismark_${species}/${id}.${species}.cut.CpG_report.txt


CpG_stat $CpG_file $outdir/bismark_${species}/${id}.${species}.CpG.txt

CpG_f=`ls ${outdir}/bismark_${species}/${id}.*CpG_report.txt.CpG_report.txt|xargs`

awk 'BEGIN{chr="Chr";pos="Pos";methyl="Methyl";cov="UMethyl"}{OFS="\t";if($3=="-"){$2=$2-1};if(pos==$2){methyl+=$4;cov+=$5}else{sum=methyl+cov;if(sum==0){rr=""}else{rr=methyl/sum};print chr,pos,methyl,cov,rr;chr=$1;pos=$2;methyl=$4;cov=$5}}' $CpG_f  > $outdir/bismark_${species}/${id}.${species}.CpG_methy.txt

grep -v Pos $outdir/bismark_${species}/${id}.${species}.CpG_methy.txt|grep -v "_" |sort -V > $outdir/bismark_${species}/${id}.${species}.CpG_methy.tmp
cat <(echo -e "${id}.cov\t${id}.ratio") <(cut -f 5 $outdir/bismark_${species}/${id}.${species}.CpG_methy.tmp) > $outdir/bismark_${species}/${id}.${species}.CpG_methy.csv


}

function fragment_len {
java -jar $picard CollectInsertSizeMetrics \
        I=$outdir/bismark_${species}/${id}.${species}.bam \
        O=$outdir/bismark_${species}/${id}.${species}.insert_size_metrics.txt \
        H=$outdir/bismark_${species}/${id}.${species}.insert_size_histogram.pdf \
        M=0.5
}

function se_remap {
mv $outdir/bismark/${id}.R1_val_1.fq.gz_unmapped_reads_1.fq.gz $outdir/bismark/${id}.${species}.unmap_r1.fq.gz
mv $outdir/bismark/${id}.R2_val_2.fq.gz_unmapped_reads_2.fq.gz $outdir/bismark/${id}.${species}.unmap_r2.fq.gz

unmap_r1=$outdir/bismark/${id}.${species}.unmap_r1.fq.gz
unmap_r2=$outdir/bismark/${id}.${species}.unmap_r2.fq.gz

se_bismark $unmap_r1 r1
se_bismark $unmap_r2 r2

$samtools merge -@ $SLURM_CPUS_PER_TASK -f \
        $outdir/bismark/${id}.${species}.merged.bam \
        $outdir/bismark/${id}.${species}.sorted.bam \
        $outdir/bismark/${id}.${species}.unmap_r1.sorted.bam \
        $outdir/bismark/${id}.${species}.unmap_r2.sorted.bam

$samtools sort -@ $SLURM_CPUS_PER_TASK $outdir/bismark/${id}.${species}.merged.bam > $outdir/bismark/${id}.${species}.merged.sort.bam
$samtools index $outdir/bismark/${id}.${species}.merged.sort.bam

cat $outdir/bismark/${id}.${species}.CpG.txt \
        $outdir/bismark/${id}.${species}.unmap_r1.CpG.txt \
        $outdir/bismark/${id}.${species}.unmap_r2.CpG.txt \
        > $outdir/bismark/${id}.${species}.merged.CpG.txt

}

function methy_stat {

CHG_methy_stat=`awk 'BEGIN{mm=0;sum=0}{if($5~"[XHZ]"){mm+=1};sum+=1}END{rr=mm/sum;print sum,mm,rr}' $outdir/bismark_${species}/CHG_context_${id}.${species}.sorted.txt |sed 's/ /\t/g'|xargs`
CHH_methy_stat=`awk 'BEGIN{mm=0;sum=0}{if($5~"[XHZ]"){mm+=1};sum+=1}END{rr=mm/sum;print sum,mm,rr}' $outdir/bismark_${species}/CHH_context_${id}.${species}.sorted.txt |sed 's/ /\t/g'|xargs`


PE_CpG_stat=`less $outdir/bismark_${species}/${id}.${species}.CpG.stat.txt |xargs|sed 's/ /\t/g'`


if [[ $run_SE == "SE" ]];then

r1_CHG_stat=`awk 'BEGIN{mm=0;sum=0}{if($5~"[XHZ]"){mm+=1};sum+=1}END{rr=mm/sum;print sum,mm,rr}' $outdir/bismark_${species}/CHG_context_${id}.${species}.unmap_r1.sorted.txt |sed 's/ /\t/g'|xargs`
r1_CHH_stat=`awk 'BEGIN{mm=0;sum=0}{if($5~"[XHZ]"){mm+=1};sum+=1}END{rr=mm/sum;print sum,mm,rr}' $outdir/bismark_${species}/CHH_context_${id}.${species}.unmap_r1.sorted.txt |sed 's/ /\t/g'|xargs`

r2_CHG_stat=`awk 'BEGIN{mm=0;sum=0}{if($5~"[XHZ]"){mm+=1};sum+=1}END{rr=mm/sum;print sum,mm,rr}' $outdir/bismark_${species}/CHG_context_${id}.${species}.unmap_r2.sorted.txt |sed 's/ /\t/g'|xargs`
r2_CHH_stat=`awk 'BEGIN{mm=0;sum=0}{if($5~"[XHZ]"){mm+=1};sum+=1}END{rr=mm/sum;print sum,mm,rr}' $outdir/bismark_${species}/CHH_context_${id}.${species}.unmap_r2.sorted.txt |sed 's/ /\t/g'|xargs`

r1_CpG_stat=(`less $outdir/bismark_${species}/${id}.${species}.unmap_r1.CpG.stat.txt |xargs`)
r2_CpG_stat=(`less $outdir/bismark_${species}/${id}.${species}.unmap_r2.CpG.stat.txt |xargs`)

cat $outdir/bismark_${species}/${id}.${species}*.cut.CpG_report.txt > $outdir/bismark_${species}/${id}.${species}.merged.cut.CpG_report.txt
CpG_file=$outdir/bismark_${species}/${id}.${species}.merged.cut.CpG_report.txt
CpG_out=$outdir/bismark_${species}/${id}.${species}.merged.CpG.txt

PE_SE_CpG_stat=(`perl $CpG_merge $CpG_file $CpG_out $methy_threshold|xargs`)

Total_CpG_seqed=${PE_SE_CpG_stat[0]}
Total_mCpG_seqed=${PE_SE_CpG_stat[1]}
CpG_methy_depth=${PE_SE_CpG_stat[2]}
Total_CpG_sites=${PE_SE_CpG_stat[3]}
Total_mCpG_sites=${PE_SE_CpG_stat[4]}
CpG_methy_sites=${PE_SE_CpG_stat[5]}

new_CpG_depth=`echo "scale=4;$Total_CpG_seqed/$genome_CpG"|bc|xargs`
new_CpG_cover=`echo "scale=4;100*$Total_CpG_sites/$genome_CpG"|bc|xargs`



echo -e "${id}\t${species}\t" \
"${PE_CpG_stat}\t" \
"${Total_CpG_seqed}\t" \
"${new_CpG_depth}\t" \
"${Total_CpG_sites}\t" \
"${new_CpG_cover}%\t" \
"${CpG_methy_depth}\t" \
"${CpG_methy_sites}\t" \
"${CHG_methy_stat}\t" \
"${r1_CHG_stat}\t" \
"${r2_CHG_stat}\t" \
"${CHH_methy_stat}\t" \
"${r1_CHH_stat}\t" \
"${r2_CHH_stat}\t" \
>> $outdir/stat/${species}.methy_stat.SEremap.${ddate}.txt
else

echo -e "${id}\t${species}\t" \
"${PE_CpG_stat}\t" \
"${CHG_methy_stat}\t" \
"${CHH_methy_stat}" \
>> $outdir/stat/${species}.methy_stat.${ddate}.txt


fi


}

function Align_stat {
Total_read_pair=`grep "Total reads processed" $outdir/trim/${id}[._]R1*_trimming_report.txt | sed 's/\s\s*/ /g' | cut -d " " -f 4 | sed 's/,//g'|xargs`
if [[ $Total_read_pair == "" ]];then
  Total_read_pair=`grep "Total reads processed" $outdir/trim/${id}[._]R2*_trimming_report.txt | sed 's/\s\s*/ /g' | cut -d " " -f 4 | sed 's/,//g'|xargs`
fi
Total_reads=$(($Total_read_pair * 2 ))
Total_bases=`echo "$Total_reads * 150"|bc`

trim_Read1=$outdir/trim/${id}.R1_val_1.fq.gz
trim_Read2=$outdir/trim/${id}.R2_val_2.fq.gz

unzip -n $outdir/trim/${sampleid}[._]R*_val_1_fastqc.zip -d $outdir/trim
unzip -n $outdir/trim/${sampleid}[._]R*_val_2_fastqc.zip -d $outdir/trim

Clean_reads=$(($(grep "Total Sequences" $outdir/trim/${id}.R1_val_1_fastqc/fastqc_data.txt |sed 's/Total Sequences\s//g'|xargs) + $(grep "Total Sequences" $outdir/trim/${id}.R2_val_2_fastqc/fastqc_data.txt |sed 's/Total Sequences\s//g'|xargs)))
Clean_read_base=$(($(grep "Total written (filtered)" $outdir/trim/${id}[._]R1*_trimming_report.txt |sed 's/.*:\s//g'|awk '{print $1}'|sed 's/,//g'|xargs) + $(grep "Total written (filtered)" $outdir/trim/${id}[._]R2*_trimming_report.txt |sed 's/.*:\s//g'|awk '{print $1}'|sed 's/,//g'|xargs)))

Filter_reads_r=`echo "scale=2;100*(1-$Clean_reads/$Total_reads)"|bc`
Filter_base_r=`echo "scale=2;100*(1-$Clean_read_base/$Total_bases)"|bc`

Read1_Base_after_adapter_trim=`zcat $Read1|awk 'BEGIN{aa=2;sum=0}{if(NR==aa){aa+=4;sum+=length($0);}}END{print sum}'|xargs`
Read2_Base_after_adapter_trim=`zcat $Read2|awk 'BEGIN{aa=2;sum=0}{if(NR==aa){aa+=4;sum+=length($0);}}END{print sum}'|xargs`
Base_after_adapter_trim=`echo "$Read1_Base_after_adapter_trim + $Read2_Base_after_adapter_trim"|bc|xargs`

reads_aligned=`grep "reads mapped:" $outdir/bismark_${species}/${id}.${species}.bam.stat|cut -f 3|xargs`
base_aligned=`grep "bases mapped (cigar):" $outdir/bismark_${species}/${id}.${species}.bam.stat|cut -f 3|xargs`
Align_rate=`echo "scale=2;100*$reads_aligned/$Clean_reads"|bc|xargs`

Mean_Insertion_size=`grep -A 1 MEAN_INSERT_SIZE $outdir/bismark_${species}/${id}.${species}.insert_size_metrics.txt|cut -f 5|tail -n 1|xargs`
Median_Insertion_size=`grep -A 1 MEDIAN_INSERT_SIZE $outdir/bismark_${species}/${id}.${species}.insert_size_metrics.txt|cut -f 1|tail -n 1|xargs`

chrM_readsN=`$samtools idxstats $outdir/bismark_${species}/${id}.${species}.bam|awk '{print $3,$1}'|sort|grep -v "_"|grep chr|awk 'BEGIN{sum=0;chrM=0}{if($2=="chrM"){chrM=$1};sum+=$1}END{print chrM/sum}'`

COVERED_REGION=$($samtools depth $outdir/bismark_${species}/${id}.${species}.bam |wc -l )
COVERAGE=`echo "scale=4;100*$COVERED_REGION/$genome_len"|bc|xargs`
CpG_detected=`awk '$5!=""' $outdir/bismark_${species}/${id}.${species}.CpG_methy.txt|wc -l`
CpG_cov=`echo "scale=2;100*$CpG_detected/$genome_CpG"|bc|xargs`

PERCENT_DUP=$(awk -F "\t" '{print $9}' $outdir/bismark_${species}/${id}.${species}.rmdup.txt| grep -A 1 "PERCENT_DUPLICATION" | sed -n 2p)
Dup_rate=`echo "scale=4;$PERCENT_DUP*100"|bc|xargs`

raw_depth=`echo "scale=6;$Total_bases / $genome_len "|bc|xargs`
#DEPTH=`echo "scale=6;$Clean_read_base / $genome_len "|bc|xargs`
Depth_after_trim=`echo "scale=6;$Base_after_adapter_trim / $genome_len "|bc|xargs`
Depth_after_align=`echo "scale=6;$base_aligned / $genome_len "|bc|xargs`

unzip -n $outdir/trim/${id}[._]R*_val_1_fastqc.zip -d $outdir/trim
unzip -n $outdir/trim/${id}[._]R*_val_2_fastqc.zip -d $outdir/trim
GC_1=`less $outdir/trim/${id}[._]R*_val_1_fastqc/fastqc_data.txt |grep %GC|cut -f 2|xargs`
GC_2=`less $outdir/trim/${id}[._]R*_val_2_fastqc/fastqc_data.txt |grep %GC|cut -f 2|xargs`
GC_content=`echo "scale=1;($GC_1+$GC_2)/2"|bc|xargs`

bis_stat=$outdir/bismark_${species}/${id}.*splitting_report.txt
CGmethy_stat=`grep "CpG" $bis_stat |sed 's/.*:\t//g'|awk -v ORS="\t" '{print}'|awk '{CpG=$1+$2; cpg_r=100*$1/CpG;print CpG,cpg_r"%"}'|sed 's/ /\t/g'`
CHmethy_stat=`grep "CH[HG]" $bis_stat |sed 's/.*:\t//g'|awk -v ORS="\t" '{print}'|awk '{CHG=$1+$3;CHH=$2+$4;chg_r=100*$1/CHG;chh_r=100*$2/CHH;print chg_r"%",chh_r"%"}'|sed 's/ /\t/g'`


if [[ $run_SE == "SE" ]];then
        r1_reAlign=`grep "reads mapped:" $outdir/bismark_${species}/${id}.${species}.unmap_r1.sort.bam.stat|cut -f 3|xargs`
        r1_baseAlign=`grep "bases mapped (cigar):" $outdir/bismark_${species}/${id}.${species}.unmap_r1.sort.bam.stat|cut -f 3|xargs`
        r2_reAlign=`grep "reads mapped:" $outdir/bismark_${species}/${id}.${species}.unmap_r2.sort.bam.stat|cut -f 3|xargs`
        r2_baseAlign=`grep "bases mapped (cigar):" $outdir/bismark_${species}/${id}.${species}.unmap_r2.sort.bam.stat|cut -f 3|xargs`

        new_align_rate=`echo "scale=2;100*($r1_reAlign+$r2_reAlign+$reads_aligned)/$Clean_reads"|bc|xargs`

        new_covs=$($samtools depth $outdir/bismark_${species}/${id}.${species}.merged.sort.bam|wc -l)
        new_coverage=`echo "scale=4;100*$new_covs/$genome_len"|bc|xargs`

        r1_dup=$(awk -F "\t" '{print $9}' $outdir/bismark_${species}/${id}.${species}.unmap_r1.rmdup.txt| grep -A 1 "PERCENT_DUPLICATION" | sed -n 2p)
  r2_dup=$(awk -F "\t" '{print $9}' $outdir/bismark_${species}/${id}.${species}.unmap_r2.rmdup.txt| grep -A 1 "PERCENT_DUPLICATION" | sed -n 2p)
  r1_dup_rate=`echo "scale=4;$r1_dup*100"|bc|xargs`
  r2_dup_rate=`echo "scale=4;$r2_dup*100"|bc|xargs`

  new_align_depth=`echo "scale=6;($r1_baseAlign + $r2_baseAlign + $base_aligned)/ $genome_len "|bc|xargs`


  echo -e "$id\t${species}\t" \
  "${Total_reads}\t" \
  "${Clean_reads}\t" \
  "${Filter_reads_r}%\t" \
  "${Filter_base_r}%\t" \
  "${reads_aligned}\t${Align_rate}%\t" \
  "${r1_reAlign}\t${r2_reAlign}\t${new_align_rate}%\t" \
  "${COVERAGE}%\t" \
  "${CpG_cov}%\t" \
  "${chrM_readsN}\t" \
  "${Mean_Insertion_size}\t" \
  "${Dup_rate}%\t${r1_dup_rate}%\t${r2_dup_rate}%\t" \
  "${raw_depth}\t${Depth_after_trim}\t${Depth_after_align}\t" \
  "${new_align_depth}\t" \
  "${GC_content}%\t" \
  "${CGmethy_stat}\t" \
  "${CHmethy_stat}\t" >> $outdir/stat/seq_stat.${species}.txt


else

if [[ ! -f ${outdir}/stat/seq_stat.${species}.txt ]];then
	echo -e "SampleID       Species Total_reads     Clean_Reads     Reads_Clean_rate        Base_Clean_rate reads_aligned   Align_rate      Genome_Cov      CpG_coverage    MeanFragment_len        DupRate Raw_depth       Clean_depth     Align_depth     CpG_seqed	Bismark_mCpG_ConversionRate     mCHG_ratio      mCHH_ratio" > ${outdir}/stat/seq_stat.${species}.txt
fi
  echo -e "$id\t${species}\t" \
  "${Total_reads}\t" \
  "${Clean_reads}\t" \
  "${Filter_reads_r}%\t" \
  "${Filter_base_r}%\t" \
  "${reads_aligned}\t${Align_rate}%\t" \
  "${COVERAGE}%\t" \
  "${CpG_cov}%\t" \
  "${chrM_readsN}\t" \
  "${Mean_Insertion_size}\t" \
  "${Dup_rate}%\t" \
  "${raw_depth}\t${Depth_after_trim}\t${Depth_after_align}\t" \
  "${GC_content}%" \
  "${CGmethy_stat}\t" \
  "${CHmethy_stat}\t" \
  >> $outdir/stat/seq_stat.${species}.txt

fi

}

trim_reads

bismark_call

fragment_len

#se_remap

#methy_stat

Align_stat








set +x
