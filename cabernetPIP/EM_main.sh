#!/usr/bin/bash

indir=""
outdir="."	
prefix="*"	
sep_type=""	 # "meta_mark", "meta_sep"
meta_type=meta20	# "meta20"
tag_file=""	# 
trim_type=basic_trim	# "basic_trim", "nextera_trim"
thread=10
read_len=20
#clip_len=
species=hg38

ddate=`date -R|awk '{print $2$3}'`

src_dir=/share/home/baiyl/merlot/TAPS-merlot/workflow/RUN_src/src_2.0
sep_sh=$src_dir/scripts/meta_process.sh


temp=`getopt -o i:o:s:p:m:t: --long indir:,outdir:,sep_type:,prefix:,meta_type:,thread:,tag_file:,trim_type:,read_len:,clip_len:,species:,seq_type:,partition: -- "$@"`
if [ $? != 0 ] ; then echo "terminating..." >&2 ; exit 1 ; fi
eval set -- "$temp" 

while true ; do
        case "$1" in
                -i|--indir) echo "Input dir is $2"; indir=$2; shift 2;;
                -o|--outdir) echo "Output dir is $2" ;outdir=$2; shift 2;;
                -s|--sep_type) sep_type=$2; echo "Separate the reads with barcodes. Default: meta_mark" ; shift 2;;
                -p|--prefix) prefix=$2; echo "Reads prefix is $prefix" ; shift 2;;
                -m|--meta_type) meta_type=$2; echo "meta type is $meta_type"; shift 2 ;;
		--tag_file) tag_file=$2; echo "Tag file is $tag_file"; shift 2 ;;
		--trim_type) trim_type=$2;echo "option b, argument \`$2'" ; shift 2;;
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

if [[ ! -n $indir ]];
then
	echo "Please give input dir"
	exit
fi
if [[ ! -n $clip_len ]];
then
	clip_len="Na"
fi
if [[ ! -n $tag_file ]];then
	tag_file="Na"
fi
if [[ ! -n $seq_type ]];then
        seq_type="nextera"
fi




if [[ ! -n $prefix ]];then prefix="*";fi

#set -x

mkdir -p $outdir
mkdir -p $outdir/run_logs
LOG_DIR=$outdir/run_logs

function sep_main {
if [[ -n $sep_type ]]
then
	echo "Start seperate reads by barcode"

	mkdir -p $outdir/sep_cell
		
	find $indir -maxdepth 1 -name "${prefix}*.f*q.gz" -print |sort --version-sort | xargs -L2 |column -t |while read reads
do
        sampleid=`echo $reads|sed 's/.*\///g'|sed 's/[._]R[12].*//g'|xargs`

        arr=(${reads// / })
        read1=${arr[0]}
        read2=${arr[1]}

	if [[ ! $read1 =~ "R1" ]];then
                echo $read1 "Read1 not right"
                exit
        fi
        if [[ ! $read2 =~ "R2" ]];then
                echo $read2 "Read2 not right"
                exit
        fi

        sample_size=$(stat -c "%s" $read1)
        sample_size=`echo "($sample_size * 2 + 2000000000)/1000000"|bc|xargs`

	echo $sampleid, $sample_size


sbatch << RUN
#!/usr/bin/bash
#SBATCH -J s_$sampleid 
#SBATCH -o $LOG_DIR/${sampleid}.sep.${ddate}.log
#SBATCH -e $LOG_DIR/${sampleid}.sep.${ddate}.log
#SBATCH --mem=$sample_size
#SBATCH --exclude=node08,node07
#SBATCH --cpus-per-task=$thread

if [[ $tag_file != "Na" ]];then
$sep_sh --read1 $read1 \
        --read2 $read2 \
        --sep_type $sep_type  \
        --prefix $sampleid \
        --meta_type $meta_type \
        --thread $thread \
	--outdir $outdir \
        --tag_file $tag_file

else
$sep_sh --read1 $read1 \
        --read2 $read2 \
        --sep_type $sep_type  \
        --prefix $sampleid \
        --meta_type $meta_type \
        --thread $thread \
	--outdir $outdir
fi

RUN

done

fi
}
function stat_head {
mkdir -p $outdir/stat

echo -e "SampleID\tSpecies\tTotal_reads\t" \
 "Clean_Reads\tReads_Clean_rate\tBase_Clean_rate\t" \
 "Align_reads\tAlign_rate\t" \
 "Genome_Cov\tCHcorrected_GenomeCov\t" \
 "MeanFragment_len\tDupRate\t" \
 "Raw_depth\tClean_depth\tAlign_depth\t" \
 "GC_content\t" \
 "Total_CpGs_dep\tmCpGs_dep\tmCpGs_dep_ratio\t" \
 "Total_CpGs_detected\tCpG_coverage\t" \
 "mCpGs_detected(sites ratio>0)\tmCpGs_ratio(sites ratio>0)\t" \
 "mCpGs_detected(sites ratio>0.5)\tmCpGs_ratio(sites ratio>0.5)\t" \
 "Total_CHGs\tmCHGs\tmCHG_ratio\t" \
 "Total_CHHs\tmCHHs\tmCHH_ratio\t" \
 > $outdir/stat/seq_stat.${species}.txt

}


function pip_run {
pip_sh=$src_dir/scripts/pip_EM.sh

if [[ $sep_type =~ "meta" ]]
then
	indir=$outdir/sep_cell
fi

find $indir -maxdepth 1 -name "${prefix}*.f*q.gz" -print |sort --version-sort | xargs -L2 |column -t |while read reads
do
        sampleid=`echo $reads|sed 's/.*\///g'|sed 's/.R[12].*//g'|xargs`
        arr=(${reads// / })
        read1=${arr[0]}
        read2=${arr[1]}
	if [[ ! $read1 =~ "R1" ]];then
                echo $read1 "Read1 not exist"
                exit
        fi
	if [[ ! $read2 =~ "R2" ]];then
                echo $read2 "Read2 not exist"
                exit
        fi

	sample_size=$(stat -c "%s" $read1)
        sample_size=`echo "($sample_size * 2 + 1000000000)/100000"|bc|xargs`

	echo "   ---"$sampleid memery: ${sample_size}M

	sbatch -J $sampleid \
        -o $LOG_DIR/${sampleid}.pip.${species}.p.${ddate}.log \
        -e $LOG_DIR/${sampleid}.pip.${species}.p.${ddate}.log \
        --mem=$sample_size \
        --exclude=node07 \
	--partition=${partition} \
        --cpus-per-task=$thread \
	$pip_sh --read1 $read1 \
	--read2 $read2 \
	$temp

#	bash $pip_sh --read1 $read1 --read2 $read2 $temp

done
}

temp=`echo $temp|sed "s/'//g"`
echo $temp


#sep_main

#stat_head

pip_run


set +x
