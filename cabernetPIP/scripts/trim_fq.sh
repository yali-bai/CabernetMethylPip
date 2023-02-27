#!/usr/bin/bash
#SBATCH -p compute
#SBATCH -D ./
#SBATCH --mem-per-cpu=20000M
#SBATCH --exclude=node07
##SBATCH --exclude=node01
##SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH -o %J_%A_%a.out
trim_galore=/share/newdata4/baiyl/basic/tools/TrimGalore-0.4.5/trim_galore

CRC_QC=/share/home/baiyl/merlot/TAPS-merlot/workflow/Analysis/CRC_2018_run/GSE97693_colon_2018_RRBS-CRC_RNA-seq-CRC/CRC_code/Methylation/bin/QC_lite_0.0.7c.pl

if [[ $# -lt 1 ]]; then
    echo "No input file with sequence file list specified."
    exit
fi;
set -x

temp=`getopt -o o:t: --long outdir:,trim_type:,read1:,read2:,read_len:,clip_len:  -- "$@"`
if [ $? != 0 ] ; then echo "terminating..." >&2 ; exit 1 ; fi
eval set -- "$temp" 

while true ; do
        case "$1" in
                -o|--outdir) echo "Output dir is $2" ;outdir=$2; shift 2;;
                -t|--trim_type) trim_type=$2;echo "option b, argument \`$2'" ; shift 2;;
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

                --) shift ; break ;;
                *) echo "internal error!" ; exit 1 ;;
        esac
done


echo "Start trimming: "
echo "Outdir is " $outdir


if [ ! -s $Read1 ];then
        echo $Read1 "Read not exist"
        exit
fi
if [ ! -s $Read2 ];then
	echo $Read2 "Read not exist"
fi

if [[ ! $read_len ]];then
	read_len=20
fi
if [[ ! $clip_len ]];then
	clip_len=9
fi
if [[ ! $outdir ]];then
	outdir="trim"
	mkdir -p $outdir
fi


echo "Start processing with ", $Read1, $Read2
sampleid=`basename $Read1 |sed 's/.R1.*//g'|xargs`

if [[ -f $outdir/${sampleid}.R1_val_1.fq.gz ]] && [[ -f $outdir/${sampleid}.R2_val_2.fq.gz ]];then
        echo "Trim already finished"
#        return
else
        rm -fr $outdir/${sampleid}*
fi


function basic_trim {
read_len=$1
$trim_galore --fastqc --paired \
	--phred33 \
	--length $read_len \
	--retain_unpaired \
	--output_dir $outdir \
	$Read1 $Read2
}

function nextera_basic_trim {
clip_len=$1

$trim_galore --fastqc --paired \
	--phred33 \
	--retain_unpaired \
	--length $read_len \
	--clip_R1 $clip_len \
	--clip_R2 $clip_len \
	--three_prime_clip_R1 $clip_len \
	--three_prime_clip_R2 $clip_len \
	--output_dir $outdir \
	$Read1 $Read2
}

function clip5_trim {
set -x
clip_len=$1
$trim_galore --fastqc --paired \
        --phred33 \
        --length $read_len \
        --retain_unpaired \
        --output_dir $outdir \
        --clip_R1 $clip_len \
        --clip_R2 $clip_len \
        $Read1 $Read2

#$trim_galore \
#        --phred33 \
#	--length $read_len \
#        --retain_unpaired \
#        --output_dir $outdir \
#        --clip_R1 $clip_len \
#	$Read1

#$trim_galore --fastqc --paired \
#$trim_galore \
#        --phred33 \
#        --retain_unpaired \
#        --length $read_len \
#        --clip_R1 $clip_len \
#        --three_prime_clip_R1 $clip_len \
#        --output_dir $outdir \
#        $Read1



set +x
}

if [[ $trim_type =~ "basic_trim" ]];then
	basic_trim $read_len
elif [[ $trim_type =~ "nextera_trim" ]];then
	nextera_basic_trim $clip_len
elif [[ $trim_type =~ "clip5_trim" ]];then
	clip5_trim $clip_len
fi




set +x
