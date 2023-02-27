#!/usr/bin/bash
python=/share/home/baiyl/basic/tools/anaconda2/bin/python
RUN_src=/share/home/baiyl/merlot/TAPS-merlot/workflow/RUN_src/src_1.9
sep_py=$RUN_src/scripts/Cell.mark_meta.6.py 


temp=`getopt -o o:s:p:m:t: --long outdir:,sep_type:,prefix:,meta_type:,thread:,tag_file:,read1:,read2: -- "$@"`
if [ $? != 0 ] ; then echo "terminating..." >&2 ; exit 1 ; fi
eval set -- "$temp" 

separate=meta_mark
while true ; do
        case "$1" in
                -o|--outdir) echo "Output dir is $2" ;outdir=$2; shift 2;;
                -s|--sep_type) separate=$2; echo "Separate the reads with barcodes. Default: meta_mark" ; shift 2;;
		-p|--prefix) prefix=$2; echo "Reads prefix is $prefix" ; shift 2;;
		-m|--meta_type) meta_type=$2; echo "meta type is $meta_type"; shift 2 ;;
		-t|--thread) thread=$2; echo "multi-thread is $thread"; shift 2 ;;
		--tag_file) tag_file=$2; echo "Tag file is $tag_file"; shift 2 ;;
		--read1) read1=$2;shift 2 ;;
		--read2) read2=$2;shift 2 ;;

                --) shift ; break ;;
                *) echo "internal error!" ; exit 1 ;;
        esac
done

min_dis=4

set -x

if [[ ! -n $thread ]];then thread=1;fi

if [[ ! -n $outdir ]];then outdir=".";fi

if [[ ! -n $meta_type ]] && [[ ! -n $tag_file ]];then
	echo "there is no tag file input or chosed."
	exit
fi

declare -A barcode_f
barcode_f["meta16"]=$RUN_src/data/meta16.tags.txt
barcode_f["meta20"]=$RUN_src/data/meta20.tags.txt
barcode_f["Cmeta20"]=$RUN_src/data/Cmeta20.tags.txt 
barcode_f["metaGAT"]=$RUN_src/data/metaGAT.tags.txt
barcode_f["metaGATi"]=$RUN_src/data/metaGATi.tags.txt

if [[ ! -n $tag_file ]];
then
	tag_file=${barcode_f[$meta_type]}
fi


if [[ $sep_type == "meta_mark" ]];
then
$python $sep_py \
        -R1 $read1 \
        -R2 $read2 \
        --prefix $prefix \
        --tag_file $tag_file \
        --threads $thread \
        --bypart --split_part $thread \
        --min_dis $min_dis \
        --outdir $outdir/sep_cell
else
$python $sep_py \
        -R1 $read1 \
        -R2 $read2 \
        --prefix $prefix \
        --tag_file $tag_file \
        --threads $thread \
        --bypart --split_part $thread \
        --min_dis $min_dis \
        --outdir $outdir/sep_cell

fi



set +x
