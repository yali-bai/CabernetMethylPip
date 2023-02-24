#!/usr/bin/bash
#SBATCH -p compute
#SBATCH -D ./
##SBATCH --mem-per-cpu=2000M
#SBATCH --exclude=node07
##SBATCH --exclude=node01
##SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH -o %A_%a.out


#~~~~~~~~~~~~~~~~~~~~ REFERENCE & SOFTWARE ~~~~~~~~~~~~~~~
python=/share/home/baiyl/basic/tools/anaconda2/bin/python
two_meta_py=/share/home/baiyl/merlot/TAPS-merlot/workflow/RUN_src/src_1.7/Cell.sep.2_meta.py
mark_meta_py=/share/home/baiyl/merlot/TAPS-merlot/workflow/RUN_src/src_1.7/Cell.mark_meta.py
mark_meta_2_py=/share/home/baiyl/merlot/TAPS-merlot/workflow/RUN_src/src_1.8/sep_scripts/Cell.mark_meta.3.py

meta_barcode_py=/share/home/baiyl/merlot/TAPS-merlot/workflow/RUN_src/src_1.7/Cell.sep.meta-barcode.py
index_barcode_py=/share/home/baiyl/merlot/TAPS-merlot/workflow/RUN_src/src_1.7/sep_scripts/Tn5.barcode.find.2.py
index_barcode_py=/share/home/baiyl/merlot/TAPS-merlot/workflow/RUN_src/src_1.7/sep_scripts/Tn5.barcode.find.5.py


meta16_file=/share/home/baiyl/merlot/TAPS-merlot/workflow/RUN_src/src_1.6/data/meta16.tags.txt
meta20_file=/share/home/baiyl/merlot/TAPS-merlot/workflow/RUN_src/src_1.6/data/meta20.tags.txt

index_tag_file=/share/home/baiyl/merlot/TAPS-merlot/workflow/RUN_src/src_1.6/index_tag.txt
barcode_2=/share/home/baiyl/merlot/TAPS-merlot/workflow/RUN_src/src_1.6/barcode_2.txt
barcode_3=/share/home/baiyl/merlot/TAPS-merlot/workflow/RUN_src/src_1.7/data/barcode_3.txt

barcode_new=/share/home/baiyl/merlot/TAPS-merlot/workflow/RUN_src/src_1.8/data/bcNextera_new_barcode.seq.txt


echo "  "
echo "  ------- This is meta tag mark --------"
if [[ $# -lt 1 ]]; then
    echo "No input file with sequence file list specified."
    exit
fi;
#~~~~~~~~~~~~~~~~~~~~~~~~~ ARGS ~~~~~~~~~~~~~~~~~~~~~~~~~
set -x
outdir=$1
Read1=$2
Read2=$3
prefix=$4
trim_type=$5
mindis=$6


echo $outdir


if [[ ! -s $Read1 ]]
then
	echo $Read1 "Read not exist"
	exit
fi

if [[ ! -s $mindis ]]
then
mindis=2
fi

barcode_f=$barcode_new


echo "Start processing with " $Read1, $Read2
meta_tag_file=$meta16_file

#~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~


function mark_sep_cell {
sep_py=$mark_meta_2_py
set -x
$python $sep_py \
        -R1 $Read1 \
        -R2 $Read2 \
	--prefix $prefix \
	--tag_file $meta_tag_file \
	--outdir $outdir
set +x
}

function barcode_sep {
sep_py=$index_barcode_py
$python $sep_py \
        -R1 $Read1 \
        -R2 $Read2 \
        --prefix $prefix \
	--barcode_seq_file $barcode_f \
	--seprate \
	--min_dis $mindis \
        --outdir $outdir

}

if [[ $trim_type == "meta_mark" ]];then
	mark_sep_cell
else
	barcode_sep
fi

set +x
