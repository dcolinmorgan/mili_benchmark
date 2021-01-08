#!/bin/bash
#$ -cwd
## bash data/MotifPipeline/ENCODE/sthlm_motif_pipeline.sh -b50 -c'A549 K562 GM12878 SKNSH HepG2 HeLa'
usage=""$sthlm_motif_pipeline" [-h] [-bco] -- function calling bedtools2 to calculate the intersect for ENCODE methylation benchmark with arbitrary buffer length

where:
	-h  show this help text
  -b  buffer length (integer)
  -c  cell lines (string array 'A549 K562 ...')
  -o  output directory path"

while getopts ":h:b:c:o:" opt; do
  case $opt in
  	h) echo "$usage"
       exit
       ;;
    b) buffer=$OPTARG
    	;;
    c) cells="$OPTARG"
      ;;
    o) outdir={OPTARG:-"data/MotifPipeline/sthlm_motif"$buffer"QC"}
      ;;
    \?) echo "Invalid option -$OPTARG" >&3
    	;;
  esac
done
printf "User has opted to buffer the CpG site with +/- "$buffer" bp \n"


# buffer=0
out_fmt='.txt'


tfdb='data/MotifPipeline/ENCODE/Homo_sapiens_motifinfo.txt'
motifdir='../rekrg/MotifScans/MotifScans/hg38_bed/'
motifs=$(ls $motifdir*)
outdir="data/MotifPipeline/sthlm_motif_"$buffer"_QC"
mkdir -p benchmark_tmp
mkdir -p $outdir

for cell in $cells
do
printf "...running wgbs intersection with array for the "$cell" cell line \n"

eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect  -wa -wb -a data/MotifPipeline/ENCODE/wgbsin/"$cell"both.txt -b data/MotifPipeline/ENCODE/methyl_array/crossQC/"$cell"_MeArrayHG38c.txt" |cut -f1,2,3,4,5,9,10 > benchmark_tmp/overlap_tmpA_$cell"_"$buffer$out_fmt
if [ $buffer != 0 ];then
  printf "...adding +/-"$buffer" bp buffer \n"
	eval "~/../rekrg/Tools/bedtools2/bin/bedtools slop  -i benchmark_tmp/overlap_tmpA_"$cell"_"$buffer".txt -g ~/../rekrg/Tools/bedtools2/genomes/human.hg38.genome -r "$buffer" -l "$buffer"" > benchmark_tmp/overlap_tmpB_$cell"_"$buffer$out_fmt
else
	cp benchmark_tmp/overlap_tmpA_$cell"_"$buffer$out_fmt benchmark_tmp/overlap_tmpB_$cell"_"$buffer$out_fmt
fi
printf "...intersecting methylation information with ChIP for the "$cell" cell line \n"

eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a benchmark_tmp/overlap_tmpB_"$cell"_"$buffer$out_fmt" -b data/MotifPipeline/remap/"$cell"_spRE2020.txt " |cut -f1,2,3,4,5,6,7,8,9,10,11,12 > benchmark_tmp/overlap_tmpC_$cell"_"$buffer$out_fmt

done
printf "intersected all cell lines successfully! \n \n"

for TF in $motifs
do
tf=$(eval "echo "$TF" | cut -d / -f6| cut -d _ -f1")

gene=$(awk -v pat=$tf '$1 ~ pat' $tfdb | cut -f 2)

for cell in $cells
do
  printf "intersecting motif with "$gene" methyl-predictions in the "$cell" cell line \n"
  eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a $TF -b ~/benchmark_tmp/overlap_tmpC_"$cell"_"$buffer$out_fmt""|cut -f1,2,3,5,8,9,10,11,12,13,14,15,16,17,18 > $outdir/$cell"_"$gene
done
done

# rm -rf benchmark_tmp

source /proj/relibs/relib00/conda/bin/activate
source activate mypy3 ## install netZooPy into this pyenv
chmod +x netZooPy/netZooPy/milipeed/benchmark/run_predScore.py 
python netZooPy/netZooPy/milipeed/benchmark/run_predScore.py -i $outdir -o $outdir/test/

# find sthlm_motif_250_QC/ -maxdepth 1 -type f -exec rm -rf {} \;

