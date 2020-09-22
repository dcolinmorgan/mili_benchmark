#!/bin/bash
#$ -cwd
## bash mili_benchmark/camb_motif_pipeline_beta.sh -b5 -c'A549 K562 GM12878 SKNSH HepG2 HeLa' -o'../../../pc/redmo/data/MotifPipeline/'
usage=""$camb_motif_pipeline" [-h] [-bco] -- function calling bedtools2 to calculate the intersect for ENCODE methylation benchmark with arbitrary buffer length

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
    o) bench={OPTARG:-"data/MotifPipeline/camb_motif"$buffer"QC"}
      ;;
    \?) echo "Invalid option -$OPTARG" >&3
      ;;
  esac
done
if [ $buffer != 0 ];then
  printf "User has opted to buffer the motif region with +/- "$buffer" bp \n"
else
    printf "intersection running for CG content analysis without added buffer around motif \n"
fi

# buffer=0
out_fmt='.txt'

# me=`basename "$0"` | cut -d . -f1
me=$(eval "basename "$0" | cut -d . -f1")

tfdb='data/MotifPipeline/ENCODE/Homo_sapiens_motifinfo.txt'
motifdir='../rekrg/MotifScans/MotifScans/hg38/'
motifs=$(ls $motifdir*)
bench="../../../pc/redmo/data/expand_base0_MotifPipeline/"
bench=$bench$me$buffer

# rm -rf bench
# mkdir bench
rm -rf $bench
mkdir -p $bench/test
# bench="benchmark_"
# bench=$bench$me$buffer

# if [ ! -d "$bench" ]; then
# mkdir benchmark_$me$buffer
for cell in $cells
do
  # printf "...running wgbs intersection with array for the "$cell" cell line \n"

  # eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect  -wa -wb -a data/MotifPipeline/ENCODE/wgbsin/"$cell"both.txt -b data/MotifPipeline/ENCODE/methyl_array/crossQC/"$cell"_MeArrayHG38c.txt" |cut -f1,2,3,4,5,9,10 > $bench/overlap_tmpA_$cell"_"$buffer$out_fmt
  # printf "...intersecting methylation information with ChIP for the "$cell" cell line \n"

# done
# fi

# for cell in $cells
# do
  for TF in $motifs
  do
    tf=$(eval "echo "$TF" | cut -d / -f6| cut -d _ -f1")

    gene=$(awk -v pat=$tf '$1 ~ pat' $tfdb | cut -f 2)
    cut -f3,4,5,7,10 $TF | tail -n +2 > $bench/overlap_tmpD_$cell"_"$buffer$out_fmt
    if [ $buffer != 0 ];then
      printf "...adding +/-"$buffer" bp buffer around motif region for "$gene" in the "$cell" cell line \n"
      eval "~/../rekrg/Tools/bedtools2/bin/bedtools slop  -i $bench/overlap_tmpD_$cell"_"$buffer$out_fmt -g ~/../rekrg/Tools/bedtools2/genomes/human.hg38.genome -r "$buffer" -l "$buffer"" > $bench/overlap_tmpC_$cell"_"$buffer$out_fmt
      eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a $bench/overlap_tmpC_"$cell"_"$buffer$out_fmt" -b data/MotifPipeline/remap/"$cell"_spRE2020.txt " >  $bench/overlap_tmpB_$cell"_"$buffer$out_fmt #$bench/$cell"_"$gene
      sed -i.bak $'s/\t\t/\t/' $bench/overlap_tmpB_$cell"_"$buffer$out_fmt
      cut -f1,2,3,4,5,11 $bench/overlap_tmpB_$cell"_"$buffer$out_fmt > $bench/overlap_tmpBB_$cell"_"$buffer$out_fmt
      eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -a $bench/overlap_tmpBB_$cell"_"$buffer$out_fmt -b $bench/overlap_tmpBB_$cell"_"0$out_fmt " >  $outdir/test/$cell"_"$gene

    else
      printf "...no buffered motif for "$gene" in the "$cell" cell line \n"
      cat $bench/overlap_tmpD_$cell"_"$buffer$out_fmt > $bench/overlap_tmpC_$cell"_"$buffer$out_fmt
      eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a $bench/overlap_tmpC_"$cell"_"$buffer$out_fmt" -b data/MotifPipeline/remap/"$cell"_spRE2020.txt " >  $bench/overlap_tmpB_$cell"_"$buffer$out_fmt #$bench/$cell"_"$gene
      sed -i.bak $'s/\t\t/\t/' $bench/overlap_tmpB_$cell"_"$buffer$out_fmt
      cut -f1,2,3,4,5,11 $bench/overlap_tmpB_$cell"_"$buffer$out_fmt > $bench/overlap_tmpBB_$cell"_"$buffer$out_fmt
    fi
    printf "intersecting buffered motif with "$gene" methyl-predictions in the "$cell" cell line \n"

    eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a $bench/overlap_tmpC_$cell"_"$buffer$out_fmt -b $bench/overlap_tmpBB_"$cell"_"$buffer$out_fmt"" > $bench/$cell"_"$gene #bench/overlap_tmpE_$cell"_"$buffer$out_fmt
    
    # printf "...reducing buffered ("$buffer"bp) motif regions to only those in common w/ 0 \n"
    # eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -a bench/overlap_tmpE_$cell"_"$buffer$out_fmt -b data/MotifPipeline/sthlm_motif_0_QCbeta/$cell"_"$gene " >  $bench/red/$cell"_"$gene
    # rm $bench/$cell"_"$gene bench/overlap_tmpD_$cell"_"$buffer$out_fmt  bench/overlap_tmE_$cell"_"$buffer$out_fmt bench/overlap_tmpC_$cell"_"$buffer$out_fmt

done
done
printf "intersected all cell line methyl-methyl-chip info successfully! \n \n"
s
# done
# find $bench/test/ -size +5G -delete

source /proj/relibs/relib00/conda/bin/activate
source activate mypy3 ## install netZooPy into this pyenv
chmod +x netZooPy/netZooPy/milipeed/benchmark/run_predScore.py 
# if [ $buffer != 0 ];then
#   python netZooPy/netZooPy/milipeed/benchmark/run_predScore.py -i $bench/red -o $bench/red/test/
  # find "data/MotifPipeline/camb_motif_"$buffer"_QCdelta/" -maxdepth 1 -type f -exec rm -rf {} \;
# else
python mili_benchmark/run_cambPredScore.py -i $bench -o $bench/test/
# fi
# rm -rf $bench
# find "data/MotifPipeline/camb_motif_20_QCdelta/" -maxdepth 2 -type f -exec rm -rf {} \;

