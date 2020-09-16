#!/bin/bash
#$ -cwd
## bash data/MotifPipeline/ENCODE/camb_motif_pipeline_delta.sh -b5 -c'A549 K562 GM12878 SKNSH HepG2 HeLa'
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
    o) outdir={OPTARG:-"data/MotifPipeline/camb_motif"$buffer"QC"}
      ;;
    \?) echo "Invalid option -$OPTARG" >&3
      ;;
  esac
done
printf "User has opted to buffer the motif region with +/- "$buffer" bp \n"


# buffer=0
out_fmt='.txt'


tfdb='data/MotifPipeline/ENCODE/Homo_sapiens_motifinfo.txt'
motifdir='../rekrg/MotifScans/MotifScans/hg38/'
motifs=$(ls $motifdir*)
outdir="../../../pc/redmo/data/MotifPipeline/camb_motif_"$buffer"_QCdelta"

mkdir -p benchmark_tmp
# mkdir -p $outdir
rm -rf $outdir
mkdir -p $outdir/red

for cell in $cells
do
  printf "...running wgbs intersection with array for the "$cell" cell line \n"

  eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect  -wa -wb -a data/MotifPipeline/ENCODE/wgbsin/"$cell"both.txt -b data/MotifPipeline/ENCODE/methyl_array/crossQC/"$cell"_MeArrayHG38c.txt" |cut -f1,2,3,4,5,9,10 > benchmark_tmp/overlap_tmpA_$cell"_"$buffer$out_fmt

done
printf "intersected all cell line methyl info successfully! \n \n"

for TF in $motifs
do
  tf=$(eval "echo "$TF" | cut -d / -f6| cut -d _ -f1")

  gene=$(awk -v pat=$tf '$1 ~ pat' $tfdb | cut -f 2)
  if [ $buffer != 0 ];then
    printf "...adding +/-"$buffer" bp buffer to motif region \n"
    cut -f3,4,5,7,10 $TF | tail -n +2 > benchmark_tmp/overlap_tmpD_$cell"_"$buffer$out_fmt
    eval "~/../rekrg/Tools/bedtools2/bin/bedtools slop  -i benchmark_tmp/overlap_tmpD_$cell"_"$buffer$out_fmt -g ~/../rekrg/Tools/bedtools2/genomes/human.hg38.genome -r "$buffer" -l "$buffer"" > benchmark_tmp/overlap_tmpB_$cell"_"$buffer$out_fmt
    
    for cell in $cells
    do
      printf "intersecting buffered motif with "$gene" methyl-predictions in the "$cell" cell line \n"
      

      eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a benchmark_tmp/overlap_tmpB_$cell"_"$buffer$out_fmt -b ~/benchmark_tmp/overlap_tmpA_"$cell"_"$buffer$out_fmt"" > benchmark_tmp/overlap_tmpC_$cell"_"$buffer$out_fmt
      
      printf "...intersecting methylation information with ChIP for the "$cell" cell line \n"
      eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a benchmark_tmp/overlap_tmpC_"$cell"_"$buffer$out_fmt" -b data/MotifPipeline/remap/"$cell"_spRE2020.txt " >  $outdir/$cell"_"$gene
      
      printf "...reducing buffered ("$buffer"bp) motif regions to only those in common w/ 0 \n"
      eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -a $outdir/$cell"_"$gene -b data/MotifPipeline/sthlm_motif_0_QCbeta/$cell"_"$gene " >  $outdir/red/$cell"_"$gene
      rm $outdir/$cell"_"$gene benchmark_tmp/overlap_tmpD_$cell"_"$buffer$out_fmt # benchmark_tmp/overlap_tmB_$cell"_"$buffer$out_fmt benchmark_tmp/overlap_tmpC_$cell"_"$buffer$out_fmt
    done
  else
    for cell in $cells
    do
      eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a $TF -b ~/benchmark_tmp/overlap_tmpB_"$cell"_"$buffer$out_fmt"" > benchmark_tmp/overlap_tmpC_$cell"_"$buffer$out_fmt
      printf "...intersecting methylation information with ChIP for the "$cell" cell line \n"

      eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a benchmark_tmp/overlap_tmpC_"$cell"_"$buffer$out_fmt" -b data/MotifPipeline/remap/"$cell"_spRE2020.txt "  >  $outdir/$cell"_"$gene
    done
  fi

done
done
find $outdir/red/ -size +5G -delete

source /proj/relibs/relib00/conda/bin/activate
source activate mypy3 ## install netZooPy into this pyenv
chmod +x netZooPy/netZooPy/milipeed/benchmark/run_predScore.py 
if [ $buffer != 0 ];then
  python netZooPy/netZooPy/milipeed/benchmark/run_predScore.py -i $outdir/red -o $outdir/red/test/
  # find "data/MotifPipeline/camb_motif_"$buffer"_QCdelta/" -maxdepth 1 -type f -exec rm -rf {} \;
else
  python netZooPy/netZooPy/milipeed/benchmark/run_predScore.py -i $outdir -o $outdir/red/test/
fi

# find "data/MotifPipeline/camb_motif_20_QCdelta/" -maxdepth 2 -type f -exec rm -rf {} \;

