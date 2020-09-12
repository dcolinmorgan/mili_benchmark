#!/bin/bash
#$ -cwd
## bash data/MotifPipeline/ENCODE/camb_motif_pipeline.sh -c'A549 K562 GM12878 SKNSH HepG2 HeLa' -g'HOXA7 RFX6 ZBTB49'
usage=""$sthlm_motif_pipeline" [-h] [-cgo] -- function calling bedtools2 to calculate the intersect for ENCODE methylation benchmark with arbitrary buffer length

where:
  -h  show this help text
  -c  cell lines (string array 'A549 K562 ...')
  -g  genes of interest (string 'HOXA7 RFX6 ...')
  -o  output directory path"

while getopts ":h:c:g:o:" opt; do
  case $opt in
    h) echo "$usage"
       exit
       ;;
    c) cells="$OPTARG"
      ;;
    g) genes="$OPTARG":-'all'
      ;;
    o) outdir={OPTARG:-"../../../d/tmp/redmo/camb_motif"}
      ;;
    \?) echo "Invalid option -$OPTARG" >&4
      ;;
  esac
done

out_fmt='.txt'


tfdb='data/MotifPipeline/ENCODE/Homo_sapiens_motifinfo.txt'
motifdir='../rekrg/MotifScans/MotifScans/hg38/'
motifs=$(ls $motifdir*)
outdir="../../../d/tmp/redmo/camb_motif"

rm -rf $outdir
mkdir -p $outdir
mkdir -p $outdir/benchmark_tmp
mkdir -p $outdir/nonCG
mkdir -p $outdir/CG


if [[ $* == *-g* ]];then
  # tf=$(eval "echo "$TF" | cut -d / -f6| cut -d _ -f1")
  for gene in $genes
  do
    TF=$(awk -v pat=$gene '$2 ~ pat' $tfdb | cut -f 1)
    tf=$(eval "echo "$TF" | cut -d / -f6| cut -d _ -f1")
    motif_fmt='_1.02.txt'
    for cell in $cells
    do
      cut -f3,4,5,7,9 $motifdir/$TF$motif_fmt | tail -n +2 > $outdir/benchmark_tmp/tmp1.txt

      eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a $outdir/benchmark_tmp/tmp1.txt -b data/MotifPipeline/remap/"$cell"_spRE2020.txt" > $outdir/benchmark_tmp/tmp2.txt
      cat $outdir/benchmark_tmp/tmp2.txt | grep -i 'cg' > $outdir/benchmark_tmp/tmp3.txt
      cut -f4,8 $outdir/benchmark_tmp/tmp3.txt > $outdir/CG/$cell"_"$gene

      cat $outdir/benchmark_tmp/tmp2.txt | grep -i -v 'cg' > $outdir/benchmark_tmp/tmp4.txt
      cut -f4,8 $outdir/benchmark_tmp/tmp4.txt > $outdir/nonCG/$cell"_"$gene

      CG=$(eval wc -l $outdir/CG/$cell"_"$gene | cut -f1 --delimiter=' ')
      nonCG=$(eval wc -l $outdir/nonCG/$cell"_"$gene | cut -f1 --delimiter=' ')
      
      echo $CG $nonCG $tf $gene $cell

      rm $outdir/benchmark_tmp/tmp3.txt $outdir/benchmark_tmp/tmp2.txt $outdir/benchmark_tmp/tmp4.txt $outdir/benchmark_tmp/tmp1.txt

    done
  done

else
  for TF in $motifs
    do
    tf=$(eval "echo "$TF" | cut -d / -f6| cut -d _ -f1")

    gene=$(awk -v pat=$tf '$1 ~ pat' $tfdb | cut -f 2)

    for cell in $cells
    do
      cut -f3,4,5,7,9 $TF | tail -n +2 > $outdir/benchmark_tmp/tmp1.txt

      eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a $outdir/benchmark_tmp/tmp1.txt -b data/MotifPipeline/remap/"$cell"_spRE2020.txt" > $outdir/benchmark_tmp/tmp2.txt
      cat $outdir/benchmark_tmp/tmp2.txt | grep -i 'cg' > $outdir/benchmark_tmp/tmp3.txt
      cut -f4,8 $outdir/benchmark_tmp/tmp3.txt > $outdir/CG/$cell"_"$gene

      cat $outdir/benchmark_tmp/tmp2.txt | grep -i -v 'cg' > $outdir/benchmark_tmp/tmp4.txt
      cut -f4,8 $outdir/benchmark_tmp/tmp4.txt > $outdir/nonCG/$cell"_"$gene

      CG=$(eval wc -l $outdir/CG/$cell"_"$gene | cut -f1 --delimiter=' ')
      nonCG=$(eval wc -l $outdir/nonCG/$cell"_"$gene | cut -f1 --delimiter=' ')
      
      echo $CG $nonCg $tf $gene $cell

      rm $outdir/benchmark_tmp/tmp3.txt $outdir/benchmark_tmp/tmp2.txt $outdir/benchmark_tmp/tmp4.txt $outdir/benchmark_tmp/tmp1.txt

    done
  done
fi

