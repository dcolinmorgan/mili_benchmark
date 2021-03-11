#!/bin/bash
#$ -cwd
## bash data/MotifPipeline/ENCODE/camb_motif_pipeline.sh -c'A549' -g'HOXA7 RFX6 ZBTB4' -m'array'
usage=""$sthlm_motif_pipeline" [-h][-cgmo]  -- function calling bedtools2 to calculate the intersect for ENCODE methylation benchmark with arbitrary buffer length

where:
  -h  show this help text
  -c  cell lines (string array 'A549 K562 ...')
  -g  genes of interest (string 'HOXA7 RFX6 ...')
  -m  methylation type (wgbs or array)
  -o  output directory path"
  
while getopts ":h:c:g:m:o:" opt; do
  case $opt in
    h) echo "$usage"
       exit
       ;;
    c) cells=$OPTARG #:-"A549 K562 GM12878 SKNSH HepG2 HeLa"
      ;;
    g) genes=$OPTARG
      ;;
    m) methyl="$OPTARG"
      ;;
    o) outdir={OPTARG:-"data/MotifPipeline/sthlm_motif"$buffer"QC"}
      ;;
    \?) echo "Invalid option -$OPTARG" >&3
      ;;
  esac
done

out_fmt='.txt'
tfdb='data/MotifPipeline/ENCODE/Homo_sapiens_motifinfo.txt'
motifdir='../rekrg/MotifScans/MotifScans/hg38/'
motifs=$(ls $motifdir*)
outdir="../../../d/tmp/redmo/camb_motif_"
outdir=$outdir$methyl

rm -rf $outdir
mkdir -p $outdir
mkdir -p $outdir/benchmark_tmp
mkdir -p $outdir/nonCG/
mkdir -p $outdir/CG/
mkdir -p $outdir/perf/

if [[ $* == *-g* ]];then
  # tf=$(eval "echo "$TF" | cut -d / -f6| cut -d _ -f1")
  for gene in $genes
  do
    TF=$(awk -v pat=$gene '$2 ~ pat' $tfdb | cut -f 1)
    tf=$(eval "echo "$TF" | cut -d / -f6| cut -d _ -f1")
    motif_fmt='_1.02.txt'
    for cell in $cells
    do
      cut -f3,4,5,7,10 $motifdir/$TF$motif_fmt | tail -n +2 > $outdir/benchmark_tmp/tmp1.txt

      eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a $outdir/benchmark_tmp/tmp1.txt -b data/MotifPipeline/remap/"$cell"_spRE2020.txt" > $outdir/benchmark_tmp/tmp22.txt
      sed -i.bak $'s/\t\t/\t/' $outdir/benchmark_tmp/tmp22.txt
      cut -f1,2,3,4,5,9 $outdir/benchmark_tmp/tmp22.txt > $outdir/benchmark_tmp/tmp222.txt
      if [[ $methyl == *"array"* ]];then
        printf "running array interX"
        eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a $outdir/benchmark_tmp/tmp222.txt -b data/MotifPipeline/ENCODE/methyl_array/crossQC/"$cell"_MeArrayHG38c.txt" > $outdir/benchmark_tmp/tmp2.txt
        # eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a $outdir/nonCG/$cell"_"$gene -b data/MotifPipeline/ENCODE/methyl_array/crossQC/"$cell"_MeArrayHG38c.txt" > $outdir/benchmark_tmp/tmp2.txt
        cat $outdir/benchmark_tmp/tmp2.txt | grep -i 'cg' > $outdir/benchmark_tmp/tmp3.txt
        cut -f1,2,3,4,5,6,10,11 $outdir/benchmark_tmp/tmp3.txt > $outdir/CG/$cell"_"$gene

        cat $outdir/benchmark_tmp/tmp2.txt | grep -i -v 'cg' > $outdir/benchmark_tmp/tmp4.txt
        cut -f1,2,3,4,5,6,10,11 $outdir/benchmark_tmp/tmp4.txt > $outdir/nonCG/$cell"_"$gene

        CG=$(eval wc -l $outdir/CG/$cell"_"$gene | cut -f1 --delimiter=' ')
        nonCG=$(eval wc -l $outdir/nonCG/$cell"_"$gene | cut -f1 --delimiter=' ')
        total=$(eval wc -l $outdir/benchmark_tmp/tmp2.txt| cut -f1 --delimiter=' ')
        echo $CG $nonCG $total $tf $gene $cell
      elif [[ $methyl == *"wgbs"* ]];then
        printf "running wgbs interX"
        eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a $outdir/benchmark_tmp/tmp222.txt -b data/MotifPipeline/ENCODE/wgbsin/"$cell"both.txt" > $outdir/benchmark_tmp/tmp2.txt
        # eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a $outdir/nonCG/$cell"_"$gene -b data/MotifPipeline/ENCODE/wgbsin/"$cell"both.txt" > $outdir/benchmark_tmp/tmp2.txt
        cat $outdir/benchmark_tmp/tmp2.txt | grep -i 'cg' > $outdir/benchmark_tmp/tmp3.txt
        cut -f1,2,3,4,5,6,10,11 $outdir/benchmark_tmp/tmp3.txt > $outdir/CG/$cell"_"$gene

        cat $outdir/benchmark_tmp/tmp2.txt | grep -i -v 'cg' > $outdir/benchmark_tmp/tmp4.txt
        cut -f1,2,3,4,5,6,10,11 $outdir/benchmark_tmp/tmp4.txt > $outdir/nonCG/$cell"_"$gene
        CG=$(eval wc -l $outdir/CG/$cell"_"$gene | cut -f1 --delimiter=' ')
        nonCG=$(eval wc -l $outdir/nonCG/$cell"_"$gene | cut -f1 --delimiter=' ')
        total=$(eval wc -l $outdir/benchmark_tmp/tmp2.txt| cut -f1 --delimiter=' ')
        echo $CG $nonCG $total $tf $gene $cell
      else
        cat $outdir/benchmark_tmp/tmp2.txt | grep -i 'cg' > $outdir/benchmark_tmp/tmp3.txt
        cut -f4,9 $outdir/benchmark_tmp/tmp3.txt > $outdir/CG/$cell"_"$gene

        cat $outdir/benchmark_tmp/tmp2.txt | grep -i -v 'cg' > $outdir/benchmark_tmp/tmp4.txt
        cut -f4,9 $outdir/benchmark_tmp/tmp4.txt > $outdir/nonCG/$cell"_"$gene
        CG=$(eval wc -l $outdir/CG/$cell"_"$gene | cut -f1 --delimiter=' ')
        nonCG=$(eval wc -l $outdir/nonCG/$cell"_"$gene | cut -f1 --delimiter=' ')
        total=$(eval wc -l $outdir/benchmark_tmp/tmp2.txt| cut -f1 --delimiter=' ')
        
        echo $CG $nonCG $total $tf $gene $cell
        echo $CG $nonCG $total $tf $gene $cell  >> $outdir/benchmark_tmp/CGcont.txt
      fi

      # rm $outdir/benchmark_tmp/tmp3.txt $outdir/benchmark_tmp/tmp2.txt $outdir/benchmark_tmp/tmp4.txt $outdir/benchmark_tmp/tmp1.txt

    done
  done

else
  for TF in $motifs
    do
    tf=$(eval "echo "$TF" | cut -d / -f6| cut -d _ -f1")

    gene=$(awk -v pat=$tf '$1 ~ pat' $tfdb | cut -f 2)

    for cell in $cells
    do
      cut -f3,4,5,7,10 $motifdir/$TF$motif_fmt | tail -n +2 > $outdir/benchmark_tmp/tmp1.txt

      eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a $outdir/benchmark_tmp/tmp1.txt -b data/MotifPipeline/remap/"$cell"_spRE2020.txt" > $outdir/benchmark_tmp/tmp22.txt
      sed -i.bak $'s/\t\t/\t/' $outdir/benchmark_tmp/tmp22.txt
      cut -f1,2,3,4,5,9 $outdir/benchmark_tmp/tmp22.txt > $outdir/benchmark_tmp/tmp222.txt
      if [[ $methyl == *"array"* ]];then
        printf "running array interX"
        eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a $outdir/benchmark_tmp/tmp222.txt -b data/MotifPipeline/ENCODE/methyl_array/crossQC/"$cell"_MeArrayHG38c.txt" > $outdir/benchmark_tmp/tmp2.txt
        # eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a $outdir/nonCG/$cell"_"$gene -b data/MotifPipeline/ENCODE/methyl_array/crossQC/"$cell"_MeArrayHG38c.txt" > $outdir/benchmark_tmp/tmp2.txt
        cat $outdir/benchmark_tmp/tmp2.txt | grep -i 'cg' > $outdir/benchmark_tmp/tmp3.txt
        cut -f1,2,3,4,5,6,10,11 $outdir/benchmark_tmp/tmp3.txt > $outdir/CG/$cell"_"$gene

        cat $outdir/benchmark_tmp/tmp2.txt | grep -i -v 'cg' > $outdir/benchmark_tmp/tmp4.txt
        cut -f1,2,3,4,5,6,10,11 $outdir/benchmark_tmp/tmp4.txt > $outdir/nonCG/$cell"_"$gene

        CG=$(eval wc -l $outdir/CG/$cell"_"$gene | cut -f1 --delimiter=' ')
        nonCG=$(eval wc -l $outdir/nonCG/$cell"_"$gene | cut -f1 --delimiter=' ')
        total=$(eval wc -l $outdir/benchmark_tmp/tmp2.txt| cut -f1 --delimiter=' ')
        echo $CG $nonCG $total $tf $gene $cell
      elif [[ $methyl == *"wgbs"* ]];then
        printf "running wgbs interX"
        eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a $outdir/benchmark_tmp/tmp222.txt -b data/MotifPipeline/ENCODE/wgbsin/"$cell"both.txt" > $outdir/benchmark_tmp/tmp2.txt
        # eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a $outdir/nonCG/$cell"_"$gene -b data/MotifPipeline/ENCODE/wgbsin/"$cell"both.txt" > $outdir/benchmark_tmp/tmp2.txt
        cat $outdir/benchmark_tmp/tmp2.txt | grep -i 'cg' > $outdir/benchmark_tmp/tmp3.txt
        cut -f1,2,3,4,5,6,10,11 $outdir/benchmark_tmp/tmp3.txt > $outdir/CG/$cell"_"$gene

        cat $outdir/benchmark_tmp/tmp2.txt | grep -i -v 'cg' > $outdir/benchmark_tmp/tmp4.txt
        cut -f1,2,3,4,5,6,10,11 $outdir/benchmark_tmp/tmp4.txt > $outdir/nonCG/$cell"_"$gene
        CG=$(eval wc -l $outdir/CG/$cell"_"$gene | cut -f1 --delimiter=' ')
        nonCG=$(eval wc -l $outdir/nonCG/$cell"_"$gene | cut -f1 --delimiter=' ')
        total=$(eval wc -l $outdir/benchmark_tmp/tmp2.txt| cut -f1 --delimiter=' ')
        echo $CG $nonCG $total $tf $gene $cell
      else
        cat $outdir/benchmark_tmp/tmp2.txt | grep -i 'cg' > $outdir/benchmark_tmp/tmp3.txt
        cut -f4,9 $outdir/benchmark_tmp/tmp3.txt > $outdir/CG/$cell"_"$gene

        cat $outdir/benchmark_tmp/tmp2.txt | grep -i -v 'cg' > $outdir/benchmark_tmp/tmp4.txt
        cut -f4,9 $outdir/benchmark_tmp/tmp4.txt > $outdir/nonCG/$cell"_"$gene
        CG=$(eval wc -l $outdir/CG/$cell"_"$gene | cut -f1 --delimiter=' ')
        nonCG=$(eval wc -l $outdir/nonCG/$cell"_"$gene | cut -f1 --delimiter=' ')
        total=$(eval wc -l $outdir/benchmark_tmp/tmp2.txt| cut -f1 --delimiter=' ')
        
        echo $CG $nonCG $total $tf $gene $cell
        echo $CG $nonCG $total $tf $gene $cell  >> $outdir/benchmark_tmp/CGcont.txt
      fi

      # rm $outdir/benchmark_tmp/tmp3.txt $outdir/benchmark_tmp/tmp2.txt $outdir/benchmark_tmp/tmp4.txt $outdir/benchmark_tmp/tmp1.txt

    done
  done
fi

source /proj/relibs/relib00/conda/bin/activate
source activate mypy3 ## install netZooPy into this pyenv
chmod +x netZooPy/netZooPy/milipeed/benchmark/run_predScore.py 
# python netZooPy/netZooPy/milipeed/benchmark/run_predScore.py -i $outdir -o $outdir/perf/
# mkdir ../../../d/tmp/redmo/camb_motif_ALL/perf/
python netZooPy/netZooPy/milipeed/benchmark/run_predScore.py -i $outdir -o $outdir/perf/
