#!/bin/bash
#$ -cwd
## bash data/MotifPipeline/ENCODE/sthlm_CGfind.sh

tfdb='data/MotifPipeline/ENCODE/Homo_sapiens_motifinfo.txt'
motifdir='../rekrg/MotifScans/MotifScans/hg38/'
motifs=$(ls $motifdir*)

for TF in $motifs
do
# tf=$(eval "echo "$TF" | cut -d / -f6| cut -d _ -f1")
# gene=$(awk -v pat=$tf '$1 ~ pat' $tfdb | cut -f 2)
tf=$(eval "echo "$TF" | cut -d / -f6| cut -d _ -f1")
gene=$(awk -v pat=$tf '$1 ~ pat' $tfdb | cut -f 2)

cut -f10 $TF | grep -i 'cg' > tmp.txt
CG=$(eval wc -l tmp.txt | cut -f1 --delimiter=' ')
full=$(eval wc -l $TF | cut -f1 --delimiter=' ')

echo $CG $full $tf $gene  >> tmp1.txt
echo $CG $full $tf $gene #>> tmp3.txt

done