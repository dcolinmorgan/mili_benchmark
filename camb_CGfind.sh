#!/bin/bash
#$ -cwd
## bash data/MotifPipeline/ENCODE/camb_CGfind.sh

tfdb='data/MotifPipeline/ENCODE/Homo_sapiens_motifinfo.txt'
motifdir='../rekrg/MotifScans/MotifScans/hg38/'
motifs=$(ls $motifdir*)
subdir='data/MotifPipeline/sthlm_motif_0_QCbeta/'
subs=$(ls $subdir*)

for TF in $motifs
do
# tf=$(eval "echo "$TF" | cut -d / -f6| cut -d _ -f1")
# gene=$(awk -v pat=$tf '$1 ~ pat' $tfdb | cut -f 2)
tf=$(eval "echo "$TF" | cut -d / -f6| cut -d _ -f1")
gene=$(awk -v pat=$tf '$1 ~ pat' $tfdb | cut -f 2)

sub=$(eval ls data/MotifPipeline/sthlm_motif_0_QCbeta/ | grep $gene)
cd data/MotifPipeline/sthlm_motif_0_QCbeta/
cat $sub > ~/tmp2.txt
cd ~

# cut -f10 $TF | grep -i 'cg' > tmp.txt
cut -f3,4,5,10 $TF | grep -i 'chr' > tmp.txt

eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -a tmp.txt -b tmp2.txt " >  tmp3.txt

cat tmp3.txt | grep -i 'cg' > tmp4.txt

cut -f4 tmp4.txt > tmp.txt

CG=$(eval wc -l tmp4.txt | cut -f1 --delimiter=' ')
full=$(eval wc -l tmp3.txt | cut -f1 --delimiter=' ')

echo $CG $full $tf $gene  >> CGcont.txt
echo $CG $full $tf $gene #>> tmp3.txt

rm -rf tmp.txt tmp2.txt tmp3.txt
done

M0952 SIX2