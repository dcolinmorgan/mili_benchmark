#!/bin/bash
#bash mili_benchmark/src/shell/mega_motif.sh

motifdir='/udd/rekrg/MotifScans/MotifScans/hg38'  ## chr start stop pwm per gene name/
motiffiles=$(ls $motifdir/*)
# echo $depth
# rm -rf ~/counts.txt
# echo $motiffiles
for tf in $motiffiles #
do
	awk '{print($3,"\t",$4,"\t",$5,"\t",$7,"\t",$2)}' $tf |tail -n +2>>mega_motif.txt
    # cut -f3,4,5,7,3 $tf
done
# done

cut -f1 mega_motif.txt -d -
tr -d ' ' < mega_motif0.txt > mega_motif00.txt
awk  '$3!=""' mega_motif1.txt > mega_motif01.txt

wc -l mega_motif01.txt


cut -f3,5,6,13 hg38_Tss_coordinates.csv > hg38_Tss_coordinates01.txt
cat hg38_Tss_coordinates01.txt | tail -n +2 >hg38_tss_coord.txt

# ../../../../udd/rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a mega_motif01.txt -b hg38_tss_coord.txt > Motif_BS_Ch_refseq.txt

../../../../udd/rekrg/Tools/bedtools2/bin/bedtools slop  -i mega_motif01.txt -g ../rekrg/Tools/bedtools2/genomes/human.hg38.genome -r 100 -l 100 > megaSlop100.txt


../../../../udd/rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a data/MotifPipeline/ENCODE/wgbsin/GM12878both.txt -b data/MotifPipeline/remap/GM12878_spRE2020.txt > GM12878_BS_Ch.txt

sed -i.bak 's/\t\t/\t/' GM12878_BS_Ch.txt
cut -f1,2,3,4,5,9 GM12878_BS_Ch.txt > GM12878_BS_Ch0.txt
cat GM12878_BS_Ch0.txt | tr "." 0 >GM12878_BS_Ch01.txt

##run in ../../d/tmp/redmo/bench/
../../../../udd/rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a mega_motif01.txt -b GM12878_BS_Ch01.txt > Motif_BS_Ch.txt

../../../../udd/rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a megaSlop100.txt -b GM12878_BS_Ch01.txt > slop100_Motif_BS_Ch.txt

# ../../../../udd/rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a megaSlop250.txt -b GM12878_BS_Ch01.txt > slop250_Motif_BS_Ch.txt

## load into python and take groupby().max() to reduce replicates

chunk_size=100000
for i in range(np.round(226877058/100000).astype('int')):
    data=pd.read_csv('../../d/tmp/redmo/bench/slop100_Motif_BS_Ch.txt',sep='\t',header=None,nrows=chunk_size,skiprows=chunk_size*i)
    data1=data.groupby([0,1,2,10]).max()
#     df.groupby(['Team', 'Pos']).
    data1[3]=data1[3].round(4)
    data1[[3,4,8,9]].to_csv('../../d/tmp/redmo/bench/slop100_Motif_BS_Chmax.txt',sep='\t',index=True,header=False,mode='a')

# orignal motif
../../../../udd/rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a Motif_BS_Chmax.txt -b hg38_tss_coord.txt > Motif_BS_Chmax_refseq.txt
awk '{print($4,"\t",$5,"\t",$6,"\t",$7,"\t",$8,"\t",$12)}' Motif_BS_Chmax_refseq.txt |uniq >>uniq_Motif_BS_Chmax_refseq.txt

# SLOP100
# ../../../../udd/rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a slop100_Motif_BS_Chmax.txt -b hg38_tss_coord.txt > slop100_Motif_BS_Chmax_refseq.txt   ## 1.5Tb
# awk '{print($4,"\t",$5,"\t",$6,"\t",$7,"\t",$8,"\t",$12)}' slop100_Motif_BS_Chmax_refseq.txt |uniq >>uniq_slop100_Motif_BS_Chmax_refseq.txt
#OR

../../../../udd/rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a slop100_Motif_BS_Chmax.txt -b hg38_tss_coord.txt | awk '{print($4,"\t",$5,"\t",$6,"\t",$7,"\t",$8,"\t",$12)}' |uniq > slop100_Motif_BS_Chmax_refseq.txt





