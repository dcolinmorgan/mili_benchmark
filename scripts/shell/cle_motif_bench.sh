#!/bin/bash
#$ -cwd
## bash mili_benchmark/scripts/shell/cle_motif_bench.sh -s100 -e200
usage=""$cle_motif_bench" [-h] [-bseo] -- function calling bedtools2 to calculate the intersect for ENCODE methylation benchmark with arbitrary buffer length
where:
  -h  show this help text
  -b  buffer length (integer)
  -s  start
  -e  end
  -o  output directory path"

while getopts ":h:b:s:e:o:" opt; do
  case $opt in
    h) echo "$usage"
       exit
       ;;
    b) buffer=$OPTARG
      ;;
    s) start="$OPTARG"
      ;;
    e) end="$OPTARG" #{OPTARG:-"data/MotifPipeline/camb_motif"$buffer"QC"}
      ;;
    o) end="$OPTARG" #{OPTARG:-"data/MotifPipeline/camb_motif"$buffer"QC"}
      ;;
    \?) echo "Invalid option -$OPTARG" >&3
      ;;
  esac
done


motifdir='/udd/rekrg/MotifScans/MotifScans/hg38'  ## chr start stop pwm per gene name/
motiffiles=$(ls -U $motifdir/* |sed -n "$start,$end p")

# echo $depth
# rm -rf ~/counts.txt
# echo $motiffiles
for tf in $motiffiles # ${motiffiles[@]:$start:$end} #
do
  # TF=$(eval "echo "$tf" | cut -d / -f6| cut -d _ -f1")
  TF=$(eval "cut -f2 $tf"|head -2 | tail -1|cut -f1 -d _|sed 's/[)(]//g')
  awk '{print($3,"\t",$4,"\t",$5,"\t",$7,"\t",$2)}' $tf |tail -n +2 > '../../d/tmp/redmo/bench/alt/'$TF'A.txt' #>mega_motif.txt
  cut -f1 '../../d/tmp/redmo/bench/alt/'$TF'A.txt' -d _ > '../../d/tmp/redmo/bench/alt/'$TF'B.txt' # cut -f3,4,5,7,3 $tf
  tr -d ' ' < '../../d/tmp/redmo/bench/alt/'$TF'B.txt' > '../../d/tmp/redmo/bench/alt/'$TF'C.txt'
  awk  '$3!=""' '../../d/tmp/redmo/bench/alt/'$TF'C.txt' > '../../d/tmp/redmo/bench/alt/'$TF'D.txt'
  tr -d ' ' < '../../d/tmp/redmo/bench/alt/'$TF'D.txt' > '../../d/tmp/redmo/bench/alt/'$TF'E.txt'

  printf ""$TF"\n"
  # cut -f3,5,6,13 ../../d/tmp/redmo/bench/hg38_Tss_coordinates.csv > ../../d/tmp/redmo/bench/hg38_Tss_coordinates01.txt
  # cat ../../d/tmp/redmo/bench/hg38_Tss_coordinates01.txt | tail -n +2 >../../d/tmp/redmo/bench/hg38_tss_coord.txt
  eval "~/../../../../udd/rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a ../../d/tmp/redmo/bench/alt/"$TF"E.txt -b ../../d/tmp/redmo/bench/hg38_tss_coord.txt > ../../d/tmp/redmo/bench/alt/"$TF"_TSS.txt"
  # if buffer!=0;then
  #   eval "~/../../../../udd/rekrg/Tools/bedtools2/bin/bedtools slop  -i mega_motif01.txt -g ../rekrg/Tools/bedtools2/genomes/human.hg38.genome -r buffer -l buffer > ../../d/tmp/redmo/bench/alt/megaSlop_"$buffer".txt"
  #   sed -i.bak 's/\t\t/\t/' ../../d/tmp/redmo/bench/alt/megaSlop_"$buffer".txt
  #     cut -f1,2,3,4,5,11 ../../d/tmp/redmo/bench/alt/megaSlop_"$buffer".txt > ../../d/tmp/redmo/bench/alt/megaSlop_"$buffer".txt
 #      eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -a $bench/overlap_tmpBB_$cell"_"$buffer$out_fmt -b $bench/overlap_tmpBB_$cell"_"0$out_fmt " >  $outdir/test/$cell"_"$gene

  #     eval "~/../../../../udd/rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a ../../d/tmp/redmo/bench/alt/Motif_TSS.txt -b ../../d/tmp/redmo/bench/GM12878_BS_Ch.txt > Motif_TSS_BS_Ch.txt"
 
  # else
    eval "~/../../../../udd/rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a ../../d/tmp/redmo/bench/alt/Motif_TSS.txt -b ../../d/tmp/redmo/bench/GM12878_BS_Ch01.txt "> "../../d/tmp/redmo/bench/alt/"$TF".txt"
  # fi
  rm -rf "../../d/tmp/redmo/bench/alt/"$TF"_TSS.txt" "../../d/tmp/redmo/bench/alt/"$TF"A.txt" "../../d/tmp/redmo/bench/alt/"$TF"B.txt" "../../d/tmp/redmo/bench/alt/"$TF"C.txt" "../../d/tmp/redmo/bench/alt/"$TF"D.txt" "../../d/tmp/redmo/bench/alt/"$TF"E.txt"
done
# done
c

# source /proj/relibs/relib00/conda/bin/activate
# source activate mypy3
# chmod +x mili_benchmark/scripts/python/milipede_red.py
# python mili_benchmark/scripts/python/milipede_red.py




# cut -f1 mega_motif.txt -d -
# tr -d ' ' < mega_motif0.txt > mega_motif00.txt
# awk  '$3!=""' mega_motif1.txt > mega_motif01.txt


# cut -f3,5,6,13 hg38_Tss_coordinates.csv > hg38_Tss_coordinates01.txt
# cat hg38_Tss_coordinates01.txt | tail -n +2 >hg38_tss_coord.txt

# # ../../../../udd/rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a mega_motif01.txt -b hg38_tss_coord.txt > Motif_BS_Ch_refseq.txt

# ../../../../udd/rekrg/Tools/bedtools2/bin/bedtools slop  -i mega_motif01.txt -g ../rekrg/Tools/bedtools2/genomes/human.hg38.genome -r 100 -l 100 > megaSlop100.txt


# ../../../../udd/rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a data/MotifPipeline/ENCODE/wgbsin/GM12878both.txt -b data/MotifPipeline/remap/GM12878_spRE2020.txt > GM12878_BS_Ch.txt

# sed -i.bak 's/\t\t/\t/' GM12878_BS_Ch.txt
# cut -f1,2,3,4,5,9 GM12878_BS_Ch.txt > GM12878_BS_Ch0.txt
# cat GM12878_BS_Ch0.txt | tr "." 0 >GM12878_BS_Ch01.txt

# ##run in ../../d/tmp/redmo/bench/
# ../../../../udd/rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a mega_motif01.txt -b GM12878_BS_Ch01.txt > Motif_BS_Ch.txt

# ../../../../udd/rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a megaSlop100.txt -b GM12878_BS_Ch01.txt > slop100_Motif_BS_Ch.txt

# # ../../../../udd/rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a megaSlop250.txt -b GM12878_BS_Ch01.txt > slop250_Motif_BS_Ch.txt


# sort -u -k1,4 -k11 test.txt

# ../../../../udd/rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a Motif_BS_Ch.txt -b hg38_tss_coord.txt > Motif_BS_Ch_refseq.txt

# ../../../../udd/rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a slop100_Motif_BS_Ch.txt -b hg38_tss_coord.txt > slop100_Motif_BS_Ch_refseq.txt





