#!/bin/bash
#$ -cwd
## bash mili_benchmark/scripts/shell/cle_motif_bench.sh -b0 -s1 -e100
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
    e) end="$OPTARG" 
      ;;
    o) output="$OPTARG" 
      ;;
    \?) echo "Invalid option -$OPTARG" >&3
      ;;
  esac
done


motifdir='/udd/rekrg/MotifScans/MotifScans/hg38'  ## chr start stop pwm per gene name/
motiffiles=$(ls -U $motifdir/* |sed -n "$start,$end p")

output='../../d/tmp/redmo/bench/'
output=$output$buffer
# rm -rf $output
# mkdir $output
if [ ! -d "$output" ]; then
  mkdir $output
fi
# echo $depth
# rm -rf ~/counts.txt
# echo $motiffiles
for tf in $motiffiles # ${motiffiles[@]:$start:$end} #
do
  # TF=$(eval "echo "$tf" | cut -d / -f6| cut -d _ -f1")
  TF=$(eval "cut -f2 $tf"|head -2 | tail -1|cut -f1 -d _|sed 's/[)(]//g')
  awk '{print($3,"\t",$4,"\t",$5,"\t",$7,"\t",$9,"\t",$2)}' $tf |tail -n +2 > $output/$TF'_A.txt' #>mega_motif.txt
  cut -f1 $output/$TF'_A.txt' -d _ > $output/$TF'_B.txt' # cut -f3,4,5,7,3 $tf
  tr -d ' ' < $output/$TF'_B.txt' > $output/$TF'_C.txt'
  awk  '$3!=""' $output/$TF'_C.txt' > $output/$TF'_D.txt'
  tr -d ' ' < $output/$TF'_D.txt' > $output/$TF'_E.txt'

  printf ""$TF"\n"
  # cut -f3,5,6,13 ../../d/tmp/redmo/bench/hg38_Tss_coordinates.csv > ../../d/tmp/redmo/bench/hg38_Tss_coordinates01.txt
  # cat ../../d/tmp/redmo/bench/hg38_Tss_coordinates01.txt | tail -n +2 >../../d/tmp/redmo/bench/hg38_tss_coord.txt
  eval "~/../../../../udd/rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a " $output/$TF"_E.txt -b ../../d/tmp/redmo/bench/hg38_tss_coord.txt" > $output/$TF"_TSS.txt"
  
  # ../../../../udd/rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a data/MotifPipeline/ENCODE/wgbsin/GM12878both.txt -b data/MotifPipeline/remap/GM12878_spRE2020.txt > GM12878_BS_Ch.txt

  # sed -i.bak 's/\t\t/\t/' GM12878_BS_Ch.txt
  # cut -f1,2,3,4,5,9 GM12878_BS_Ch.txt > GM12878_BS_Ch0.txt
  # cat GM12878_BS_Ch0.txt | tr "." 0 >GM12878_BS_Ch01.txt

  if [ $buffer != 0 ];then
    printf "adding +/-"$buffer" bp buffer around motif region for "$TF" \n"
    eval "~/../../../../udd/rekrg/Tools/bedtools2/bin/bedtools slop  -i " $output/$TF"_TSS.txt -g ../rekrg/Tools/bedtools2/genomes/human.hg38.genome -r "$buffer" -l "$buffer" "> $output/$TF"_megaSlop_"$buffer".txt"
    
    printf "...intersecting buffered motif with "$TF" WGBS-predictions with "$buffer" bp buffer \n"
    eval "~/../../../../udd/rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a " $output/$TF"_megaSlop_"$buffer".txt -b ../../d/tmp/redmo/bench/GM12878_BS_Ch01.txt "> $output/$TF".txt"
    rm -rf $output"megaSlop_"$buffer"_B.txt"
  else
    printf "...intersecting buffered motif with "$TF" WGBS-predictions with 0bp buffer \n"
    eval "~/../../../../udd/rekrg/Tools/bedtools2/bin/bedtools intersect -wao -a " $output/$TF"_TSS.txt -b ../../d/tmp/redmo/bench/GM12878_BS_Ch01.txt "> $output/$TF".txt"
  fi
  rm -rf $output/$TF"_TSS.txt" $output/$TF"_megaSlop_"$buffer".txt" $output/$TF"_A.txt" $output/$TF"_B.txt" $output/$TF"_C.txt" $output/$TF"_D.txt" $output/$TF"_E.txt"
done
# done


source /proj/relibs/relib00/conda/bin/activate
source activate mypy3
chmod +x mili_benchmark/scripts/python/milipede_red.py
python mili_benchmark/scripts/python/milipede_red.py $output test


####TEST #####



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





