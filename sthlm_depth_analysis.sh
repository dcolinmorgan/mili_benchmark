#!/bin/bash
#$ -cwd
## bash mili_benchmark/sthlm_depth_analysis.sh -i'data/MotifPipeline' -b'0 5 10 20 50 100 250 500 1000 5000 10000 50000' -c'A549 K562 GM12878 SKNSH HepG2 HeLa'
usage=""$sthlm_depth_analysis" [-h] [-bco] -- function calling bedtools2 to calculate the intersect for ENCODE methylation benchmark with arbitrary buffer length

where:
  -h  show this help text
  -i  indir to analyze
  -b  buffer length (integer)
  -c  cell lines (string array 'A549 K562 ...')"

while getopts ":h:i:b:c:" opt; do
  case $opt in
    h) echo "$usage"
       exit
       ;;
    i) indir=$OPTARG
      ;;
    b) buffers=$OPTARG
      ;;
    c) cells="$OPTARG"
      ;;
    \?) echo "Invalid option -$OPTARG" >&3
      ;;
  esac
done

# indir='data/MotifPipeline'
# buffer=[0,5,10,20,50,100,250,500,1000,5000,10000,50000]
# cells=[A549,K562,GM12878,SKNSH,HepG2,HeLa]
rm -rf $indir/buffer_analysis.txt
for buffer in $buffers
do
	for cell in $cells
	do
		if [ $buffer != 0 ];then
			printf "$cell$buffer \t" >> $indir/buffer_analysis.txt
			eval "cat "$indir"/sthlm_motif_"$buffer"_QCbeta/red/"$cell"*" |uniq -u | wc -l >> $indir/buffer_analysis.txt
		else
			printf "$cell$buffer \t" >> $indir/buffer_analysis.txt
			eval "cat "$indir"/sthlm_motif_"$buffer"_QCbeta/"$cell"*" |uniq -u | wc -l >> $indir/buffer_analysis.txt
		fi
	done
done