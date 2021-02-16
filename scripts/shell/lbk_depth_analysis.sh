#!/bin/bash
#$ -cwd
## bash mili_benchmark/lbk_depth_analysis.sh -i'../../../d/tmp/redmo/data/MotifPipeline/lbk_motif_pipeline_gamma0' -b'0' -c'GM12878 SKNSH HepG2 HeLa'
## bash mili_benchmark/lbk_depth_analysis.sh -i'data/MotifPipeline/sthlm_motif_0_QCbeta' -b'0' -c'A549 K562 GM12878 SKNSH HepG2 HeLa'
## bash mili_benchmark/lbk_depth_analysis.sh -i'../../../d/tmp/redmo/data/MotifPipeline/cbus_motif_pipeline_delta0' -b'0' -c'A549 K562 GM12878 SKNSH HepG2 HeLa'

## bash mili_benchmark/lbk_depth_analysis.sh -i'data/MotifPipeline/compare/sthlm_motif_0_QCbeta/gene' -b'0' -c'A549 K562 GM12878 SKNSH HepG2 HeLa'
## bash mili_benchmark/lbk_depth_analysis.sh -i'../../d/tmp/redmo/data/MotifPipeline/cbus_motif_pipeline_delta0/gene' -b'0' -c'A549 K562 GM12878 SKNSH HepG2 HeLa'
## bash mili_benchmark/lbk_depth_analysis.sh -i'../../d/tmp/redmo/data/MotifPipeline/camb_motif_pipeline_gamma100/gene' -b'0' -c'A549 K562 GM12878 SKNSH HepG2 HeLa'



usage=""$camb_depth_analysis" [-h] [-bco] -- function calling bedtools2 to calculate the intersect for ENCODE methylation benchmark with arbitrary buffer length

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
# rm -rf $indir/buffer_analysis.txt
# for buffer in $buffers
# do
	for cell in $cells
	do
		# if [ $buffer != 0 ];then
			# printf "$cell$buffer \t" >> $indir/buffer_analysis.txt
			# eval "cat "$indir"/sthlm_motif_"$buffer"_QCbeta/red/"$cell"*" |uniq -u| wc -l >> $indir/buffer_analysis.txt
      # mkdir $indir/$cell/
			# eval "cat "$indir"/camb_motif_pipeline_gamma"$buffer"/"$cell"*" |uniq -u| wc -l >> $indir/buffer_analysis.txt
			# eval "cat "$indir"/camb_motif_pipeline_gamma"$buffer"/"$cell"*" |uniq -u | wc -l >> $indir/buffer_analysis.txt
		# else
		  # printf "$cell$buffer \t" >> $indir/depth_analysis.txt
			# eval "cat "$indir"/sthlm_motif_"$buffer"_QCbeta/"$cell"*" |uniq -u | wc -l >> $indir/buffer_analysis.txt
		  # eval "cat "$indir"/"$cell"*" |uniq -u |wc -l > $indir/depth_analysis_$cell
      eval "cat "$indir"/"$cell"*" |cut -f6|sort|uniq |wc -l > $indir/gene_analysis_$cell

			# eval "cat "$indir"/camb_motif_"$buffer"_QCbeta/"$cell"*" |uniq -u | wc -l >> $indir/buffer_analysis.txt
		# fi
	# done
  done
      eval "cat "$indir"/*" |cut -f6|sort|uniq|wc -l  > $indir/gene_analysis_all
