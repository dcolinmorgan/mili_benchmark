#!/bin/bash
#bash count_coverage.sh

motifdir='data/MotifPipeline/compare'  ## chr start stop pwm per gene name/
depths=$(ls $motifdir/|grep sthlm_motif_)
cd $motifdir
for  depth in $depths #
do
motiffiles=$(ls $depth/red/*)
# echo $depth
rm -rf ~/counts.txt
# echo $motiffiles
for tf in $motiffiles #
do
    # tf=HeLa_ESR1
    A=$(eval "wc -l" $tf)
    B=$(eval "uniq "$tf" | wc -l")
    C=$(basename $tf)

    echo $C $B $A $depth>> ~/count.txt
    echo $C $B $A $depth
done
done