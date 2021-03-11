# bash data/MotifPipeline/ENCODE/encode_bench_mega.sh
# rm -rf data/MotifPipeline/sthlm_motif_0_QCbeta/A549_bench.txt
motifdir='data/MotifPipeline/sthlm_motif_0_QCbeta'
motifs=$(ls $motifdir/A549*)

for f in $motifs
do
	echo "$(awk '{print $0, split(FILENAME, parts, "/"),parts[4]}' $f)" > $f #test/$c
	printf "formatted $f \n"

done

cat data/MotifPipeline/sthlm_motif_0_QCbeta/A549* > data/MotifPipeline/sthlm_motif_0_QCbeta/A549_bench.txt