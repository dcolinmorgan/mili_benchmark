#!/bin/bash
#$ -cwd
## bash mili_benchmark/gene_count.sh


cg_gene_loc='data/MotifPipeline/ENCODE/methyl_array/cg_gene_loc.txt'

motifdirA='data/MotifPipeline/compare/sthlm_motif_0_QCbeta'
motifsA=$(ls -f $motifdirA/*)
outdirA=$motifdirA/gene
mkdir -p $outdirA

motifdirB='../../d/tmp/redmo/data/MotifPipeline/cbus_motif_pipeline_delta0'
motifsB=$(ls -f $motifdirB/*)
outdirB=$motifdirB/gene
mkdir -p $outdirB

motifdirC='../../d/tmp/redmo/data/MotifPipeline/camb_motif_pipeline_gamma100'
motifsC=$(ls -f $motifdirC/*)
outdirC=$motifdirC/gene
mkdir -p $outdirC

# for TF in $motifsA
# do
# 	tf=$(eval "echo "$TF" | cut -d / -f5")

# 	printf "running "$tf" in A \n"
# 	eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect  -wa -wb -a $cg_gene_loc -b $TF" > $outdirA/$tf
# done

# for TF in $motifsB
# do
# 	tf=$(eval "echo "$TF" | cut -d / -f9")

# 	printf "running "$tf" in B \n"
# 	eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect  -wa -wb -a $cg_gene_loc -b $TF" > $outdirB/$tf
# done

for TF in $motifsC
do
	tf=$(eval "echo "$TF" | cut -d / -f9")

	printf "running "$tf" in C \n"
	eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect  -wa -wb -a $cg_gene_loc -b $TF" > $outdirC/$tf
done