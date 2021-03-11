#!/bin/bash
#$ -cwd
qsub -N "cmpg0" -v b="0",c='A549 K562 GM12878 SKNSH HepG2 HeLa',o='../../../d/tmp/redmo/data/subPipeTest/' -l h_vmem=10G mili_benchmark/camb_motif_pipeline_gamma.sh
qsub -N "cmpg5" -v b="5",c='A549 K562 GM12878 SKNSH HepG2 HeLa',o='../../../d/tmp/redmo/data/subPipeTest/' -l h_vmem=10G mili_benchmark/camb_motif_pipeline_gamma.sh
qsub -N "cmpg10" -v b="10",c='A549 K562 GM12878 SKNSH HepG2 HeLa',o='../../../d/tmp/redmo/data/subPipeTest/' -l h_vmem=10G mili_benchmark/camb_motif_pipeline_gamma.sh
qsub -N "cmpg20" -v b="20",c='A549 K562 GM12878 SKNSH HepG2 HeLa',o='../../../d/tmp/redmo/data/subPipeTest/' -l h_vmem=10G mili_benchmark/camb_motif_pipeline_gamma.sh
qsub -N "cmpg50" -v b="50",c='A549 K562 GM12878 SKNSH HepG2 HeLa',o='../../../d/tmp/redmo/data/subPipeTest/' -l h_vmem=10G mili_benchmark/camb_motif_pipeline_gamma.sh
qsub -N "cmpg100" -v b="100",c='A549 K562 GM12878 SKNSH HepG2 HeLa',o='../../../d/tmp/redmo/data/subPipeTest/' -l h_vmem=10G mili_benchmark/camb_motif_pipeline_gamma.sh
qsub -N "cmpg250" -v b="250",c='A549 K562 GM12878 SKNSH HepG2 HeLa',o='../../../d/tmp/redmo/data/subPipeTest/' -l h_vmem=10G mili_benchmark/camb_motif_pipeline_gamma.sh
qsub -N "cmpg500" -v b="500",c='A549 K562 GM12878 SKNSH HepG2 HeLa',o='../../../d/tmp/redmo/data/subPipeTest/' -l h_vmem=10G mili_benchmark/camb_motif_pipeline_gamma.sh
qsub -N "cmpg1000" -v b="1000",c='A549 K562 GM12878 SKNSH HepG2 HeLa',o='../../../d/tmp/redmo/data/subPipeTest/' -l h_vmem=10G mili_benchmark/camb_motif_pipeline_gamma.sh
qsub -N "cmpg5000" -v b="5000",c='A549 K562 GM12878 SKNSH HepG2 HeLa',o='../../../d/tmp/redmo/data/subPipeTest/' -l h_vmem=10G  mili_benchmark/camb_motif_pipeline_gamma.sh
qsub -N "cmpg10000" -v b="10000",c="A549 K562 GM12878 SKNSH HepG2 HeLa",o="../../../d/tmp/redmo/data/subPipeTest/" -l h_vmem=10G  mili_benchmark/camb_motif_pipeline_gamma.sh
