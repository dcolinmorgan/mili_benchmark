# mili_benchmark

-work to predict TF binding via various methylation marks, for use in [MILIPEED](https://github.com/dcolinmorgan/netZooPy/tree/milipeed) project within [netZooPy](https://github.com/netZoo/netZooPy)


Implementation of the methylation benchmark as input for MILIPEDE (Methyl-Informed LInk Priming to Integrate Epigenetic Determinants of Expression) network reconstruction approach.\
Authors: Daniel Morgan; Kimberly Glass

<space>\
<space>

this benchmark takes four primary inputs, including the locations of:\
(1) potential transcription factor binding sites, which can be defined by position weight matrices mapped onto the DNA\
(2) 850k Methyl array data\
(3) WGBS methylation data\
(4) ChIP-seq data\

The output are single TF-cell specific intersections of these four data types which are read into jupyternotebooks for analysis.

IMPORTANT: bedtools must be installed in order to run. Please see https://bedtools.readthedocs.io/ for more information. 

<space>\
<space>

Clone current version & run from shell
--------------------------------------------------
```bash
		git clone https://github.com/dcolinmorgan/mili_benchmark ~/src/mili_benchmark
    bash camb_motif_pipeline_gamma.sh -b0 -c'A549 K562 GM12878 SKNSH HepG2 HeLa' -o'outdirXX'

```

This benchmark can be summarized by four main steps:

(1) *STEP 1*: Intersect WGBS with meArray per cell line

   (i) LiftOver meArray to hg38\
   (ii) Separate hg38 meme file into individual motif files (652)\
   (iii) Run FIMO with threshold=0.00001 and auto-background\
   (iv) Convert fimo-output to bedfile output\
 
<space>\
<space>
  
(2) *STEP 2*: For every motif bedfile & for every cell line & for various window sizes:

   (i) Intersect & compare full motif to ChIP (baseline)\
   (ii) Intersect & compare motifs with & w/o CG to ChIP (check for bias)\
  
   (iii) *STEP 3*: Intersect full motif with methyl-data (with various +/-bp windows)\
      (a) Buffer around methyl event (H0 A)\
      (b) Buffer around motif (H0 B)\
   (iv) *STEP 4*: Intersect this file using motif region window with ChIP-seq data\
      (a) 0 where no intersection\
   (v) Count multiple CpGs per motif region (varies per TF, ~15) with PWM hit\
      (a) Compare pairwise distances\
      (b) Compare mean/median/max/min CpG\
   (vi) Compare full WGBS and only subset within array to ChIP, calculate AUROC\
      (a) threshold WGBS reads >10, calculate AUROC\

<space>\
<space>
  
Calculate & summarize AUROCs from intersection files
--------------------------------------------------

```python
  mili_benchmark/run_predScore.py -i outdirXX -o outdirXX/test

```
  
Following this, [the jupyter notebook](https://github.com/dcolinmorgan/mili_benchmark/blob/master/channing_methyl_benchmark.ipynb) processes AUROCs and figures

