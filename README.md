# mili_benchmark

Implementation of the methylation benchmark as input for MILIPEDE (Methyl-Informed LInk Priming to Integrate Epigenetic Determinants of Expression) network reconstruction approach.\


Investigation to predict TF binding via various methylation marks, for use in [MILIPEED](https://github.com/dcolinmorgan/netZooPy/tree/milipeed) project within [netZooPy](https://github.com/netZoo/netZooPy)

Authors: Daniel Morgan; Kimberly Glass

<space>\
<space>

Clone current version & run from shell
--------------------------------------------------
```bash
    git clone https://github.com/dcolinmorgan/mili_benchmark ~/src/mili_benchmark
    bash camb_motif_pipeline_gamma.sh -b0 -c'A549 K562 GM12878 SKNSH HepG2 HeLa' -o'outdirXX'

```

This benchmark can be summarized by three main steps (numbered) and takes four primary inputs, including the locations of:
  1. potential transcription factor binding sites, which can be defined by position weight matrices mapped onto the DNA
  1. 850k Methyl array data
  1. WGBS methylation data
  1. ChIP-seq data

The output are TF-cell specific intersections of these four data types which are read into jupyternotebooks for analysis. Several other steps provide checks against bias for supplemental figures (bullet points)

IMPORTANT: bedtools must be installed in order to run. Please see https://bedtools.readthedocs.io/ for more information. 

1. Intersect WGBS with meArray per cell line
    1. LiftOver meArray to hg38
    1. Separate hg38 meme file into individual motif files (652)
    1. Run FIMO with threshold=0.00001 and auto-background
    1. Convert fimo-output to bedfile output
  
     __For every motif bedfile & for every cell line & for various window sizes:__ \
          i. Intersect & compare full motif to ChIP (baseline)\
          ii. Intersect & compare motifs with & w/o CG to ChIP (check for bias)\
      
1. Intersect full motif with methyl-data (with various +/-bp windows)
    1. Buffer around methyl event (H0 A)
    1. Buffer around motif (H0 B)
1. Intersect this file using motif region window with ChIP-seq data
    1. 0 where no intersection

<space>\
<space>
  
Calculate & summarize AUROCs from intersection files
--------------------------------------------------

```python
   python mili_benchmark/run_predScore.py -i outdirXX -o outdirXX/test

```

Following this, [the jupyter notebook](https://github.com/dcolinmorgan/mili_benchmark/blob/master/channing_methyl_benchmark.ipynb) processes AUROCs and figures

Among other things, these checks are performed herewithin:
1. Count multiple CpGs per motif region (varies per TF, ~15) with PWM hit
    1. Compare pairwise distances
    1. Compare mean/median/max/min CpG
1. Compare full WGBS and only subset within array to ChIP, calculate AUROC
    1. threshold WGBS reads >10, calculate AUROC


>Workflow figure from manuscript
>--------------------------------------------------
>![Figure 1. Intersection schema between data modalities](https://github.com/dcolinmorgan/mili_benchmark/blob/master/figures/motif_interx%20X%20link_calls_v5.png)\
> __Figure 1. Intersection schema between data modalities.__ Schematic workflow of bedtools2 intersection calls to calculate the prediction accuracy. The original motif is used as the template, onto which methylation information is supplemented/overwritten to predict ChIP-seq binding activity, where possible. Intersections: 1. Motif to WGBS, 2. Motif-WGBS to methyl array, 3. Motif+WGBS+methyl to ChIP. This narrow view is then expanded by adding buffers before and after motif sites (H0A) and methyl sites (H0B).
