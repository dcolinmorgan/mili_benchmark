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
  1. potential [transcription factor binding sites](https://genome.ucsc.edu/), which can be defined by position weight matrices mapped onto the DNA
  1. [850k Methyl array data](https://www.encodeproject.org/matrix/?type=Experiment&status=released&award.project=ENCODE&files.platform.term_name=Illumina+Infinium+Methylation+EPIC+BeadChip&biosample_ontology.term_name=A549&biosample_ontology.term_name=K562&biosample_ontology.term_name=GM12878&biosample_ontology.term_name=HeLa-S3&biosample_ontology.term_name=HepG2&biosample_ontology.term_name=SK-N-SH&assay_title=DNAme+array)
  1. [WGBS methylation data](https://www.encodeproject.org/matrix/?type=Experiment&status=released&assay_slims=DNA+methylation&biosample_ontology.classification=cell+line&assay_title=WGBS&biosample_ontology.term_name=A549&biosample_ontology.term_name=K562&biosample_ontology.term_name=GM12878&biosample_ontology.term_name=HeLa-S3&biosample_ontology.term_name=HepG2&biosample_ontology.term_name=SK-N-SH)
  1. [ChIP-seq data](http://remap.univ-amu.fr/)

The output are TF-cell specific intersections of these four data types which are read into jupyternotebooks for analysis. Several other steps provide checks against bias for supplemental figures (bullet points)

IMPORTANT: bedtools must be installed in order to run. Please see https://bedtools.readthedocs.io/ for more information. 

1. Separate hg38 meme file into individual motif files (730)
1. Run FIMO with threshold=0.00001 and auto-background
1. Convert fimo-output to bedfile output
1. For every motif bedfile & for every cell line score motif locations 1 where ∩Ch observed, otherwise 0 within following operations: <br>
    1. nonCG motif + ∩Ch
        1. Input: Intersect motif locations devoid of CG to ChIP
    1. CG motif + ∩Ch 
        1.Input: Intersect CG containing motif to ChIP
        1. Parameters: motif length, CG count per motif
        1. Output: Figure SM1-SM3, Confirm baseline motif AUROC [23]
    1. CG motif + ∩WB +∩Ch
        1. Input: Intersect full motif with WGBS methyl-data
        1. Parameter: sequence read depth
        1. Output: Figure S4
    1. CG motif + ∩WB +∩M+∩Ch
        1. Input: Intersect full motif with methyl array data
        1. Parameter: add +/- 0-10kb buffer sizes (Figure SM4)
        1. Output: Main analysis (Figures 2-5, S1-S3)
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
>![Figure 1. Intersection schema between data modalities](https://github.com/dcolinmorgan/mili_benchmark/blob/master/figures/motif_interx%20X%20link_calls_v6.png)\
> __Figure 1. Intersection schema between data modalities.__ Schematic workflow of bedtools2 intersection calls to calculate the prediction accuracy. The original motif is used as the template, onto which methylation information is supplemented/overwritten to predict ChIP-seq binding activity, where possible. Intersections: 1. Motif to WGBS, 2. Motif-WGBS to methyl array, 3. Motif+WGBS+methyl to ChIP. This narrow view is then expanded by adding buffers before and after motif sites (H0A) and methyl sites (H0B).
