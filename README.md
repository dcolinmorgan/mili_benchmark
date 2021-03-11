# Leveraging methylation information to infer TF binding

__Background__: Position weight matrices (PWM) have been used to identify potential locations of transcription factor binding in the genome (motifs). However, this in silico model assumes many simplifications of environmental regulatory processes and does not routinely consider the epigenetic context necessary for TF binding. Furthermore, benchmark studies using ChIP-seq data have demonstrated only modest accuracy when using PWM-based motifs to infer in vivo TF binding. Here we investigate scoring motif locations with methylation information to more accurately infer TF binding. <br>
__Results__: We intersected TF binding motif locations identified using PWMs with methylation information from both whole genome bisulfite and Illumina EPIC array data. We  used the methylation information for six cell lines to score potential TF binding locations and then compared with experimental data (ChIP-seq) to assess whether methylation information can be used to infer binding activity. Our results demonstrate that most TF binding is inhibited in the presence of methylation. This confirms several smaller case studies on the topic of individual TF binding (including NRF1, CTCF), as well as a larger study which infers binding from SELEX data. Finally, we explore how our approach can be expanded to allow for the inference of TF binding locations when methylation information is only proximally available, i.e. no exact bp intersection but allowing for a windowed intersection. <br>
__Conclusions__: Incorporating methylation data improves TF binding inference over standard motif binding profiles, i.e., PWM p-values, and does so in a cell-type specific manner. Although most TFs do not bind to methylated promoter regions, our analysis highlights several exceptions to this rule, indicating that the role of methylation in TF binding is likely cell-type and context specific. <br>
__Authors__: Daniel Morgan; Kimberly Glass

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

Following this, [the jupyter notebook](https://github.com/dcolinmorgan/mili_benchmark/blob/master/notebook/v7_channing_methyl_benchmark.ipynb) processes AUROCs and figures

Among other things, these checks are performed herewithin:
1. Count multiple CpGs per motif region (varies per TF, ~15) with PWM hit
    1. Compare pairwise distances
    1. Compare mean/median/max/min CpG
1. Compare full WGBS and only subset within array to ChIP, calculate AUROC
    1. threshold WGBS reads >10, calculate AUROC


>Workflow figure from manuscript
>--------------------------------------------------
>![Figure 1. Intersection schema between data modalities](https://github.com/dcolinmorgan/mili_benchmark/blob/master/figures/motif_interx_X_link_calls_v6.png)\
> __Figure 1. Intersection schema between data modalities.__ Schematic workflow of bedtools2 intersection calls to calculate the prediction accuracy. The original motif is used as the template, onto which methylation information is supplemented/overwritten to predict ChIP-seq binding activity, where possible. Intersections: 1. Motif to WGBS, 2. Motif-WGBS to methyl array, 3. Motif+WGBS+methyl to ChIP. This narrow view is then expanded by adding buffers before and after motif sites (H0A) and methyl sites (H0B).
