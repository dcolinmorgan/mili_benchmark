# Leveraging methylation information to infer TF binding

>__Conclusions__: Incorporating methylation data improves TF binding inference over standard motif binding profiles, i.e., PWM p-values, and does so in a cell-type specific manner. Although most TFs do not bind to methylated promoter regions, our analysis highlights several exceptions to this rule, indicating that the role of methylation in TF binding is likely cell-type and context specific. <br>
>__Authors__: Daniel Morgan; Kimberly Glass

<space>\
<space>

Clone current version & run [camb_motif_pipeline_gamma.sh](https://github.com/dcolinmorgan/mili_benchmark/blob/master/src/shell/camb_motif_pipeline_gamma.sh)
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

1. Separate hg38 MEME file into individual motif files (730 total PWM files)
1. Run FIMO with threshold=0.0001
1. Convert fimo-output to the BED format
1. For every motif BED file and for every cell line assign motif locations (~) 1 where ChIP for the corresponding TF assayed in that cell line is also observed (otherwise 0). This data to assess prediction within the following motif subsets:
    *  nonCG motif  <br />
        Input: motif locations devoid of CG  <br />
        Output: white in Figure S1  <br />
    *  CG motif  <br />
        Input:CG containing motif locations  <br />
            Output: Figure S1,S2,S7, Confirm baseline motif AUROC (Glass et al. 2013) <br />
    *  CG motif  ∩ WGBS  <br />
        Input: motif locations with WGBS methyl-data  <br />
        Parameters tested: sequence read depth  <br />
        Output: Figure S4 <br />
    *  CG motif ∩ WGBS ∩ array  <br />
        Input: motif locations with both WGBS and methyl array data  <br />
        Output: Figures 2, S3,S5,S6  <br />
        Overlay Yin 2017 data (Figure 3)  <br />
        Additional parameters tested: add +/- 0 : 10kb buffer sizes (Figure 4) 
<space>\
<space>
  
Calculate & summarize AUROCs from intersection files
--------------------------------------------------

```python
   python mili_benchmark/src/python/run_predScore.py -i outdirXX -o outdirXX/test

```

Following this, [the jupyter notebook](https://github.com/dcolinmorgan/mili_benchmark/blob/master/notebook/v8_channing_methyl_benchmark.ipynb) processes AUROCs and figures

Among other things, these checks are performed herewithin:
1. Count multiple CpGs per motif region (varies per TF, ~15) with PWM hit
    1. Compare pairwise distances
    1. Compare mean/median/max/min CpG
1. Compare full WGBS and only subset within array to ChIP, calculate AUROC
    1. threshold WGBS reads >10, calculate AUROC


>Workflow figure from manuscript
>--------------------------------------------------
>![Figure 1. Intersection schema between data modalities](https://github.com/dcolinmorgan/mili_benchmark/blob/master/figures/v3_pdf/fig1.svg)\
> __Figure 1. Intersection schema between data modalities.__ Schematic workflow of bedtools2 intersection calls to calculate the prediction accuracy. The original motif is used as the template, onto which methylation information is supplemented/overwritten to predict ChIP-seq binding activity, where possible. Intersections: 1. Motif to WGBS, 2. Motif-WGBS to methyl array, 3. Motif+WGBS+methyl to ChIP. This narrow view is then expanded by adding buffers before and after motif sites (H0A) and methyl sites (H0B).
