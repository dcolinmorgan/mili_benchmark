library(minfi)
# BiocManager::install("minfiData")
library(minfiData)
library(IlluminaHumanMethylationEPICmanifest)#,
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
RGSet<-minfi::read.metharray.exp('data/MotifPipeline/ENCODE/methyl_array/850k/')
# phenoData <- pData(RGSet)
manifest <- getManifest(RGSet)
head(getProbeInfo(manifest))
MSet <- preprocessRaw(RGSet) 
MSet
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
RSet
beta <- getBeta(RSet)
GRset <- mapToGenome(RSet)
GRset

##QC
RGSet <- preprocessQuantile(RGSet, fixOutliers = TRUE,removeBadSamples = TRUE, badSampleCutoff = 10.5,quantileNormalize = TRUE, stratified = TRUE, mergeManifest = FALSE, sex = NULL)
# MSet.swan <- preprocessSWAN(RGSet)
# RSet2 <- ratioConvert(MSet.swan)

snps <- getSnpInfo(RGSet)
head(snps,10)
RGSet <- addSnpInfo(GRset) #GRset
##filter out snps
RGSet <- dropLociWithSnps(RGSet, snps=c("SBE","CpG"), maf=0)

###filter out cross-reactive probes
devtools::install_github("markgene/maxprobes")
RGSet<-maxprobes::dropXreactiveLoci(RGSet)


##else comment above out
beta <- getBeta(RGSet)
M <- getM(GRset)
CN <- getCN(GRset)
sampleNames <- sampleNames(GRset)
probeNames <- featureNames(GRset)
pheno <- pData(GRset)
gr <- granges(GRset)
annotation <- getAnnotation(RGSet)
qc <- getQC(MSet)
head(qc)
plotQC(qc)
densityPlot(beta, sampGroups = GRset@colData@rownames,main = 'preQC_beta',legend = None)
densityBeanPlot(beta, sampGroups = GRset@colData@rownames,main = 'preQC_beta')
controlStripPlot(RGSet, controls="BISULFITE CONVERSION II")
qcReport(RGSet, pdf= "qcReport.pdf",sampNames = MSet@colData@rownames,sampGroups = MSet@colData@rownames)

write.table(beta[,1],'data/MotifPipeline/ENCODE/methyl_array/A-549b.txt',sep ='\t',quote = FALSE,dec = '.',col.names = FALSE)
write.table(beta[,2],'data/MotifPipeline/ENCODE/methyl_array/GM12878b.txt',sep ='\t',quote = FALSE,dec = '.',col.names = FALSE)
write.table(beta[,3],'data/MotifPipeline/ENCODE/methyl_array/HeLa-S3b.txt',sep ='\t',quote = FALSE,dec = '.',col.names = FALSE)
write.table(beta[,4],'data/MotifPipeline/ENCODE/methyl_array/Hep-G2b.txt',sep ='\t',quote = FALSE,dec = '.',col.names = FALSE)
write.table(beta[,5],'data/MotifPipeline/ENCODE/methyl_array/K562b.txt',sep ='\t',quote = FALSE,dec = '.',col.names = FALSE)
write.table(beta[,6],'data/MotifPipeline/ENCODE/methyl_array/SKNSHb.txt',sep ='\t',quote = FALSE,dec = '.',col.names = FALSE)
write.table(RGSet@rowRanges@ranges@start,'data/MotifPipeline/ENCODE/methyl_array/startb.txt',sep ='\t',quote = FALSE,dec = '.',row.names=FALSE,col.names = FALSE)
write.table(annotation@listData$chr,'data/MotifPipeline/ENCODE/methyl_array/chrb.txt',sep ='\t',quote = FALSE,dec = '.',row.names=FALSE,col.names = FALSE)
write.table(annotation$Relation_to_Island,'data/MotifPipeline/ENCODE/methyl_array/R2island.txt',sep ='\t',quote = FALSE,dec = '.',row.names=FALSE,col.names = FALSE)
write.table(annotation@listData$UCSC_RefGene_Group,'data/MotifPipeline/ENCODE/methyl_array/gene_bodyA.txt',sep ='\t',quote = FALSE,dec = '.',row.names=FALSE,col.names = FALSE)
write.table(annotation@listData$GencodeBasicV12_Group,'data/MotifPipeline/ENCODE/methyl_array/gene_bodyB.txt',sep ='\t',quote = FALSE,dec = '.',row.names=FALSE,col.names = FALSE)
write.table(annotation@listData$Name,'data/MotifPipeline/ENCODE/methyl_array/cpg.txt',sep ='\t',quote = FALSE,dec = '.',row.names=FALSE,col.names = FALSE)

##bedtools intersect to add gene annotation to cg location info
write.table(annotation@listData$Gene,'data/MotifPipeline/ENCODE/methyl_array/gene.txt',sep ='\t',quote = FALSE,dec = '.',row.names=FALSE,col.names = FALSE)

##python convert empty rows to 0s

# ##bash to format as bed file cd data/MotifPipeline/ENCODE/methyl_array/
# paste chrb.txt startb.txt A-549b.txt R2island.txt gene_bodyE.txt gene_bodyF.txt| awk '{print $1,$2,$2+1,$4,$5,$6,$7}' OFS='\t'> A-549_MeArrayc.txt
# paste chrb.txt startb.txt GM12878b.txt R2island.txt gene_bodyE.txt gene_bodyF.txt| awk '{print $1,$2,$2+1,$4,$5,$6,$7}' OFS='\t'> GM12878_MeArrayc.txt
# paste chrb.txt startb.txt HeLa-S3b.txt R2island.txt gene_bodyE.txt gene_bodyF.txt| awk '{print $1,$2,$2+1,$4,$5,$6,$7}' OFS='\t'> HeLa-S3_MeArrayc.txt
# paste chrb.txt startb.txt Hep-G2b.txt R2island.txt gene_bodyE.txt gene_bodyF.txt| awk '{print $1,$2,$2+1,$4,$5,$6,$7}' OFS='\t'> Hep-G2_MeArrayc.txt
# paste chrb.txt startb.txt K562b.txt R2island.txt gene_bodyE.txt gene_bodyF.txt| awk '{print $1,$2,$2+1,$4,$5,$6,$7}' OFS='\t'> K-562_MeArrayc.txt
# paste chrb.txt startb.txt SKNSHb.txt R2island.txt gene_bodyE.txt gene_bodyF.txt| awk '{print $1,$2,$2+1,$4,$5,$6,$7}' OFS='\t'> SK-N-SH_MeArrayc.txt
# # 
# #
# # ##map from 19 to 38 cd to data
# ./data/liftOver data/MotifPipeline/ENCODE/methyl_array/A-549_MeArrayc.txt data/hg19ToHg38.over.chain data/MotifPipeline/ENCODE/methyl_array/A549_MeArrayHG38d.txt unmapped
# ./data/liftOver data/MotifPipeline/ENCODE/methyl_array/GM12878_MeArrayc.txt data/hg19ToHg38.over.chain data/MotifPipeline/ENCODE/methyl_array/GM12878_MeArrayHG38d.txt unmapped
# ./data/liftOver data/MotifPipeline/ENCODE/methyl_array/HeLa-S3_MeArrayc.txt data/hg19ToHg38.over.chain data/MotifPipeline/ENCODE/methyl_array/HeLa_MeArrayHG38d.txt unmapped
# ./data/liftOver data/MotifPipeline/ENCODE/methyl_array/Hep-G2_MeArrayc.txt data/hg19ToHg38.over.chain data/MotifPipeline/ENCODE/methyl_array/HepG2_MeArrayHG38d.txt unmapped
# ./data/liftOver data/MotifPipeline/ENCODE/methyl_array/K-562_MeArrayc.txt data/hg19ToHg38.over.chain data/MotifPipeline/ENCODE/methyl_array/K562_MeArrayHG38d.txt unmapped
# ./data/liftOver data/MotifPipeline/ENCODE/methyl_array/SK-N-SH_MeArrayc.txt data/hg19ToHg38.over.chain data/MotifPipeline/ENCODE/methyl_array/SKNSH_MeArrayHG38d.txt unmapped
# # 
# 
# 
# cut -f1,2,3,4,5,7 data/MotifPipeline/ENCODE/methyl_array/A549_MeArrayHG38d.txt > data/MotifPipeline/ENCODE/methyl_array/A549_MeArrayHG38e.txt
# cut -f1,2,3,4,5,7 data/MotifPipeline/ENCODE/methyl_array/GM12878_MeArrayHG38d.txt > data/MotifPipeline/ENCODE/methyl_array/GM12878_MeArrayHG38e.txt
# cut -f1,2,3,4,5,7 data/MotifPipeline/ENCODE/methyl_array/HeLa_MeArrayHG38d.txt > data/MotifPipeline/ENCODE/methyl_array/HeLa_MeArrayHG38e.txt
# cut -f1,2,3,4,5,7 data/MotifPipeline/ENCODE/methyl_array/HepG2_MeArrayHG38d.txt > data/MotifPipeline/ENCODE/methyl_array/HepG2_MeArrayHG38e.txt
# cut -f1,2,3,4,5,7 data/MotifPipeline/ENCODE/methyl_array/K562_MeArrayHG38d.txt > data/MotifPipeline/ENCODE/methyl_array/K562_MeArrayHG38e.txt
# cut -f1,2,3,4,5,7 data/MotifPipeline/ENCODE/methyl_array/SKNSH_MeArrayHG38d.txt > data/MotifPipeline/ENCODE/methyl_array/SKNSH_MeArrayHG38e.txt

# ##build shuffled versions
# cut -f4 A-549_MeArrayHG38c.txt|cut -f4 |shuf> shuf_A549.txt
# paste A-549_MeArrayHG38c.txt shuf_A549.txt | awk '{print $1,$2,$3,$4,$5,$6,$7}' OFS='\t' > A-549_shufMeArrayHG38c.txt
# cut -f4 GM12878_MeArrayHG38c.txt|cut -f4 |shuf> shuf_GM12878.txt
# paste GM12878_MeArrayHG38c.txt shuf_GM12878.txt | awk '{print $1,$2,$3,$4,$5,$6,$7}' OFS='\t' > GM12878_shufMeArrayHG38c.txt
# cut -f4 HeLa-S3_MeArrayHG38c.txt|cut -f4 |shuf> shuf_HeLa-S3.txt
# paste HeLa-S3_MeArrayHG38c.txt shuf_HeLa-S3.txt | awk '{print $1,$2,$3,$4,$5,$6,$7}' OFS='\t' > HeLa-S3_shufMeArrayHG38c.txt
# cut -f4 Hep-G2_MeArrayHG38c.txt|cut -f4 |shuf> shuf_Hep-G2.txt
# paste Hep-G2_MeArrayHG38c.txt shuf_Hep-G2.txt | awk '{print $1,$2,$3,$4,$5,$6,$7}' OFS='\t' > Hep-G2_shufMeArrayHG38c.txt
# cut -f4 K-562_MeArrayHG38c.txt|cut -f4 |shuf> shuf_K-562.txt
# paste K-562_MeArrayHG38c.txt shuf_K-562.txt | awk '{print $1,$2,$3,$4,$5,$6,$7}' OFS='\t' > K-562_shufMeArrayHG38c.txt
# cut -f4 SK-N-SH_MeArrayHG38c.txt|cut -f4 |shuf> shuf_SK-N-SH.txt
# paste SK-N-SH_MeArrayHG38c.txt shuf_SK-N-SH.txt | awk '{print $1,$2,$3,$4,$5,$6,$7}' OFS='\t' > SK-N-SH_shufMeArrayHG38c.txt
# 
#

### create milipeed motifs
# paste cpg.txt A-549b.txt | awk '{print $1,$2}' OFS='\t'> A-549_Metif.txt
# paste cpg.txt GM12878b.txt | awk '{print $1,$2}' OFS='\t'> GM12878_Metif.txt
# paste cpg.txt HeLa-S3b.txt | awk '{print $1,$2}' OFS='\t'> HeLa_Metif.txt
# paste cpg.txt Hep-G2b.txt | awk '{print $1,$2}' OFS='\t'> HepG2_Metif.txt
# paste cpg.txt K562b.txt | awk '{print $1,$2}' OFS='\t'> K562_Metif.txt
# paste cpg.txt SKNSHb.txt | awk '{print $1,$2}' OFS='\t'> SKNSH_Metif.txt

### create gold standard chip motif for comparison
# cd data/MotifPipeline/ENCODE/methyl_array/
# paste chrb.txt startb.txt A-549b.txt | awk '{print $1,$2,$2+1,$3,$4}' OFS='\t'>> A-549_Metifb.txt
# eval "~/../rekrg/Tools/bedtools2/bin/bedtools slop  -i data/MotifPipeline/remap/A549_spRE2020.txt -g ~/../rekrg/Tools/bedtools2/genomes/human.hg38.genome -r 1000 -l 1000" > data/MotifPipeline/remap/slopA549.txt
# eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a data/MotifPipeline/remap/slopA549.txt -b data/MotifPipeline/ENCODE/methyl_array/A-549_Metifb.txt " > data/MotifPipeline/A549GOLD_1000kb.txt
# eval "~/../rekrg/Tools/bedtools2/bin/bedtools intersect -wa -wb -a data/MotifPipeline/remap/A549_spRE2020.txt -b data/MotifPipeline/ENCODE/methyl_array/A-549_Metifb.txt " > data/MotifPipeline/A549GOLD.txt

#