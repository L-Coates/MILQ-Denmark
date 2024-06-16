#This script includes the steps involved in assessing the variation in infant stool [DNA]
#and read depth by extraction batch, bead basher type (old vs. new), plate batch, or sequencing batch. 
#This script also includes the steps involved in determining if beta-diversity
#varied by DNA extraction batch, PCR/plate batch, sequencing batch, read depth,
#or [DNA]. 


#STEP 1: Determine if [DNA] varied by extraction batch or bead basher type

#read in the metadata table which has the DNA yield and batch and bead basher information
metadata <- read.table("../qiime_metadata.tsv", header = T)
print(colnames(metadata))

#only look at the Denmark samples and the zymo mock samples. 
metadata.v2 <- rbind(metadata[metadata$sample.type=="stool" & metadata$country=="Denmark",], metadata[metadata$sample.type=="mock",])
#looking at the number of samples in each DNA extraction batch among this subset
table(metadata.v2$dna.extraction.batch)
#there are only five samples in extract batch 23 and only 9 in extraction batch 24
#while there are 22 - 46 samples in all the other batches. Therefore I will compare
#results with and without batches #23 and #24 in the dataset.  
metadata.v3 <- metadata.v2[metadata.v2$dna.extraction.batch !="23",]
metadata.v3 <- metadata.v3[metadata.v3$dna.extraction.batch !="24",]
unique(metadata.v3$dna.extraction.batch)

#plot [DNA] by extraction batch
library(ggplot2)
#all batches
metadata.v2$dna.extraction.batch <- as.character(metadata.v2$dna.extraction.batch)
ggplot(aes(y=`nanodrop.dna.concentration.ng.per.uL`, x=`dna.extraction.batch`), data=metadata.v2)+geom_boxplot()
#batch 24 looks like it might have lower DNA yield than the rest of the batches, 
#but it is still within the range of [DNA] of the other batches. 

metadata.v3$dna.extraction.batch <- as.character(metadata.v3$dna.extraction.batch)
ggplot(aes(y=`nanodrop.dna.concentration.ng.per.uL`, x=`dna.extraction.batch`), data=metadata.v3)+geom_boxplot()


#one way ANOVA
library(car)
#seeing if [DNA] meets the assumptions of ANOVA. 
#homogeneity of variances
metadata.v3$dna.extraction.batch <- as.character(metadata.v3$dna.extraction.batch)
leveneTest(nanodrop.dna.concentration.ng.per.uL ~dna.extraction.batch, data=metadata.v3)
#there is not homogeneity of variance so [DNA] needs to be transformed.
bestNormalize::bestNormalize(metadata.v3$nanodrop.dna.concentration.ng.per.uL)
#Yeo-Johnson was best but box cox was second best and more common
metadata.v3$box.cox.transformed.nanodrop.dna.concentration <- (bestNormalize::boxcox(metadata.v3$nanodrop.dna.concentration.ng.per.uL))$x.t
#now test the assumptions of ANOVA again. 
leveneTest(box.cox.transformed.nanodrop.dna.concentration ~dna.extraction.batch, data=metadata.v3)
#good on homogeneity of variance. 
#looking at distribution of residuals
model<-aov(box.cox.transformed.nanodrop.dna.concentration~dna.extraction.batch, data=metadata.v3)
summary(model)
#no difference in by DNA extraction batch. 
pairwise.t.test(x=metadata.v3$box.cox.transformed.nanodrop.dna.concentration, g=metadata.v3$dna.extraction.batch, p.adjust.method="fdr")
#no significant difference between dna extraction batches according to pairwise t-test 

#plot [DNA] by bead basher with batches #23 and #24
ggplot(aes(y=`nanodrop.dna.concentration.ng.per.uL`, x=`bead.basher.used`), data=metadata.v2)+geom_boxplot()
#no difference by bead basher

#plot [DNA] by bead basher without batches #23 and #24
ggplot(aes(y=`nanodrop.dna.concentration.ng.per.uL`, x=`bead.basher.used`), data=metadata.v3)+geom_boxplot()
#no difference by bead basher

#look specifically at the samples that were designated bead basher controls (across all batches)
replicates <- metadata[metadata$replicate=="yes",]
replicates <- replicates[-c(grep(replicates$extraction.id, pattern="330L|278W")),]
replicates$duplicate.id <- replicates$extraction.id
replicates$duplicate.id <- gsub(pattern="b$", replacement="", replicates$duplicate.id)

ggplot(aes(y=`nanodrop.dna.concentration.ng.per.uL`, x=`bead.basher.used`, group=duplicate.id), data=replicates)+geom_point()+geom_line()
replicate.oldbasher <- replicates[replicates$bead.basher.used=="old",]["nanodrop.dna.concentration.ng.per.uL"]
replicate.oldbasher <- replicate.oldbasher$nanodrop.dna.concentration.ng.per.uL

replicate.newbasher <- replicates[replicates$bead.basher.used=="new",]["nanodrop.dna.concentration.ng.per.uL"]
replicate.newbasher <- replicate.newbasher$nanodrop.dna.concentration.ng.per.uL

wilcox.test(x=replicate.newbasher, replicate.oldbasher, paired=T) #not significant difference. 

#STEP 2: Determine if read depth varied by: (a) DNA extraction batch, (b) bead basher (old vs. new) for the replicates,
#(c) plate batch, or (d) sequencing batch

#read depth by DNA extraction batch (including all batches)
ggplot(aes(y=`read_depth_afterDADA2`, x=`dna.extraction.batch`), data=metadata.v2)+geom_boxplot()

leveneTest(read_depth_afterDADA2 ~ dna.extraction.batch, data=metadata.v2)
#non-homogeneity of variances 
bestNormalize::bestNormalize(metadata.v2$read_depth_afterDADA2)
#box cox was best normalization
metadata.v2$box.cox.transformed.read_depth_afterDADA2 <- (bestNormalize::boxcox(metadata.v2$read_depth_afterDADA2))$x.t
leveneTest(box.cox.transformed.read_depth_afterDADA2 ~ dna.extraction.batch, data=metadata.v2)
#still non-homogeneity of variances so using non-parametric test.
pairwise.wilcox.test(x=metadata.v2$read_depth_afterDADA2, g=metadata.v2$dna.extraction.batch, p.adjust.method="fdr")
#quite a few of the DNA extraction batches tended to have different read depth from each other.  

#plot read depth by bead basher for the replicates
ggplot(aes(y=`read_depth_afterDADA2`, x=`bead.basher.used`, group=duplicate.id), data=replicates)+geom_point()+geom_line()

replicate.oldbasher <- replicates[replicates$bead.basher.used=="old",]["read_depth_afterDADA2"]
replicate.oldbasher <- replicate.oldbasher$read_afdada2

replicate.newbasher <- replicates[replicates$bead.basher.used=="new",]["read_depth_afterDADA2"]
replicate.newbasher <- replicate.newbasher$read.depth.after.dada2

t.test(x=replicates[replicates$bead.basher.used=="old","read_depth_afterDADA2"], y=replicates[replicates$bead.basher.used=="new", "read_depth_afterDADA2"], paired=T) #not a significant difference

#plot read depth by plate batch
unique(metadata.v2$plate.batch)
metadata.v2$plate.batch <- as.factor(metadata.v2$plate.batch)
ggplot(aes(y=`read_depth_afterDADA2`, x=`plate.batch`), data=metadata.v2)+geom_boxplot()

leveneTest(read_depth_afterDADA2 ~ plate.batch, data=metadata.v2)
#non-homogeneity of variances 
bestNormalize::bestNormalize(metadata.v2$read_depth_afterDADA2)
#box cox was best normalization
metadata.v2$box.cox.transformed.read_depth_afterDADA2 <- (bestNormalize::boxcox(metadata.v2$read_depth_afterDADA2))$x.t
leveneTest(box.cox.transformed.read_depth_afterDADA2 ~ plate.batch, data=metadata.v2)
#still non-homogeneity of variances

#plate batches 2 and 4 are significantly different from the rest. 
pairwise.wilcox.test(x=metadata.v2$read_depth_afterDADA2, g=metadata.v2$plate.batch, p.adjust.method="fdr")
#quite a few of the plates tended to have different read depth from each other.  

#plot read depth by sequencing batch
metadata.v2$sequencing.batch <- as.factor(metadata.v2$sequencing.batch)
ggplot(aes(y=`read_depth_afterDADA2`, x=`sequencing.batch`), data=metadata.v2)+geom_boxplot()
wilcox.test(x=(metadata.v2[metadata.v2$sequencing.batch=="1","read_depth_afterDADA2"]), y=(metadata.v2[metadata.v2$sequencing.batch=="2","read_depth_afterDADA2"]))
#no difference. 

#STEP 3: Determine if beta-diversity varied by DNA extraction batch, [DNA], sequencing
#batch, PCR/plate batch, or read depth. 

library(qiime2R)
#read in unweighted unifrac distance
unweighted.unifrac.dis.matrix <- read_qza("../diversity/diversity_core_metrics_samples_327/unweighted_unifrac_distance_matrix.qza")
unweighted.unifrac.dis.matrix <- unweighted.unifrac.dis.matrix$data

#select the samples that are in the distance matrix. 
metadata.v3 <- metadata[c(which(metadata$sample.id %in% attributes(unweighted.unifrac.dis.matrix)$Labels)),]
str(metadata.v3)
metadata.v3$dna.extraction.batch <- as.character(metadata.v3$dna.extraction.batch)
metadata.v3$sequencing.batch <- as.character(metadata.v3$sequencing.batch)
metadata.v3$plate.batch <- as.character(metadata.v3$plate.batch)

#generate a child ID variable from the sample.id variable
metadata.v3$cid <- gsub(pattern="[.][234][a-z]*", metadata.v3$extraction.id, replacement = "")
table(metadata.v3$cid)

#determine if unweighted unifrac varied with any of the batches. 
set.seed(123)
adonis.model.unwunifrac.extractbatch <- adonis2(formula=unweighted.unifrac.dis.matrix~cid+dna.extraction.batch, data=metadata.v3)
#beta-diversity did not vary with extraction batch. 
set.seed(123)
adonis.model.unwunifrac.platebatch <- adonis2(formula=unweighted.unifrac.dis.matrix~cid+plate.batch, data=metadata.v3)
#not different with plate/PCR batch. 
set.seed(123)
adonis.model.unwunifrac.sequencingbatch <- adonis2(formula=unweighted.unifrac.dis.matrix~cid+sequencing.batch, data=metadata.v3)
#not different with sequencing batch

#looking at weighted unifrac
weighted.unifrac.dis.matrix <- read_qza("../diversity/diversity_core_metrics_samples_327/weighted_unifrac_distance_matrix.qza")
weighted.unifrac.dis.matrix <- weighted.unifrac.dis.matrix$data

#determine if weighted unifrac varied with any of the batches, [DNA] read depth. 
set.seed(123)
adonis.model.wunifrac.extractbatch <- adonis2(formula=weighted.unifrac.dis.matrix~cid+dna.extraction.batch, data=metadata.v3)
#not different with DNA extraction batch.
set.seed(123)
adonis.model.wunifrac.platebatch <- adonis2(formula=weighted.unifrac.dis.matrix~cid+plate.batch, data=metadata.v3)
#not different with plate/PCR batch. 
set.seed(123)
adonis.model.wunifrac.sequencingbatch <- adonis2(formula=weighted.unifrac.dis.matrix~cid+sequencing.batch, data=metadata.v3)
#not different with sequencing batch

#beta-diversity did not associate with these batches, therefore we won't include
#those co-variates in the subsequent models. 

