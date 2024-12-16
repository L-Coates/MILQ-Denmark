#This script describes the steps taken and the code used to determine if beta-diversity
#differed by diarrhea.   

#load needed libraries
library(tidyverse)
library(phyloseq)
library(vegan)
library(qiime2R)
library(ggpubr)

#STEP 1: Calculating unweighted and weighted unifrac distances in QIIME2, as well
#as performing the ordination (PCoA) and the visualization with emperor.
#Also performing biplots for the (un)weighted unifrac PCoA plots in QIIME2 and
#need to generate a feature metadata table in which the taxonomic assignments/classifications
#are outlined so they can be annotated on the arrows in the biplot. 
feature.table <- read_qza("../feature-table-rarefied-6963.qza")
feature.table <- feature.table$data
taxonomy <- read_qza("../taxonomy-classification.qza")
taxonomy2 <- as.data.frame(taxonomy$data) %>% column_to_rownames("Feature.ID") %>% separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>% as.matrix()
#remove "Confidence" column
taxonomy2 <- taxonomy2[,-8]

#trim taxonomy to include only ASVs in the ASV table
taxonomy3 <- taxonomy2[c(rownames(taxonomy2)%in%rownames(feature.table)),]

#make the taxonomy into table that can be read by qiime2. 
taxonomy.table <- as.data.frame(taxonomy3)
taxonomy.table$featureid <- rownames(taxonomy.table)
taxonomy.table <- taxonomy.table[,c(8,1:7)]
write.table(taxonomy.table, "../metadata-feature-taxonomy.txt", sep="\t", row.names = F)


#STEP 2: reading in the subset of the metadata file that only has the 327 Denmark stool samples
#that will be used for analysis, and also the morbidity and medication history
#and combine with information on parity, infant gender, birth mode, and exclusive breasting
#during first 3.5 months. 
metadata <- read.csv("../morbidity_medication_prevalences/infant.morbidities.medications.csv", header = T)
metadata <- metadata[,-1]

gender.mod.parity <- read.csv("../sociodemographics_and_sample_collection_info/birthmode.gender.parity.csv", header = T)
gender.mod.parity <- gender.mod.parity[,c(grep(pattern="^mid|visit$|parity|infant[.]gender|mode[.]of[.]delivery|first|household", colnames(gender.mod.parity), ignore.case = T))]
gender.mod.parity <- distinct(gender.mod.parity)

eBF <- read.csv("../complementary_feeding_vs_Granulicatella/Exclusive.breastfeeding.before.4months.csv")

metadata.v2 <- merge(x=metadata, y=gender.mod.parity, by=c("mid","visit"), all=F)
metadata.v2 <- merge(x=metadata.v2, y=eBF, by="mid", all=F)


#STEP 3: plot the unweighted and weighted UniFrac distances
unweighted.unifrac <- read_qza("diversity_core_metrics_samples_327/unweighted_unifrac_pcoa_results.qza")
weighted.unifrac <- read_qza("diversity_core_metrics_samples_327/weighted_unifrac_pcoa_results.qza")

#make metadata dataframe compatible with PCoA results
metadata.v3 <- metadata.v2
metadata.v3$SampleID <- metadata.v3$sample.id
#the column "SampleID" should be moved to the first column of the data frame
metadata.v3 <- metadata.v3[,c(56,1:55)]

unweighted.unifrac$data$Vectors %>%
    select(SampleID, PC1, PC2) %>%
    right_join(metadata.v3) %>%
    ggplot(aes(x=PC1, y=PC2, color=(as.character(ExclusiveBreastfeeding.before.4months)), label=extraction.id)) +
    geom_point()+
    geom_text(hjust=0, vjust=0)+
    xlab(paste("PC1", round(unweighted.unifrac$data$ProportionExplained[1]*100, digits=1), "%", sep=" "))+
    ylab(paste("PC2", round(unweighted.unifrac$data$ProportionExplained[2]*100, digits=1), "%", sep=" "))+
    theme_q2r()

weighted.unifrac$data$Vectors %>%
    select(SampleID, PC1, PC2) %>%
    right_join(metadata.v3) %>%
    ggplot(aes(x=PC1, y=PC2, color=infant.gender, label=mid)) +
    geom_point()+
    geom_text(hjust=0, vjust=0)+
    xlab(paste("PC1", round(weighted.unifrac$data$ProportionExplained[1]*100, digits=1), "%", sep=" "))+
    ylab(paste("PC2", round(weighted.unifrac$data$ProportionExplained[2]*100, digits=1), "%", sep=" "))+
    theme_q2r()

#principal coordinates analysis of all the samples to determine if there is
#apparent separation.  

#STEP 4: test for separation by morbidity and other covariates of interest across all three visits. 
metadata.v3[,c(grep(pattern="^mid|dna[.]extraction|sequencing|plate|sample[.]type|country|bead[.]basher|infant[.]|visit$|mode", colnames(metadata.v3)))] <- lapply(metadata.v3[,c(grep(pattern="^mid|dna[.]extraction|sequencing|plate|sample[.]type|country|bead[.]basher|infant[.]|visit$|mode", colnames(metadata.v3)))],FUN=as.factor)
str(metadata.v3)

#first unweighted unifrac
unweighted.unifrac.dis.matrix <- read_qza("diversity_core_metrics_samples_327/unweighted_unifrac_distance_matrix.qza")
unweighted.unifrac.dis.matrix <- unweighted.unifrac.dis.matrix$data

#full model with morbidity by visit
set.seed(123)
adonis.model1.output.unw.unifrac<-adonis2(formula = unweighted.unifrac.dis.matrix ~ mid+stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size +ExclusiveBreastfeeding.before.4months+ infant.diarrhea+infant.fever+infant.vomit, data=metadata.v3)
#significant effect of individual only.

#trying similar model but without insignificant covariates and including diarrhea. 
set.seed(123)
adonis.model2.diarrhea.output.unw.unifrac<-adonis2(formula = unweighted.unifrac.dis.matrix ~ mid+ infant.diarrhea, data=metadata.v3)
#diarrhea not significant

#trying similar model but without insignificant covariates and including fever. 
set.seed(123)
adonis.model2.fever.output.unw.unifrac<-adonis2(formula = unweighted.unifrac.dis.matrix ~ mid+ infant.fever, data=metadata.v3)
#fever not significant

#trying similar model but without insignificant covariates and including vomit. 
set.seed(123)
adonis.model2.vomit.output.unw.unifrac<-adonis2(formula = unweighted.unifrac.dis.matrix ~ mid+ infant.vomit, data=metadata.v3)
#vomit not significant


#looking at weighted unifrac
weighted.unifrac.dis.matrix <- read_qza("diversity_core_metrics_samples_327/weighted_unifrac_distance_matrix.qza")
weighted.unifrac.dis.matrix <- weighted.unifrac.dis.matrix$data

#full model
set.seed(123)
adonis.model1.output.w.unifrac<-adonis2(formula = weighted.unifrac.dis.matrix ~ mid+stool_age_days + infant.gender + mode.of.delivery+ infant.probiotic+infant.antibiotic+infant.laxative+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months + infant.diarrhea+infant.fever+infant.vomit, data=metadata.v3)
#significant effect of individual only. 

#leave out the covariates that weren't significantly associated with weighted unifrac
#distance but including diarrhea
set.seed(123)
adonis.model2.output.diarrhea.w.unifrac<-adonis2(formula = weighted.unifrac.dis.matrix ~ mid+infant.diarrhea, data=metadata.v3)
#no relationship with diarrhea

set.seed(123)
adonis.model2.output.fever.w.unifrac<-adonis2(formula = weighted.unifrac.dis.matrix ~ mid+infant.fever, data=metadata.v3)
#fever not significant

set.seed(123)
adonis.model2.output.vomit.w.unifrac<-adonis2(formula = weighted.unifrac.dis.matrix ~ mid+infant.vomit, data=metadata.v3)
#vomit not significant

#STEP 5: Now look within each visit separately to determine if beta-diversity
#varied with morbidity. 

#reading in weighted and unweighted unifrac distances for each visit
#visit 2
weighted.unifrac.dis.matrix.visit2 <- read_qza("diversity_core_metrics_samples_327_visit2/weighted_unifrac_distance_matrix.qza")
weighted.unifrac.dis.matrix.visit2 <- weighted.unifrac.dis.matrix.visit2$data

unweighted.unifrac.dis.matrix.visit2 <- read_qza("diversity_core_metrics_samples_327_visit2/unweighted_unifrac_distance_matrix.qza")
unweighted.unifrac.dis.matrix.visit2 <- unweighted.unifrac.dis.matrix.visit2$data

#visit 3
weighted.unifrac.dis.matrix.visit3 <- read_qza("diversity_core_metrics_samples_327_visit3/weighted_unifrac_distance_matrix.qza")
weighted.unifrac.dis.matrix.visit3 <- weighted.unifrac.dis.matrix.visit3$data

unweighted.unifrac.dis.matrix.visit3 <- read_qza("diversity_core_metrics_samples_327_visit3/unweighted_unifrac_distance_matrix.qza")
unweighted.unifrac.dis.matrix.visit3 <- unweighted.unifrac.dis.matrix.visit3$data

#visit 4
weighted.unifrac.dis.matrix.visit4 <- read_qza("diversity_core_metrics_samples_327_visit4/weighted_unifrac_distance_matrix.qza")
weighted.unifrac.dis.matrix.visit4 <- weighted.unifrac.dis.matrix.visit4$data

unweighted.unifrac.dis.matrix.visit4 <- read_qza("diversity_core_metrics_samples_327_visit4/unweighted_unifrac_distance_matrix.qza")
unweighted.unifrac.dis.matrix.visit4 <- unweighted.unifrac.dis.matrix.visit4$data

#visit 2 models
metadata.visit2 <- metadata.v3[metadata.v3$visit=="2",]

#unweighted unifrac in visit 2
set.seed(123)
adonis.visit2.model1.output.unw.unifrac<-adonis2(formula = unweighted.unifrac.dis.matrix.visit2 ~ stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months + infant.diarrhea+infant.fever+infant.vomit, data=metadata.visit2)
#unweighted unifrac in visit 2 was different by antibiotic use and by infant gender.

#trying cut down model with only antibiotic use, gender, and one of the morbidities at a time. 
set.seed(123)
adonis.visit2.model2.output.diarrhea.unw.unifrac<-adonis2(formula = unweighted.unifrac.dis.matrix.visit2 ~ infant.antibiotic+infant.gender+infant.diarrhea, data=metadata.visit2)
#association with antibiotic still, but no association with diarrhea or infant gender. 

set.seed(123)
adonis.visit2.model2.output.fever.unw.unifrac<-adonis2(formula = unweighted.unifrac.dis.matrix.visit2 ~ infant.antibiotic+infant.gender+infant.fever, data=metadata.visit2)
#association with antibiotic still, but no association with fever or gender

set.seed(123)
adonis.visit2.model2.output.vomit.unw.unifrac<-adonis2(formula = unweighted.unifrac.dis.matrix.visit2 ~ infant.antibiotic+infant.gender+infant.vomit, data=metadata.visit2)
#association with antibiotic still, but no association with vomit or gender

#weighted unifrac in visit 2
set.seed(123)
adonis.visit2.model1.output.w.unifrac<-adonis2(formula = weighted.unifrac.dis.matrix.visit2 ~ stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months + infant.diarrhea+infant.fever+infant.vomit, data=metadata.visit2)
#weighted unifrac in visit 2 was only associated with infant gender.

#cut down models with one morbidity at a time
set.seed(123)
adonis.visit2.model2.output.diarrhea.w.unifrac<-adonis2(formula = weighted.unifrac.dis.matrix.visit2 ~ infant.gender+infant.diarrhea, data=metadata.visit2)
#infant gender still associated with beta-diversity, but not diarrhea

set.seed(123)
adonis.visit2.model2.output.fever.w.unifrac<-adonis2(formula = weighted.unifrac.dis.matrix.visit2 ~ infant.gender+infant.fever, data=metadata.visit2)
#infant gender still associated with beta-diversity, but not fever

set.seed(123)
adonis.visit2.model2.output.vomit.w.unifrac<-adonis2(formula = weighted.unifrac.dis.matrix.visit2 ~ infant.gender+infant.vomit, data=metadata.visit2)
#infant gender still associated with beta-diversity, but not vomit

#visit 3 models
metadata.visit3 <- metadata.v3[metadata.v3$visit=="3",]

set.seed(123)
adonis.visit3.model1.output.unw.unifrac<-adonis2(formula = unweighted.unifrac.dis.matrix.visit3 ~ stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months + infant.diarrhea+infant.fever+infant.vomit, data=metadata.visit3)
#no associations.

#trying cut down model with one of the morbidities at a time. 
set.seed(123)
adonis.visit3.model2.output.diarrhea.unw.unifrac<-adonis2(formula = unweighted.unifrac.dis.matrix.visit3 ~ infant.diarrhea, data=metadata.visit3)
#no association with diarrhea

set.seed(123)
adonis.visit3.model2.output.fever.unw.unifrac<-adonis2(formula = unweighted.unifrac.dis.matrix.visit3 ~ infant.fever, data=metadata.visit3)
#no association with fever

set.seed(123)
adonis.visit3.model2.output.vomit.unw.unifrac<-adonis2(formula = unweighted.unifrac.dis.matrix.visit3 ~ infant.vomit, data=metadata.visit3)
#no association with vomit

#now for weighted unifrac in visit 3
set.seed(123)
adonis.visit3.model1.output.w.unifrac<-adonis2(formula = weighted.unifrac.dis.matrix.visit3 ~ stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months + infant.diarrhea+infant.fever+infant.vomit, data=metadata.visit3)
#weighted unifrac in visit 3 was associated with exclusive breastfeeding before 4 months

set.seed(123)
adonis.visit3.model2.output.diarrhea.w.unifrac <-adonis2(formula = weighted.unifrac.dis.matrix.visit3 ~ ExclusiveBreastfeeding.before.4months+infant.diarrhea, data=metadata.visit3)
#but exclusive breastfeeding before 4 months isn't significant and infant diarrhea isn't either

set.seed(123)
adonis.visit3.model2.output.fever.w.unifrac<-adonis2(formula = weighted.unifrac.dis.matrix.visit3 ~ ExclusiveBreastfeeding.before.4months+infant.fever, data=metadata.visit3)
#exclusive breastfeeding before 4 months isn't significant and infant diarrhea isn't either

set.seed(123)
adonis.visit3.model2.output.vomit.w.unifrac<-adonis2(formula = weighted.unifrac.dis.matrix.visit3 ~ ExclusiveBreastfeeding.before.4months+infant.vomit, data=metadata.visit3)
#exclusive breastfeeding before 4 months isn't significant and infant diarrhea isn't either


#visit 4 models
metadata.visit4 <- metadata.v3[metadata.v3$visit=="4",]
set.seed(123)
adonis.visit4.model1.output.unw.unifrac<-adonis2(formula = unweighted.unifrac.dis.matrix.visit4 ~ stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months + infant.diarrhea+infant.fever+infant.vomit, data=metadata.visit4)
#only exclusive breastfeeding before 4 months was associated with unweighted unifrac distance

#trying cut down model with one of the morbidities at a time. 
set.seed(123)
adonis.visit4.model2.output.diarrhea.unw.unifrac<-adonis2(formula = unweighted.unifrac.dis.matrix.visit4 ~ ExclusiveBreastfeeding.before.4months+infant.diarrhea, data=metadata.visit4)
#exclusive breastfeeding before 4 months was associated with unweighted unifrac distance, but diarrhea was not. 

set.seed(123)
adonis.visit4.model2.output.fever.unw.unifrac<-adonis2(formula = unweighted.unifrac.dis.matrix.visit4 ~ ExclusiveBreastfeeding.before.4months + infant.fever, data=metadata.visit4)
#exclusive breastfeeding before 4 months was associated with unweighted unifrac distance, but diarrhea was not. 

set.seed(123)
adonis.visit4.model2.output.vomit.unw.unifrac<-adonis2(formula = unweighted.unifrac.dis.matrix.visit4 ~ ExclusiveBreastfeeding.before.4months + infant.vomit, data=metadata.visit4)
#exclusive breastfeeding before 4 months was associated with unweighted unifrac distance, but diarrhea was not. 

#now for weighted unifrac in visit 4.
set.seed(123)
adonis.visit4.model1.output.w.unifrac<-adonis2(formula = weighted.unifrac.dis.matrix.visit4 ~ stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months + infant.diarrhea+infant.fever+infant.vomit, data=metadata.visit4)
#no variables associated with weighted unifrac in visit 4. 

set.seed(123)
adonis.visit4.model2.output.diarrhea.w.unifrac<-adonis2(formula = weighted.unifrac.dis.matrix.visit4 ~ infant.diarrhea, data=metadata.visit4)
#no association

set.seed(123)
adonis.visit4.model2.output.fever.w.unifrac<-adonis2(formula = weighted.unifrac.dis.matrix.visit4 ~ infant.fever, data=metadata.visit4)
#no association

set.seed(123)
adonis.visit4.model2.output.vomit.w.unifrac<-adonis2(formula = weighted.unifrac.dis.matrix.visit4 ~ infant.vomit, data=metadata.visit4)
#no association

#STEP 6: generate PCoA of visit 2 unweighted unifrac with antibiotic coloring,
#a PCoA of visit 2 weighted unifrac with coloring for infant gender, and a PCoA 
#of visit 4 unweighted unifrac with coloring for exclusive breastfeeding in first 4 months

unweighted.unifrac.visit2 <- read_qza("diversity_core_metrics_samples_327_visit2/unweighted_unifrac_pcoa_results.qza")
weighted.unifrac.visit2 <- read_qza("diversity_core_metrics_samples_327_visit2/weighted_unifrac_pcoa_results.qza")
unweighted.unifrac.visit4 <- read_qza("diversity_core_metrics_samples_327_visit4/unweighted_unifrac_pcoa_results.qza")

pcoa.unwunifrac.visit2 <- unweighted.unifrac.visit2$data$Vectors %>%
    select(SampleID, PC1, PC2) %>%
    right_join(metadata.v3[metadata.v3$visit=="2",]) %>%
    ggplot(aes(x=PC1, y=PC2, color=infant.antibiotic)) +
    geom_point()+
    theme_bw()+
    xlab(paste("PC1", round(unweighted.unifrac.visit2$data$ProportionExplained[1]*100, digits=1), "%", sep=" "))+
    ylab(paste("PC2", round(unweighted.unifrac.visit2$data$ProportionExplained[2]*100, digits=1), "%", sep=" "))+
    theme(legend.position="top")+
    scale_color_manual(values=c("black","red"), name="unweighted UniFrac in visit 2: infant antibiotic use")

pcoa.wunifrac.visit2 <- weighted.unifrac.visit2$data$Vectors %>%
    select(SampleID, PC1, PC2) %>%
    right_join(metadata.v3[metadata.v3$visit=="2",]) %>%
    ggplot(aes(x=PC1, y=PC2, shape=infant.gender, color=infant.gender)) +
    geom_point()+
    theme_bw()+
    xlab(paste("PC1", round(weighted.unifrac.visit2$data$ProportionExplained[1]*100, digits=1), "%", sep=" "))+
    ylab(paste("PC2", round(weighted.unifrac.visit2$data$ProportionExplained[2]*100, digits=1), "%", sep=" "))+
    theme(legend.position="top")+
    scale_color_manual(name="weighted UniFrac in visit 2: infant gender", values=c("blue", "#FF33FF"), labels=c("male", "female"))+
    scale_shape_manual(name = "weighted UniFrac in visit 2: infant gender", values=c(15, 17), labels=c("male", "female"))

pcoa.unwunifrac.visit4 <- unweighted.unifrac.visit4$data$Vectors %>%
    select(SampleID, PC1, PC2) %>%
    right_join(metadata.v3[metadata.v3$visit=="4",]) %>%
    ggplot(aes(x=PC1, y=PC2, color=(as.character(ExclusiveBreastfeeding.before.4months)))) +
    geom_point()+
    theme_bw()+
    xlab(paste("PC1", round(unweighted.unifrac.visit4$data$ProportionExplained[1]*100, digits=1), "%", sep=" "))+
    ylab(paste("PC2", round(unweighted.unifrac.visit4$data$ProportionExplained[2]*100, digits=1), "%", sep=" "))+
    theme(legend.position="top")+
    scale_color_manual(values=c("purple", "green"), name = "unweighted UniFrac in visit 4: exclusive breastfeeding before 4 months of age")

ggarrange(pcoa.unwunifrac.visit2,pcoa.wunifrac.visit2, pcoa.unwunifrac.visit4, ncol=1, nrow=3, labels="AUTO")
ggsave("PCoA.unifrac.by.visit.tiff",device="tiff", plot=last_plot(), height=10, width=7, dpi=600)


#STEP 7: Determine if beta-diversity in visit 2 was associated with morbidity
#later on. 

set.seed(123)
adonis.visit2.predict.diarrhea.model1.unw.unifrac <- adonis2(formula = unweighted.unifrac.dis.matrix.visit2 ~ stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months + infant.diarrhea.in.visit3.or.visit4+infant.fever+infant.vomit, data=metadata.visit2)
#only association with infant antibiotic and infant gender 

#cut down model only retaining antibiotic use, gender and diarrhea in visit 3 or 4. 
set.seed(123)
adonis.visit2.predict.diarrhea.model2.unw.unifrac <- adonis2(formula = unweighted.unifrac.dis.matrix.visit2 ~ infant.antibiotic+infant.gender+infant.diarrhea.in.visit3.or.visit4, data=metadata.visit2)
#no association with morbidity

set.seed(123)
adonis.visit2.predict.fever.model1.unw.unifrac <- adonis2(formula = unweighted.unifrac.dis.matrix.visit2 ~ stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months + infant.fever.in.visit3.or.visit4+infant.diarrhea+infant.vomit, data=metadata.visit2)
#infant gender and antibiotic use has associations 

#cut down model only retaining gender and diarrhea in visit 3 or 4. 
set.seed(123)
adonis.visit2.predict.fever.model2.unw.unifrac <- adonis2(formula = unweighted.unifrac.dis.matrix.visit2 ~ infant.gender+infant.antibiotic+infant.fever.in.visit3.or.visit4, data=metadata.visit2)
#no association with fever.

set.seed(123)
adonis.visit2.predict.vomit.model1.unw.unifrac <- adonis2(formula = unweighted.unifrac.dis.matrix.visit2 ~ stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months + infant.vomit.in.visit3.or.visit4+infant.fever+infant.diarrhea, data=metadata.visit2)
#infant gender and antibiotic use were associated

#cut down model only retaining gender and diarrhea in visit 3 or 4. 
set.seed(123)
adonis.visit2.predict.vomit.model2.unw.unifrac <- adonis2(formula = unweighted.unifrac.dis.matrix.visit2 ~ infant.gender+infant.antibiotic+infant.vomit.in.visit3.or.visit4, data=metadata.visit2)
#no association with vomit


#weighted unifrac
set.seed(123)
adonis.visit2.predict.diarrhea.model1.w.unifrac <- adonis2(formula = weighted.unifrac.dis.matrix.visit2 ~ stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months + infant.diarrhea.in.visit3.or.visit4+infant.fever+infant.vomit, data=metadata.visit2)
#infant gender was associated

#weighted unifrac
set.seed(123)
adonis.visit2.predict.diarrhea.model2.w.unifrac <- adonis2(formula = weighted.unifrac.dis.matrix.visit2 ~ infant.gender+ infant.diarrhea.in.visit3.or.visit4, data=metadata.visit2)
#no association with diarrhea

set.seed(123)
adonis.visit2.predict.fever.model1.w.unifrac <- adonis2(formula = weighted.unifrac.dis.matrix.visit2 ~ stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months + infant.diarrhea+infant.fever.in.visit3.or.visit4+infant.vomit, data=metadata.visit2)
#infant gender was associated

#weighted unifrac
set.seed(123)
adonis.visit2.predict.fever.model2.w.unifrac <- adonis2(formula = weighted.unifrac.dis.matrix.visit2 ~ infant.gender+ infant.fever.in.visit3.or.visit4, data=metadata.visit2)
#no association with fever

set.seed(123)
adonis.visit2.predict.vomit.model1.w.unifrac <- adonis2(formula = weighted.unifrac.dis.matrix.visit2 ~ stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months + infant.diarrhea+infant.fever+infant.vomit.in.visit3.or.visit4, data=metadata.visit2)
#infant gender was associated

#weighted unifrac
set.seed(123)
adonis.visit2.predict.vomit.model2.w.unifrac <- adonis2(formula = weighted.unifrac.dis.matrix.visit2 ~ infant.gender+ infant.fever.in.visit3.or.visit4, data=metadata.visit2)
#no association with vomit
