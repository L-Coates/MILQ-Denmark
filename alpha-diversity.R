#This script describes the steps taken and the code used to determine if 
#alpha-diversity differed with diarrhea.  

#load needed libraries
library(qiime2R)
library(tidyverse)
library(phyloseq)
library(vegan)
library(lme4)
library(nlme)
library(lmerTest)
library(ggplot2)
library(ggpubr)
library(ANCOMBC)

#STEP 1: reading in the subset of the metadata file that only has the 327 Denmark stool samples
#that will be used for analysis, and also the morbidity and medication history
#and combine with parity, infant gender, and birth mode information. 
metadata <- read.csv("../morbidity_medication_prevalences/infant.morbidities.medications.csv", header = T)
metadata <- metadata[,-1]

gender.mod.parity <- read.csv("../sociodemographics_and_sample_collection_info/birthmode.gender.parity.csv", header = T)
gender.mod.parity <- gender.mod.parity[,c(grep(pattern="^mid|visit$|parity|infant[.]gender|mode[.]of[.]delivery|household|first", colnames(gender.mod.parity), ignore.case = T))]
gender.mod.parity <- distinct(gender.mod.parity)

eBF <- read.csv("../complementary_feeding_vs_Granulicatella/Exclusive.breastfeeding.before.4months.csv")

metadata.v2 <- merge(x=metadata, y=gender.mod.parity, by=c("mid", "visit"), all=F)
metadata.v2 <- merge(x=metadata.v2, y=eBF, by="mid", all=F)

#STEP 2: read in the alpha diversity metrics 
shannon <- read_qza("diversity_core_metrics_samples_327/shannon_vector.qza")
shannon <- as.data.frame(shannon$data)
shannon$sample.id <- rownames(shannon)
faith <- read_qza("diversity_core_metrics_samples_327/faith_pd_vector.qza")
faith <- as.data.frame(faith$data)
colnames(faith) <- c("sample.id", "faith.diversity")
observed.otu <- read_qza("diversity_core_metrics_samples_327/observed_features_vector.qza")
observed.otu <- observed.otu$data
observed.otu$sample.id <- rownames(observed.otu)
evenness <- read_qza("diversity_core_metrics_samples_327/evenness_vector.qza")
evenness <- evenness$data
evenness$sample.id <- rownames(evenness)

#STEP 3: combine the alpha diversity measures and the sample id's. 
metadata.v3 <- merge(x=metadata.v2, y=shannon, by="sample.id", all=F)
metadata.v3 <- merge(x=metadata.v3, y=faith, by="sample.id", all=F)
metadata.v3 <- merge(x=metadata.v3, y=observed.otu, by="sample.id", all=F)
metadata.v3 <- merge(x=metadata.v3, y=evenness, by="sample.id", all=F)

str(metadata.v3)
metadata.v3$dna.extraction.batch <- as.character(metadata.v3$dna.extraction.batch)
metadata.v3$sequencing.batch <- as.character(metadata.v3$sequencing.batch)
metadata.v3$plate.batch <- as.character(metadata.v3$plate.batch)

#STEP 4: run (generalized) linear mixed effects modeling on alpha diversity in 
#relation to the fixed effects of interest and setting individual as random effect. 

#full models 
set.seed(123)
shannon.md1 <-summary(lmer(shannon_entropy~stool_age_days+infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months + infant.diarrhea+infant.fever+infant.vomit+(1|mid), data=metadata.v3))
print(shannon.md1)
#shannon alpha diversity was positively associated with age, negatively associated with antibiotic use, and negatively associated with diarrhea
set.seed(123)
shannon.md2.diarrhea <-summary(lmer(shannon_entropy~stool_age_days+infant.antibiotic + infant.diarrhea+(1|mid), data=metadata.v3))
print(shannon.md2.diarrhea)
#diarrhea and age still associated with stool alpha diversity

set.seed(123)
shannon.md2.fever <-summary(lmer(shannon_entropy~stool_age_days+infant.antibiotic + infant.diarrhea+infant.fever+(1|mid), data=metadata.v3))
print(shannon.md2.fever)
#fever not associated with shannon diversity

set.seed(123)
shannon.md2.vomit <-summary(lmer(shannon_entropy~stool_age_days+infant.antibiotic + infant.diarrhea+infant.vomit+(1|mid), data=metadata.v3))
print(shannon.md2.vomit)
#vomit not associated with shannon diversity

set.seed(123)
faith.md1<-summary(lmer(faith.diversity~stool_age_days+infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months + infant.diarrhea+infant.fever+infant.vomit+(1|mid), data=metadata.v3))
print(faith.md1)
#faith phylogenetic diversity was positively associated with age, positively with girls,
#and positively with exclusive breastfeeding before 4 months of age.

set.seed(123)
faith.md2.diarrhea <-summary(lmer(faith.diversity~stool_age_days+infant.gender+ExclusiveBreastfeeding.before.4months+infant.diarrhea+(1|mid), data=metadata.v3))
print(faith.md2.diarrhea)
#no association with diarrhea

set.seed(123)
faith.md2.fever<-summary(lmer(faith.diversity~stool_age_days+infant.gender+ExclusiveBreastfeeding.before.4months+infant.fever+(1|mid), data=metadata.v3))
print(faith.md2.fever)
#fever not associated with faith PD.

set.seed(123)
faith.md2.vomit <-summary(lmer(faith.diversity~stool_age_days+infant.gender+ExclusiveBreastfeeding.before.4months+infant.vomit+(1|mid), data=metadata.v3))
print(faith.md2.vomit)
#vomit not associated with faith PD. 

set.seed(123)
observed_features.md1 <-summary(glmer(observed_features~stool_age_days+infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months + infant.diarrhea+infant.fever+infant.vomit+(1|mid), family=poisson, data=metadata.v3))
print(observed_features.md1)
#there was a failed model convergence so trying a transformation of the observed features
#and then running linear mixed effects modeling. 
library(bestNormalize)
bestNormalize(metadata.v3$observed_features)
metadata.v3$observed_features.boxcox.transformed <- boxcox(metadata.v3$observed_features, standardize=T)$x.t

set.seed(123)
observed_features.md1<-summary(lmer(observed_features.boxcox.transformed~stool_age_days+infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months + infant.diarrhea+infant.fever+infant.vomit+(1|mid), data=metadata.v3))
print(observed_features.md1)
#number of observed features was positively associated with age, negatively associated 
#with antibiotic use, and almost negatively associated with diarrhea (p-value = 0.062)

set.seed(123)
observed_features.md2.diarrhea <-summary(lmer(observed_features~stool_age_days+infant.antibiotic+ infant.diarrhea+(1|mid), data=metadata.v3))
print(observed_features.md2.diarrhea)
#observed OTUs were negatively associated with diarrhea and positively with age

set.seed(123)
observed_features.md2.fever <-summary(lmer(observed_features~stool_age_days+infant.antibiotic+ infant.diarrhea+infant.fever+(1|mid), data=metadata.v3))
print(observed_features.md2.fever)
#not associated with fever. 

set.seed(123)
observed_features.md2.vomit <-summary(lmer(observed_features~stool_age_days+infant.antibiotic+ infant.diarrhea+infant.vomit+(1|mid), data=metadata.v3))
print(observed_features.md2.vomit)
#not associated with vomit

#evenness
hist(metadata.v3$pielou_evenness)
set.seed(123)
evenness.md1<-summary(lmer(pielou_evenness~stool_age_days+infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months + infant.diarrhea+infant.fever+infant.vomit+(1|mid), data=metadata.v3))
print(evenness.md1)
#evenness was positively associated with age, positively with household size ( > 3),
#and negatively associated with diarrhea. 

set.seed(123)
evenness.md2.diarrhea <-summary(lmer(pielou_evenness~stool_age_days+household.size + infant.diarrhea+(1|mid), data=metadata.v3))
print(evenness.md2.diarrhea)
#diarrhea negatively associated with evenness, household size and age positively associated with evenness.

set.seed(123)
evenness.md2.fever <-summary(lmer(pielou_evenness~stool_age_days+household.size + infant.diarrhea+infant.fever+(1|mid), data=metadata.v3))
print(evenness.md2.fever)
#fever not associated with evenness.

set.seed(123)
evenness.md2.vomit <-summary(lmer(pielou_evenness~stool_age_days+household.size + infant.diarrhea+infant.vomit+(1|mid), data=metadata.v3))
print(evenness.md2.vomit)
#vomit not associated with evenness.


#STEP 5: make graphs of alpha diversity with age, antibiotic use, diarrhea, etc.
age.vs.shannon <-ggplot(aes(x=stool_age_days, y=shannon_entropy, group=mid), data=metadata.v3)+geom_point(alpha=0.3)+geom_line()+xlab("")+ylab("Shannon entropy")
age.vs.faith <-ggplot(aes(x=stool_age_days, y=faith.diversity, group=mid), data=metadata.v3)+geom_point(alpha=0.3)+geom_line()+xlab("")+ylab("Faith's phylogenetic diversity")
age.vs.observedASV <-ggplot(aes(x=stool_age_days, y=observed_features.boxcox.transformed, group=mid), data=metadata.v3)+geom_point(alpha=0.3)+geom_line()+xlab("age at stool collection (days)")+ylab("observed ASVs, Box Cox transformed")
age.vs.evenness <-ggplot(aes(x=stool_age_days, y=pielou_evenness, group=mid), data=metadata.v3)+geom_point(alpha=0.3)+geom_line()+xlab("age at stool collection (days)")+ylab("Pielou's evenness")
ggarrange(age.vs.shannon, age.vs.faith, age.vs.observedASV,age.vs.evenness, ncol=2, nrow=2, labels="AUTO")
ggsave("age.vs.alphadiversity.jpeg", plot=last_plot(), height=10, width=13, dpi=600)

#try another figure for age vs. alpha diversity that includes the regression line
#and confidence interval. 
shannon.full.model <-lmer(shannon_entropy~stool_age_days+infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months + infant.diarrhea+infant.fever+infant.vomit+(1|mid), data=metadata.v3)
faith.full.model <- lmer(faith.diversity~stool_age_days+infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months + infant.diarrhea+infant.fever+infant.vomit+(1|mid), data=metadata.v3)
observed_features.full.model <- lmer(observed_features.boxcox.transformed~stool_age_days+infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months + infant.diarrhea+infant.fever+infant.vomit+(1|mid), data=metadata.v3)
evenness.full.model <- lmer(pielou_evenness~stool_age_days+infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months + infant.diarrhea+infant.fever+infant.vomit+(1|mid), data=metadata.v3)

age.vs.shannon.RegLine <- predict_response(shannon.full.model, terms = "stool_age_days", ci_level = 0.95) |> plot(show_data=T, use_theme=T, show_title=F)+
    theme_bw()+
    xlab("")+
    ylab("Shannon entropy")

age.vs.faith.RegLine <- predict_response(faith.full.model, terms = "stool_age_days", ci_level = 0.95) |> plot(show_data=T, use_theme=T, show_title=F)+
    theme_bw()+
    xlab("")+
    ylab("Faith's phylogenetic diversity")

age.vs.observedOTU.RegLine <- predict_response(observed_features.full.model, terms = "stool_age_days", ci_level = 0.95) |> plot(show_data=T, use_theme=T, show_title=F)+
    theme_bw()+
    xlab("")+
    ylab("Observed ASVs, Box Cox transformed")

age.vs.evenness.RegLine <- predict_response(evenness.full.model, terms = "stool_age_days", ci_level = 0.95) |> plot(show_data=T, use_theme=T, show_title=F)+
    theme_bw()+
    xlab("")+
    ylab("Pielou's evenness")




table(metadata.v3$infant.antibiotic) #10 stool samples were collected during a visit when an infant had antibiotic.
antibiotics.vs.shannon<- ggplot(aes(x=infant.antibiotic, y=shannon_entropy), data=metadata.v3)+geom_boxplot()+xlab("antibiotic use")+ ylab("Shannon entropy")+scale_x_discrete(labels=c("no\nn=317", "yes\nn=10"))
antibiotics.vs.observedASV <- ggplot(aes(x=infant.antibiotic, y=observed_features.boxcox.transformed), data=metadata.v3)+geom_boxplot()+xlab("antibiotic use") + ylab("observed ASVs, Box Cox transformed")+scale_x_discrete(labels=c("no\nn=317", "yes\nn=10")) 
ggarrange(antibiotics.vs.shannon, antibiotics.vs.observedASV, ncol=2, nrow=1, labels="AUTO")
ggsave("infant.antibiotic.use.vs.alphadiversity.jpeg", plot=last_plot(), dpi=600, height=5, width=8)

faith.gender <- ggplot(aes(x=infant.gender, y=faith.diversity), data=metadata.v3)+geom_boxplot()+xlab("infant gender")+ ylab("Faith's phylogenetic diversity")+scale_x_discrete(labels=c("boy\nn=47", "girl\nn=62"))
faith.breastfeeding <- ggplot(aes(x=ExclusiveBreastfeeding.before.4months, y=faith.diversity), data=metadata.v3)+geom_boxplot()+xlab("exclusive breastfeeding the first 4 months of life")+ ylab("")+scale_x_discrete(labels=c("no\nn=29","yes\nn=80"))
ggarrange(faith.gender, faith.breastfeeding, ncol=2, nrow=1, labels="AUTO")
ggsave("infant.gender.breastfeeding.vs.faithPD.jpeg", plot=last_plot(), dpi=600, height=5, width=8)

ggplot(aes(x=household.size, y=pielou_evenness), data=metadata.v3)+geom_boxplot()+xlab("household size")+ ylab("Pielou's evenness")+scale_x_discrete(labels=c("<= 3\nn=72", "> 3\nn=37"))
ggsave("household.size.vs.evenness.jpeg", plot=last_plot(), dpi=600, height=4, width=4)

diarrhea.vs.shannon <-ggplot(aes(x=infant.diarrhea, y=shannon_entropy, color=as.factor(infant.diarrhea)), data=metadata.v3)+
    geom_boxplot()+
    theme_bw()+
    xlab("diarrhea")+
    ylab("Shannon entropy") +
    scale_x_discrete(labels=c("no\nn=287", "yes\nn=40"))+
    scale_color_manual(values=c(no="black", yes="#F8766D"))+
    theme(legend.position="none")

diarrhea.vs.observedASV <-ggplot(aes(x=infant.diarrhea, y=observed_features.boxcox.transformed, color=as.factor(infant.diarrhea)), data=metadata.v3)+
    geom_boxplot()+
    theme_bw()+
    xlab("diarrhea")+
    ylab("observed ASVs, Box Cox transformed") +
    scale_x_discrete(labels=c("no\nn=287", "yes\nn=40")) +
    scale_color_manual(values=c(no="black", yes="#F8766D"))+
    theme(legend.position="none")

diarrhea.vs.evenness <-ggplot(aes(x=infant.diarrhea, y=pielou_evenness, color=as.factor(infant.diarrhea)), data=metadata.v3)+
    geom_boxplot()+
    theme_bw()+
    xlab("diarrhea")+
    ylab("Peilou's evenness")+
    scale_x_discrete(labels=c("no\nn=287", "yes\nn=40"))+
    scale_color_manual(values=c(no="black", yes="#F8766D"))+
    theme(legend.position="none")

ggarrange(diarrhea.vs.shannon, diarrhea.vs.evenness,diarrhea.vs.observedASV, ncol=3, nrow=1, labels="AUTO")
ggsave("diarrhea.vs.alphadiversity.tiff", device="tiff", plot=last_plot(), height=5, width=10, dpi=600)



