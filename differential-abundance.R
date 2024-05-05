#This script describes the steps taken to assess differential abundance
#by diarrhea and conditions of interest using different approaches. 

#load needed libraries
library(qiime2R)
library(tidyverse)
library(phyloseq)
library(vegan)
library(ANCOMBC)
library(DT)
library(ggpubr)

#Part 1: read in the metadata, (unrarefied) feature table, and the taxonomy table and make a phyloseq object

#Step 1: read in metadata files
metadata <- read.csv("../morbidity_medication_prevalences/infant.morbidities.medications.csv", header = T)
metadata <- metadata[,-1]

gender.mod.parity <- read.csv("../sociodemographics_and_sample_collection_info/birthmode.gender.parity.csv", header = T)
gender.mod.parity <- gender.mod.parity[,c(grep(pattern="^mid|visit$|parity|infant[.]gender|mode[.]of[.]delivery|first|household", colnames(gender.mod.parity), ignore.case = T))]
gender.mod.parity <- distinct(gender.mod.parity)

eBF <- read.csv("../complementary_feeding_vs_Granulicatella/Exclusive.breastfeeding.before.4months.csv")

#Step 2: combine metadata files. 
metadata.v2 <- merge(x=metadata, y=gender.mod.parity, by=c("mid", "visit"), all=F)
metadata.v2 <- merge(x=metadata.v2, y=eBF, by="mid", all=F)
str(metadata.v2)

metadata.v2 <- metadata.v2[,c(3:55,1:2)]

#Step 3: read in the feature table. 
feature.table <- read_qza("../feature-table-no-chloroplast-eukarya-mitochondria.qza")
feature.table <- feature.table$data

#Step 4: read in the taxonomy. 
taxonomy <- read_qza("../taxonomy-classification.qza")
taxonomy2 <- as.data.frame(taxonomy$data) %>% column_to_rownames("Feature.ID") %>% separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>% as.matrix()
#remove "Confidence" column
taxonomy2 <- taxonomy2[,-8]

#trim taxonomy to include only ASVs in the ASV table
taxonomy3 <- taxonomy2[c(rownames(taxonomy2)%in%rownames(feature.table)),]

#Part 2: run ANCOM-BC2 on all three visits (i.e. the longitudinal dataset of 327 samples).

#Step 1: create phyloseq object
phylobj <- phyloseq(otu_table(feature.table, taxa_are_rows=T),tax_table(taxonomy3), sample_data(metadata.v2 %>% as.data.frame() %>% column_to_rownames("sample.id")))
#get rid of the taxa that are not present at all
phylobj.filtered <- filter_taxa(phylobj, function(x) sum(x)>0, TRUE)
#look at the number of ASVs among the feature table
dim(otu_table(phylobj.filtered))
#there are 1782 ASVs
dim(sample_data(phylobj.filtered)) #327 samples. 

#Step 2: run ANCOM-BC2 at ASV level with complete model (i.e. all covariates and morbidities)
set.seed(123)
ancom.md1.output.ASV <-ancombc2(data=phylobj.filtered, tax_level=NULL, fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md1.output.ASV.res <- ancom.md1.output.ASV$res 
lapply(ancom.md1.output.ASV.res[,c(grep(pattern="diff_", colnames(ancom.md1.output.ASV.res)))], FUN=table)
#29 ASVs were different by infant age

#Step 3: removing the co-variates that did not associate with any ASVs (i.e. only 
#keep infant age as covariate), and run a similar model with each morbidity separately. 
set.seed(123)
ancom.md2.output.diarrhea.ASV <-ancombc2(data=phylobj.filtered, tax_level=NULL, fix_formula="stool_age_days+infant.diarrhea", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.output.diarrhea.ASV.res <- ancom.md2.output.diarrhea.ASV$res 
lapply(ancom.md2.output.diarrhea.ASV.res[,c(grep(pattern="diff_", colnames(ancom.md2.output.diarrhea.ASV.res)))], FUN=table)
#25 ASVs were different by infant age and none were different by diarrhea. . 

set.seed(123)
ancom.md2.output.fever.ASV <-ancombc2(data=phylobj.filtered, tax_level=NULL, fix_formula="stool_age_days+infant.fever", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.output.fever.ASV.res <- ancom.md2.output.fever.ASV$res 
lapply(ancom.md2.output.fever.ASV.res[,c(grep(pattern="diff_", colnames(ancom.md2.output.fever.ASV.res)))], FUN=table)
#29 ASVs were different by infant age and none were different by fever. 

set.seed(123)
ancom.md2.output.vomit.ASV <-ancombc2(data=phylobj.filtered, tax_level=NULL, fix_formula="stool_age_days+infant.vomit", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.output.vomit.ASV.res <- ancom.md2.output.vomit.ASV$res 
lapply(ancom.md2.output.vomit.ASV.res[,c(grep(pattern="diff_", colnames(ancom.md2.output.vomit.ASV.res)))], FUN=table)
#31 ASVs were different by infant age and none were different by vomiting. 

#look at the species level with all co-variates
set.seed(123)
ancom.md1.output.species <-ancombc2(data=phylobj.filtered, tax_level="Species", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md1.output.species.res <- ancom.md1.output.species$res 
lapply(ancom.md1.output.species.res[,c(grep(pattern="diff_", colnames(ancom.md1.output.species.res)))], FUN=table)
#14 species differed by age and one by diarrhea -- s__uncultured_bacterium_17. 
print(ancom.md1.output.species.res[ancom.md1.output.species.res$diff_infant.diarrheayes=="TRUE",])
print(ancom.md1.output.species.res[ancom.md1.output.species.res$diff_stool_age_days=="TRUE",])


#make a figure of the log fold change of these species in relation to age in days. 
species.diff.by.age <- ancom.md1.output.species.res[ancom.md1.output.species.res$diff_stool_age_days=="TRUE",]
ggplot(aes(x=reorder(taxon,-lfc_stool_age_days), y=lfc_stool_age_days), data=species.diff.by.age)+geom_col()+ylab("infant age at stool collection (days)")+xlab("")+theme(axis.text.x = element_text(angle = 45, hjust=1))


#cut down models at species level
set.seed(123)
ancom.md2.output.diarrhea.species <-ancombc2(data=phylobj.filtered, tax_level="Species", fix_formula="stool_age_days+infant.diarrhea", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.output.diarrhea.species.res <- ancom.md2.output.diarrhea.species$res
lapply(ancom.md2.output.diarrhea.species.res[,c(grep(pattern="diff_", colnames(ancom.md2.output.diarrhea.species.res)))], FUN=table)
#two species different by diarrhea, 12 species by age. 
print(ancom.md2.output.diarrhea.species.res[ancom.md2.output.diarrhea.species.res$diff_infant.diarrheayes=="TRUE",])
#a species of Granulicatella and an uncultured bacterium were both positively
#associated with diarrhea. 

set.seed(123)
ancom.md2.output.fever.species <-ancombc2(data=phylobj.filtered, tax_level="Species", fix_formula="stool_age_days+infant.fever", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.output.fever.species.res <- ancom.md2.output.fever.species$res
lapply(ancom.md2.output.fever.species.res[,c(grep(pattern="diff_", colnames(ancom.md2.output.fever.species.res)))], FUN=table)
#14 species different by age, none with fever. 

set.seed(123)
ancom.md2.output.vomit.species <-ancombc2(data=phylobj.filtered, tax_level="Species", fix_formula="stool_age_days+infant.vomit", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.output.vomit.species.res <- ancom.md2.output.vomit.species$res
lapply(ancom.md2.output.vomit.species.res[,c(grep(pattern="diff_", colnames(ancom.md2.output.vomit.species.res)))], FUN=table)
#16 species different by age, none by vomit. 

#look at the genus level
set.seed(123)
ancom.md1.output.genus <-ancombc2(data=phylobj.filtered, tax_level="Genus", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md1.output.genus.res <- ancom.md1.output.genus$res 
lapply(ancom.md1.output.genus.res[,c(grep(pattern="diff_", colnames(ancom.md1.output.genus.res)))], FUN=table)
#15 genera were different by age, and one genus -- Granulicatella -- was different with diarrhea.
print(ancom.md1.output.genus.res[ancom.md1.output.genus.res$diff_infant.diarrheayes=="TRUE",])
print(ancom.md1.output.genus.res[ancom.md1.output.genus.res$diff_stool_age_days=="TRUE",])

#make a figure of the log fold change of these genera in relation to age in days. 
genus.diff.by.age <- ancom.md1.output.genus.res[ancom.md1.output.genus.res$diff_stool_age_days=="TRUE",]
genus.diff.by.age$taxon <- gsub(pattern="Genus: g__", replacement="", genus.diff.by.age$taxon) 
#pull out the taxonomy for just the 15 genera different by age. 
taxonomy4 <- data.frame(taxonomy3)
rownames(taxonomy4)=NULL
taxonomy4 <- distinct(taxonomy4)
taxonomy4$Genus <- gsub(pattern=" g__", replacement="", taxonomy4$Genus)
genus.diff.by.age.taxonomy <- taxonomy4[c(which(taxonomy4$Genus %in% genus.diff.by.age$taxon)),]
ggplot(aes(x=reorder(taxon,-lfc_stool_age_days), y=lfc_stool_age_days, fill=lfc_stool_age_days), data=genus.diff.by.age)+geom_col()+geom_errorbar(aes(ymin = lfc_stool_age_days - se_stool_age_days, ymax = lfc_stool_age_days + se_stool_age_days), width = 0.2, position = position_dodge(0.05), color = "black")+ylab("log(fold change) of normalized genus counts with infant age (days)")+xlab("")+theme(axis.text.x = element_text(angle = 45, hjust=1))+scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0)+labs(fill="log(fold change)")
ggsave("stool.genus.logFC.by.age.jpeg", dpi=600, plot=last_plot(), width=10, height=6.5)

#cut down model to remove covariates that were not significantly associated with a genus
set.seed(123)
ancom.md2.output.diarrhea.genus <-ancombc2(data=phylobj.filtered, tax_level="Genus", fix_formula="stool_age_days+infant.diarrhea", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.output.diarrhea.genus.res <- ancom.md2.output.diarrhea.genus$res 
lapply(ancom.md2.output.diarrhea.genus.res[,c(grep(pattern="diff_", colnames(ancom.md2.output.diarrhea.genus.res)))], FUN=table)
#one genus different by diarrhea, 16 genera different by age. 
print(ancom.md2.output.diarrhea.genus.res[ancom.md2.output.diarrhea.genus.res$diff_infant.diarrheayes=="TRUE",])
#Granulicatella was different by diarrhea. 
#checking that this genus was not sensitive to pseudo count for the condition of interest (infant diarrhea)
pseudo_sens = ancom.md2.output.diarrhea.genus$pseudo_sens_tab

set.seed(123)
ancom.md2.output.fever.genus <-ancombc2(data=phylobj.filtered, tax_level="Genus", fix_formula="stool_age_days+infant.fever", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.output.fever.genus.res <- ancom.md2.output.fever.genus$res 
lapply(ancom.md2.output.fever.genus.res[,c(grep(pattern="diff_", colnames(ancom.md2.output.fever.genus.res)))], FUN=table)
#no genera were different with fever

set.seed(123)
ancom.md2.output.vomit.genus <-ancombc2(data=phylobj.filtered, tax_level="Genus", fix_formula="stool_age_days+infant.vomit", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.output.vomit.genus.res <- ancom.md2.output.vomit.genus$res 
lapply(ancom.md2.output.vomit.genus.res[,c(grep(pattern="diff_", colnames(ancom.md2.output.vomit.genus.res)))], FUN=table)
#no genera were different with vomit. 

#look at the family level
set.seed(123)
ancom.md1.output.family <-ancombc2(data=phylobj.filtered, tax_level="Family", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md1.output.family.res <- ancom.md1.output.family$res 
lapply(ancom.md1.output.family.res[,c(grep(pattern="diff_", colnames(ancom.md1.output.family.res)))], FUN=table)
print(ancom.md1.output.family.res[ancom.md1.output.family.res$diff_infant.diarrheayes=="TRUE",])
print(ancom.md1.output.family.res[ancom.md1.output.family.res$diff_stool_age_days=="TRUE",])
#11 families were different by age, and one family -- Carnobacteriaceae -- was different with diarrhea (positively).

#look at the family level with the cut-down model leaving out the co-variates
#that didn't associate with any family-level microbial taxa
set.seed(123)
ancom.md2.output.diarrhea.family <-ancombc2(data=phylobj.filtered, tax_level="Family", fix_formula="stool_age_days + infant.diarrhea", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.output.diarrhea.family.res <- ancom.md2.output.diarrhea.family$res 
lapply(ancom.md2.output.diarrhea.family.res[,c(grep(pattern="diff_", colnames(ancom.md2.output.diarrhea.family.res)))], FUN=table)
print(ancom.md2.output.diarrhea.family.res[ancom.md2.output.diarrhea.family.res$diff_infant.diarrheayes=="TRUE",])
print(ancom.md2.output.diarrhea.family.res[ancom.md2.output.diarrhea.family.res$diff_stool_age_days=="TRUE",])
#14 families were different by age, and one family -- Carnobacteriaceae -- was different with diarrhea.

set.seed(123)
ancom.md2.output.fever.family <-ancombc2(data=phylobj.filtered, tax_level="Family", fix_formula="stool_age_days+infant.fever", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.output.fever.family.res <- ancom.md2.output.fever.family$res 
lapply(ancom.md2.output.fever.family.res[,c(grep(pattern="diff_", colnames(ancom.md2.output.fever.family.res)))], FUN=table)
print(ancom.md2.output.fever.family.res[ancom.md2.output.fever.family.res$diff_infant.feveryes=="TRUE",])
print(ancom.md2.output.fever.family.res[ancom.md2.output.fever.family.res$diff_stool_age_days=="TRUE",])
#none differed with fever

set.seed(123)
ancom.md2.output.vomit.family <-ancombc2(data=phylobj.filtered, tax_level="Family", fix_formula="stool_age_days+infant.vomit", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.output.vomit.family.res <- ancom.md2.output.vomit.family$res 
lapply(ancom.md2.output.vomit.family.res[,c(grep(pattern="diff_", colnames(ancom.md2.output.vomit.family.res)))], FUN=table)
print(ancom.md2.output.vomit.family.res[ancom.md2.output.vomit.family.res$diff_infant.vomityes=="TRUE",])
print(ancom.md2.output.vomit.family.res[ancom.md2.output.vomit.family.res$diff_stool_age_days=="TRUE",])
#none differed with vomit

#look at the order level
set.seed(123)
ancom.md1.output.order <-ancombc2(data=phylobj.filtered, tax_level="Order", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md1.output.order.res <- ancom.md1.output.order$res 
lapply(ancom.md1.output.order.res[,c(grep(pattern="diff_", colnames(ancom.md1.output.order.res)))], FUN=table)
print(ancom.md1.output.order.res[ancom.md1.output.order.res$diff_infant.diarrheayes=="TRUE",])
print(ancom.md1.output.order.res[ancom.md1.output.order.res$diff_stool_age_days=="TRUE",])
#9 orders were different by age, and none were different with diarrhea or other co-variates


#look at the order level with the cut-down model leaving out the co-variates
#that didn't associate with any order-level microbial taxa
set.seed(123)
ancom.md2.output.diarrhea.order <-ancombc2(data=phylobj.filtered, tax_level="Order", fix_formula="stool_age_days + infant.diarrhea", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.output.diarrhea.order.res <- ancom.md2.output.diarrhea.order$res 
lapply(ancom.md2.output.diarrhea.order.res[,c(grep(pattern="diff_", colnames(ancom.md2.output.diarrhea.order.res)))], FUN=table)
print(ancom.md2.output.diarrhea.order.res[ancom.md2.output.diarrhea.order.res$diff_infant.diarrheayes=="TRUE",])
print(ancom.md2.output.diarrhea.order.res[ancom.md2.output.diarrhea.order.res$diff_stool_age_days=="TRUE",])
#9 orders were different with age, and none were different with diarrhea

set.seed(123)
ancom.md2.output.fever.order <-ancombc2(data=phylobj.filtered, tax_level="Order", fix_formula="stool_age_days + infant.fever", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.output.fever.order.res <- ancom.md2.output.fever.order$res 
lapply(ancom.md2.output.fever.order.res[,c(grep(pattern="diff_", colnames(ancom.md2.output.fever.order.res)))], FUN=table)
print(ancom.md2.output.fever.order.res[ancom.md2.output.fever.order.res$diff_infant.feverayes=="TRUE",])
print(ancom.md2.output.fever.order.res[ancom.md2.output.fever.order.res$diff_stool_age_days=="TRUE",])
#9 orders different with age, none with fever.

set.seed(123)
ancom.md2.output.vomit.order <-ancombc2(data=phylobj.filtered, tax_level="Order", fix_formula="stool_age_days + infant.vomit", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.output.vomit.order.res <- ancom.md2.output.vomit.order$res 
lapply(ancom.md2.output.vomit.order.res[,c(grep(pattern="diff_", colnames(ancom.md2.output.vomit.order.res)))], FUN=table)
print(ancom.md2.output.vomit.order.res[ancom.md2.output.vomit.order.res$diff_infant.vomityes=="TRUE",])
print(ancom.md2.output.vomit.order.res[ancom.md2.output.vomit.order.res$diff_stool_age_days=="TRUE",])
#10 orders were different with age, and 1 with mode of birth, and none with vomiting. 


#Part 3: Within each visit, determine if there are associations between stool microbes
#and morbidity. 

#Step 1: make a phyloseq object for each visit set. 
metadata.v2$visit <- as.character(metadata.v2$visit)
metadata.v2.visit2 <- metadata.v2[metadata.v2$visit =="2",]
rownames(metadata.v2.visit2) = NULL
metadata.v2.visit3 <- metadata.v2[metadata.v2$visit =="3",]
rownames(metadata.v2.visit3) = NULL
metadata.v2.visit4 <- metadata.v2[metadata.v2$visit =="4",]
rownames(metadata.v2.visit4) = NULL

feature.table.visit2 <- feature.table[,c(which(colnames(feature.table) %in% metadata.v2.visit2$sample.id))]
dim(feature.table.visit2)
feature.table.visit3 <- feature.table[,c(which(colnames(feature.table) %in% metadata.v2.visit3$sample.id))]
dim(feature.table.visit3)
feature.table.visit4 <- feature.table[,c(which(colnames(feature.table) %in% metadata.v2.visit4$sample.id))]
dim(feature.table.visit4)

phylobj.visit2 <- phyloseq(otu_table(feature.table.visit2, taxa_are_rows=T),tax_table(taxonomy3), sample_data((metadata.v2.visit2) %>% as.data.frame() %>% column_to_rownames("sample.id")))
phylobj.visit3 <- phyloseq(otu_table(feature.table.visit3, taxa_are_rows=T),tax_table(taxonomy3), sample_data((metadata.v2.visit3) %>% as.data.frame() %>% column_to_rownames("sample.id")))
phylobj.visit4 <- phyloseq(otu_table(feature.table.visit4, taxa_are_rows=T),tax_table(taxonomy3), sample_data((metadata.v2.visit4) %>% as.data.frame() %>% column_to_rownames("sample.id")))

#get rid of the taxa that are not present at all in the visit
phylobj.visit2.filtered <- filter_taxa(phylobj.visit2, function(x) sum(x)>0, TRUE)
phylobj.visit3.filtered <- filter_taxa(phylobj.visit3, function(x) sum(x)>0, TRUE)
phylobj.visit4.filtered <- filter_taxa(phylobj.visit4, function(x) sum(x)>0, TRUE)

#look at the number of ASVs among the feature table for each visit
dim(otu_table(phylobj.visit2.filtered))
#794 ASVs
dim(otu_table(phylobj.visit3.filtered))
#803 ASVs
dim(otu_table(phylobj.visit4.filtered))
#1131 ASVs

#confirm there are 109 samples in for each visit
dim(sample_data(phylobj.visit2.filtered))
dim(sample_data(phylobj.visit3.filtered)) 
dim(sample_data(phylobj.visit4.filtered)) 


#Step 2: determine if any taxa were different with morbidities in visit 2. 
set.seed(123)
ancom.md1.visit2.output.ASV <-ancombc2(data=phylobj.visit2.filtered, tax_level=NULL, fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md1.visit2.output.ASV.res <- ancom.md1.visit2.output.ASV$res 
lapply(ancom.md1.visit2.output.ASV.res[,c(grep(pattern="diff_", colnames(ancom.md1.visit2.output.ASV.res)))], FUN=table)
#1 ASV different with age, one with laxative, one with mode of birth, none with morbidity

#cut down model for diarrhea
set.seed(123)
ancom.md2.visit2.diarrhea.output.ASV <-ancombc2(data=phylobj.visit2.filtered, tax_level=NULL, fix_formula="stool_age_days+infant.laxative + mode.of.delivery+infant.diarrhea", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit2.diarrhea.output.ASV.res <- ancom.md2.visit2.diarrhea.output.ASV$res 
lapply(ancom.md2.visit2.diarrhea.output.ASV.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit2.diarrhea.output.ASV.res)))], FUN=table)
#still no ASVs different with diarrhea in visit 2. 

#cut down model for fever
set.seed(123)
ancom.md2.visit2.fever.output.ASV <-ancombc2(data=phylobj.visit2.filtered, tax_level=NULL, fix_formula="stool_age_days+infant.laxative + mode.of.delivery+infant.fever", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit2.fever.output.ASV.res <- ancom.md2.visit2.fever.output.ASV$res 
lapply(ancom.md2.visit2.fever.output.ASV.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit2.fever.output.ASV.res)))], FUN=table)
#one ASV different by age, 3 by laxative use, one by mode of birth, and none by fever. 

#cut down model for vomiting
set.seed(123)
ancom.md2.visit2.vomit.output.ASV <-ancombc2(data=phylobj.visit2.filtered, tax_level=NULL, fix_formula="stool_age_days+infant.laxative + mode.of.delivery+infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit2.vomit.output.ASV.res <- ancom.md2.visit2.vomit.output.ASV$res 
lapply(ancom.md2.visit2.vomit.output.ASV.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit2.vomit.output.ASV.res)))], FUN=table)
#one ASV different by age, two by laxative, and 1 by mode of birth, and none by morbidity. 

#species level
set.seed(123)
ancom.md1.visit2.output.species <-ancombc2(data=phylobj.visit2.filtered, tax_level="Species", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md1.visit2.output.species.res <- ancom.md1.visit2.output.species$res 
lapply(ancom.md1.visit2.output.species.res[,c(grep(pattern="diff_", colnames(ancom.md1.visit2.output.species.res)))], FUN=table)
#one species different with probiotic, one with laxative,
#and none with morbidity. 
print(ancom.md1.visit2.output.species.res[ancom.md1.visit2.output.species.res$diff_infant.probioticyes==TRUE,])
#Lactobacillus rhamnosus was associated with probiotic use in visit 2. 
table(metadata.v2.visit2$infant.probiotic=="yes") #24 infants had probiotic

#cut down model for diarrhea
set.seed(123)
ancom.md2.visit2.diarrhea.output.species <-ancombc2(data=phylobj.visit2.filtered, tax_level="Species", fix_formula="infant.probiotic+infant.laxative +infant.diarrhea", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit2.diarrhea.output.species.res <- ancom.md2.visit2.diarrhea.output.species$res 
lapply(ancom.md2.visit2.diarrhea.output.species.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit2.diarrhea.output.species.res)))], FUN=table)
#one species was different with laxative, none with diarrhea

#cut down model for fever
set.seed(123)
ancom.md2.visit2.fever.output.species <-ancombc2(data=phylobj.visit2.filtered, tax_level="Species", fix_formula="infant.probiotic+infant.laxative +infant.fever", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit2.fever.output.species.res <- ancom.md2.visit2.fever.output.species$res 
lapply(ancom.md2.visit2.fever.output.species.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit2.fever.output.species.res)))], FUN=table)
#one species was different with laxative, one with probiotic, none with fever

#cut down model for vomit
set.seed(123)
ancom.md2.visit2.vomit.output.species <-ancombc2(data=phylobj.visit2.filtered, tax_level="Species", fix_formula="infant.probiotic+infant.laxative +infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit2.vomit.output.species.res <- ancom.md2.visit2.vomit.output.species$res 
lapply(ancom.md2.visit2.vomit.output.species.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit2.vomit.output.species.res)))], FUN=table)
#one species different with laxative, and none with vomiting. 

#genus level
set.seed(123)
ancom.md1.visit2.output.genus <-ancombc2(data=phylobj.visit2.filtered, tax_level="Genus", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md1.visit2.output.genus.res <- ancom.md1.visit2.output.genus$res 
lapply(ancom.md1.visit2.output.genus.res[,c(grep(pattern="diff_", colnames(ancom.md1.visit2.output.genus.res)))], FUN=table)
#several genera differed with laxative, none with morbidity. 

set.seed(123)
ancom.md2.visit2.diarrhea.output.genus <-ancombc2(data=phylobj.visit2.filtered, tax_level="Genus", fix_formula="infant.laxative+infant.diarrhea", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit2.diarrhea.output.genus.res <- ancom.md2.visit2.diarrhea.output.genus$res 
lapply(ancom.md2.visit2.diarrhea.output.genus.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit2.diarrhea.output.genus.res)))], FUN=table)
#none with diarrhea.

set.seed(123)
ancom.md2.visit2.fever.output.genus <-ancombc2(data=phylobj.visit2.filtered, tax_level="Genus", fix_formula="infant.laxative+infant.fever", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit2.fever.output.genus.res <- ancom.md2.visit2.fever.output.genus$res 
lapply(ancom.md2.visit2.fever.output.genus.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit2.fever.output.genus.res)))], FUN=table)
#none with fever.

set.seed(123)
ancom.md2.visit2.vomit.output.genus <-ancombc2(data=phylobj.visit2.filtered, tax_level="Genus", fix_formula="infant.laxative+infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit2.vomit.output.genus.res <- ancom.md2.visit2.vomit.output.genus$res 
lapply(ancom.md2.visit2.vomit.output.genus.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit2.vomit.output.genus.res)))], FUN=table)
#none with vomiting.

set.seed(123)
ancom.md1.visit2.output.family <-ancombc2(data=phylobj.visit2.filtered, tax_level="Family", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md1.visit2.output.family.res <- ancom.md1.visit2.output.family$res 
lapply(ancom.md1.visit2.output.family.res[,c(grep(pattern="diff_", colnames(ancom.md1.visit2.output.family.res)))], FUN=table)
#three families different with laxative, and none with morbidity. 

set.seed(123)
ancom.md2.visit2.diarrhea.output.family <-ancombc2(data=phylobj.visit2.filtered, tax_level="Family", fix_formula="infant.laxative+infant.diarrhea", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit2.diarrhea.output.family.res <- ancom.md2.visit2.diarrhea.output.family$res 
lapply(ancom.md2.visit2.diarrhea.output.family.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit2.diarrhea.output.family.res)))], FUN=table)
#none with diarrhea.

set.seed(123)
ancom.md2.visit2.fever.output.family <-ancombc2(data=phylobj.visit2.filtered, tax_level="Family", fix_formula="infant.laxative+infant.fever", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit2.fever.output.family.res <- ancom.md2.visit2.fever.output.family$res 
lapply(ancom.md2.visit2.fever.output.family.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit2.fever.output.family.res)))], FUN=table)
#none with fever.

set.seed(123)
ancom.md2.visit2.vomit.output.family <-ancombc2(data=phylobj.visit2.filtered, tax_level="Family", fix_formula="infant.laxative+infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit2.vomit.output.family.res <- ancom.md2.visit2.vomit.output.family$res 
lapply(ancom.md2.visit2.vomit.output.family.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit2.vomit.output.family.res)))], FUN=table)
#none with vomiting. 

set.seed(123)
ancom.md1.visit2.output.order <-ancombc2(data=phylobj.visit2.filtered, tax_level="Order", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md1.visit2.output.order.res <- ancom.md1.visit2.output.order$res 
lapply(ancom.md1.visit2.output.order.res[,c(grep(pattern="diff_", colnames(ancom.md1.visit2.output.order.res)))], FUN=table)
#two orders different with laxative

set.seed(123)
ancom.md2.visit2.diarrhea.output.order <-ancombc2(data=phylobj.visit2.filtered, tax_level="Order", fix_formula="infant.laxative+infant.diarrhea", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit2.diarrhea.output.order.res <- ancom.md2.visit2.diarrhea.output.order$res 
lapply(ancom.md2.visit2.diarrhea.output.order.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit2.diarrhea.output.order.res)))], FUN=table)
#no orders differentially abundant. 

set.seed(123)
ancom.md2.visit2.fever.output.order <-ancombc2(data=phylobj.visit2.filtered, tax_level="Order", fix_formula="infant.laxative+infant.fever", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit2.fever.output.order.res <- ancom.md2.visit2.fever.output.order$res 
lapply(ancom.md2.visit2.fever.output.order.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit2.fever.output.order.res)))], FUN=table)
#no orders differentially abundant. 

set.seed(123)
ancom.md2.visit2.vomit.output.order <-ancombc2(data=phylobj.visit2.filtered, tax_level="Order", fix_formula="infant.laxative+infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit2.vomit.output.order.res <- ancom.md2.visit2.vomit.output.order$res 
lapply(ancom.md2.visit2.vomit.output.order.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit2.vomit.output.order.res)))], FUN=table)
#no orders differentially abundant. 

#Step 3: determine if any taxa were different with morbidities in visit 3. 
set.seed(123)
ancom.md1.visit3.output.ASV <-ancombc2(data=phylobj.visit3.filtered, tax_level=NULL, fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md1.visit3.output.ASV.res <- ancom.md1.visit3.output.ASV$res 
lapply(ancom.md1.visit3.output.ASV.res[,c(grep(pattern="diff_", colnames(ancom.md1.visit3.output.ASV.res)))], FUN=table)
#three ASVs different with age, none with morbidity. 

set.seed(123)
ancom.md2.visit3.diarrhea.output.ASV <-ancombc2(data=phylobj.visit3.filtered, tax_level=NULL, fix_formula="stool_age_days +infant.diarrhea", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit3.diarrhea.output.ASV.res <- ancom.md2.visit3.diarrhea.output.ASV$res 
lapply(ancom.md2.visit3.diarrhea.output.ASV.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit3.diarrhea.output.ASV.res)))], FUN=table)
#none by diarrhea

set.seed(123)
ancom.md2.visit3.fever.output.ASV <-ancombc2(data=phylobj.visit3.filtered, tax_level=NULL, fix_formula="stool_age_days +infant.fever", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit3.fever.output.ASV.res <- ancom.md2.visit3.fever.output.ASV$res 
lapply(ancom.md2.visit3.fever.output.ASV.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit3.fever.output.ASV.res)))], FUN=table)
#none by fever

set.seed(123)
ancom.md2.visit3.vomit.output.ASV <-ancombc2(data=phylobj.visit3.filtered, tax_level=NULL, fix_formula="stool_age_days +infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit3.vomit.output.ASV.res <- ancom.md2.visit3.vomit.output.ASV$res 
lapply(ancom.md2.visit3.vomit.output.ASV.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit3.vomit.output.ASV.res)))], FUN=table)
#none by  vomiting. 

set.seed(123)
ancom.md1.visit3.output.species <-ancombc2(data=phylobj.visit3.filtered, tax_level="Species", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md1.visit3.output.species.res <- ancom.md1.visit3.output.species$res 
lapply(ancom.md1.visit3.output.species.res[,c(grep(pattern="diff_", colnames(ancom.md1.visit3.output.species.res)))], FUN=table)
#three species different with age, none with morbidity. 

set.seed(123)
ancom.md2.visit3.diarrhea.output.species <-ancombc2(data=phylobj.visit3.filtered, tax_level="Species", fix_formula="stool_age_days+infant.diarrhea", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit3.diarrhea.output.species.res <- ancom.md2.visit3.diarrhea.output.species$res 
lapply(ancom.md2.visit3.diarrhea.output.species.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit3.diarrhea.output.species.res)))], FUN=table)
#none with diarrhea

set.seed(123)
ancom.md2.visit3.fever.output.species <-ancombc2(data=phylobj.visit3.filtered, tax_level="Species", fix_formula="stool_age_days+infant.fever", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit3.fever.output.species.res <- ancom.md2.visit3.fever.output.species$res 
lapply(ancom.md2.visit3.fever.output.species.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit3.fever.output.species.res)))], FUN=table)
#none with fever

set.seed(123)
ancom.md2.visit3.vomit.output.species <-ancombc2(data=phylobj.visit3.filtered, tax_level="Species", fix_formula="stool_age_days+infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit3.vomit.output.species.res <- ancom.md2.visit3.vomit.output.species$res 
lapply(ancom.md2.visit3.vomit.output.species.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit3.vomit.output.species.res)))], FUN=table)
#none with vomit

set.seed(123)
ancom.md1.visit3.output.genus <-ancombc2(data=phylobj.visit3.filtered, tax_level="Genus", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md1.visit3.output.genus.res <- ancom.md1.visit3.output.genus$res 
lapply(ancom.md1.visit3.output.genus.res[,c(grep(pattern="diff_", colnames(ancom.md1.visit3.output.genus.res)))], FUN=table)
#none

set.seed(123)
ancom.md2.visit3.diarrhea.output.genus <-ancombc2(data=phylobj.visit3.filtered, tax_level="Genus", fix_formula="infant.diarrhea", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit3.diarrhea.output.genus.res <- ancom.md2.visit3.diarrhea.output.genus$res 
lapply(ancom.md2.visit3.diarrhea.output.genus.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit3.diarrhea.output.genus.res)))], FUN=table)
#none

set.seed(123)
ancom.md2.visit3.fever.output.genus <-ancombc2(data=phylobj.visit3.filtered, tax_level="Genus", fix_formula="infant.fever", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit3.fever.output.genus.res <- ancom.md2.visit3.fever.output.genus$res 
lapply(ancom.md2.visit3.fever.output.genus.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit3.fever.output.genus.res)))], FUN=table)
#none

set.seed(123)
ancom.md2.visit3.vomit.output.genus <-ancombc2(data=phylobj.visit3.filtered, tax_level="Genus", fix_formula="infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit3.vomit.output.genus.res <- ancom.md2.visit3.vomit.output.genus$res 
lapply(ancom.md2.visit3.vomit.output.genus.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit3.vomit.output.genus.res)))], FUN=table)
#none

set.seed(123)
ancom.md1.visit3.output.family <-ancombc2(data=phylobj.visit3.filtered, tax_level="Family", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md1.visit3.output.family.res <- ancom.md1.visit3.output.family$res 
lapply(ancom.md1.visit3.output.family.res[,c(grep(pattern="diff_", colnames(ancom.md1.visit3.output.family.res)))], FUN=table)
#none

set.seed(123)
ancom.md2.visit3.diarrhea.output.family <-ancombc2(data=phylobj.visit3.filtered, tax_level="Family", fix_formula="infant.diarrhea", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit3.diarrhea.output.family.res <- ancom.md2.visit3.diarrhea.output.family$res 
lapply(ancom.md2.visit3.diarrhea.output.family.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit3.diarrhea.output.family.res)))], FUN=table)
#none

set.seed(123)
ancom.md2.visit3.fever.output.family <-ancombc2(data=phylobj.visit3.filtered, tax_level="Family", fix_formula="infant.fever", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit3.fever.output.family.res <- ancom.md2.visit3.fever.output.family$res 
lapply(ancom.md2.visit3.fever.output.family.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit3.fever.output.family.res)))], FUN=table)
#none

set.seed(123)
ancom.md2.visit3.vomit.output.family <-ancombc2(data=phylobj.visit3.filtered, tax_level="Family", fix_formula="infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit3.vomit.output.family.res <- ancom.md2.visit3.vomit.output.family$res 
lapply(ancom.md2.visit3.vomit.output.family.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit3.vomit.output.family.res)))], FUN=table)
#none

set.seed(123)
ancom.md1.visit3.output.order <-ancombc2(data=phylobj.visit3.filtered, tax_level="Order", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md1.visit3.output.order.res <- ancom.md1.visit3.output.order$res 
lapply(ancom.md1.visit3.output.order.res[,c(grep(pattern="diff_", colnames(ancom.md1.visit3.output.order.res)))], FUN=table)
#one order different with age, one with antibiotics, none with morbidity. 

set.seed(123)
ancom.md2.visit3.diarrhea.output.order <-ancombc2(data=phylobj.visit3.filtered, tax_level="Order", fix_formula="stool_age_days + infant.antibiotic+infant.diarrhea", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit3.diarrhea.output.order.res <- ancom.md2.visit3.diarrhea.output.order$res 
lapply(ancom.md2.visit3.diarrhea.output.order.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit3.diarrhea.output.order.res)))], FUN=table)
#none

set.seed(123)
ancom.md2.visit3.fever.output.order <-ancombc2(data=phylobj.visit3.filtered, tax_level="Order", fix_formula="stool_age_days + infant.antibiotic+infant.fever", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit3.fever.output.order.res <- ancom.md2.visit3.fever.output.order$res 
lapply(ancom.md2.visit3.fever.output.order.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit3.fever.output.order.res)))], FUN=table)
#one order differed with fever in visit 3. 
print(ancom.md2.visit3.fever.output.order.res[ancom.md2.visit3.fever.output.order.res$diff_infant.feveryes=="TRUE",])
#Bacteroidales was negatively associated with fever in visit 3. 


set.seed(123)
ancom.md2.visit3.vomit.output.order <-ancombc2(data=phylobj.visit3.filtered, tax_level="Order", fix_formula="stool_age_days + infant.antibiotic+infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit3.vomit.output.order.res <- ancom.md2.visit3.vomit.output.order$res 
lapply(ancom.md2.visit3.vomit.output.order.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit3.vomit.output.order.res)))], FUN=table)
#one order differed with antibiotic use. 

#since bacteroidiales were lower with fever, and other have found vaccine response
#to be higher with lower Bacteroidetes in the gut, I am going to look at 
#bacteroidetes (phylum level) association with fever. 

set.seed(123)
ancom.md1.visit3.output.phylum <-ancombc2(data=phylobj.visit3.filtered, tax_level="Phylum", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md1.visit3.output.phylum.res <- ancom.md1.visit3.output.phylum$res 
lapply(ancom.md1.visit3.output.phylum.res[,c(grep(pattern="diff_", colnames(ancom.md1.visit3.output.phylum.res)))], FUN=table)
#one phylum associated with age, none associated with morbidity. 

#cut down model. 
set.seed(123)
ancom.md2.visit3.fever.output.phylum <-ancombc2(data=phylobj.visit3.filtered, tax_level="Phylum", fix_formula="stool_age_days + infant.fever", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit3.fever.output.phylum.res <- ancom.md2.visit3.fever.output.phylum$res 
lapply(ancom.md2.visit3.fever.output.phylum.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit3.fever.output.phylum.res)))], FUN=table)
#one differed with age, and one with fever. 
print(ancom.md2.visit3.fever.output.phylum.res[ancom.md2.visit3.fever.output.phylum.res$diff_infant.feveryes=="TRUE",])
#Bacteroidota was negatively associated with fever in visit 3. 

#Step 4: determine if any taxa were different with morbidities in visit 4. 

set.seed(123)
ancom.md1.visit4.output.ASV <-ancombc2(data=phylobj.visit4.filtered, tax_level=NULL, fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md1.visit4.output.ASV.res <- ancom.md1.visit4.output.ASV$res 
lapply(ancom.md1.visit4.output.ASV.res[,c(grep(pattern="diff_", colnames(ancom.md1.visit4.output.ASV.res)))], FUN=table)
#none

set.seed(123)
ancom.md2.visit4.diarrhea.output.ASV <-ancombc2(data=phylobj.visit4.filtered, tax_level=NULL, fix_formula="infant.diarrhea", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit4.diarrhea.output.ASV.res <- ancom.md2.visit4.diarrhea.output.ASV$res 
lapply(ancom.md2.visit4.diarrhea.output.ASV.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit4.diarrhea.output.ASV.res)))], FUN=table)
#none

set.seed(123)
ancom.md2.visit4.fever.output.ASV <-ancombc2(data=phylobj.visit4.filtered, tax_level=NULL, fix_formula="infant.fever", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit4.fever.output.ASV.res <- ancom.md2.visit4.fever.output.ASV$res 
lapply(ancom.md2.visit4.fever.output.ASV.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit4.fever.output.ASV.res)))], FUN=table)
#none

set.seed(123)
ancom.md2.visit4.vomit.output.ASV <-ancombc2(data=phylobj.visit4.filtered, tax_level=NULL, fix_formula="infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit4.vomit.output.ASV.res <- ancom.md2.visit4.vomit.output.ASV$res 
lapply(ancom.md2.visit4.vomit.output.ASV.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit4.vomit.output.ASV.res)))], FUN=table)
#none

set.seed(123)
ancom.md1.visit4.output.species <-ancombc2(data=phylobj.visit4.filtered, tax_level="Species", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md1.visit4.output.species.res <- ancom.md1.visit4.output.species$res 
lapply(ancom.md1.visit4.output.species.res[,c(grep(pattern="diff_", colnames(ancom.md1.visit4.output.species.res)))], FUN=table)
#one species different with probiotic use, none with morbidity. 

set.seed(123)
ancom.md2.visit4.diarrhea.output.species <-ancombc2(data=phylobj.visit4.filtered, tax_level="Species", fix_formula="infant.probiotic + infant.diarrhea", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit4.diarrhea.output.species.res <- ancom.md2.visit4.diarrhea.output.species$res 
lapply(ancom.md2.visit4.diarrhea.output.species.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit4.diarrhea.output.species.res)))], FUN=table)
#one species different with probiotic use, none with diarrhea. 

set.seed(123)
ancom.md2.visit4.fever.output.species <-ancombc2(data=phylobj.visit4.filtered, tax_level="Species", fix_formula="infant.probiotic + infant.fever", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit4.fever.output.species.res <- ancom.md2.visit4.fever.output.species$res 
lapply(ancom.md2.visit4.fever.output.species.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit4.fever.output.species.res)))], FUN=table)
#one species different with probiotic, none with fever. 

set.seed(123)
ancom.md2.visit4.vomit.output.species <-ancombc2(data=phylobj.visit4.filtered, tax_level="Species", fix_formula="infant.probiotic + infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit4.vomit.output.species.res <- ancom.md2.visit4.vomit.output.species$res 
lapply(ancom.md2.visit4.vomit.output.species.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit4.vomit.output.species.res)))], FUN=table)
#one species different with probiotic use, none with vomit. 

set.seed(123)
ancom.md1.visit4.output.genus <-ancombc2(data=phylobj.visit4.filtered, tax_level="Genus", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md1.visit4.output.genus.res <- ancom.md1.visit4.output.genus$res 
lapply(ancom.md1.visit4.output.genus.res[,c(grep(pattern="diff_", colnames(ancom.md1.visit4.output.genus.res)))], FUN=table)
#two genera different with age, none with morbidity. 

set.seed(123)
ancom.md2.visit4.diarrhea.output.genus <-ancombc2(data=phylobj.visit4.filtered, tax_level="Genus", fix_formula="stool_age_days + infant.diarrhea", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit4.diarrhea.output.genus.res <- ancom.md2.visit4.diarrhea.output.genus$res 
lapply(ancom.md2.visit4.diarrhea.output.genus.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit4.diarrhea.output.genus.res)))], FUN=table)
#none with diarrhea.

set.seed(123)
ancom.md2.visit4.fever.output.genus <-ancombc2(data=phylobj.visit4.filtered, tax_level="Genus", fix_formula="stool_age_days + infant.fever", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit4.fever.output.genus.res <- ancom.md2.visit4.fever.output.genus$res 
lapply(ancom.md2.visit4.fever.output.genus.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit4.fever.output.genus.res)))], FUN=table)
#none with fever. 

set.seed(123)
ancom.md2.visit4.vomit.output.genus <-ancombc2(data=phylobj.visit4.filtered, tax_level="Genus", fix_formula="stool_age_days + infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit4.vomit.output.genus.res <- ancom.md2.visit4.vomit.output.genus$res 
lapply(ancom.md2.visit4.vomit.output.genus.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit4.vomit.output.genus.res)))], FUN=table)
#none with vomiting. 

set.seed(123)
ancom.md1.visit4.output.family <-ancombc2(data=phylobj.visit4.filtered, tax_level="Family", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md1.visit4.output.family.res <- ancom.md1.visit4.output.family$res 
lapply(ancom.md1.visit4.output.family.res[,c(grep(pattern="diff_", colnames(ancom.md1.visit4.output.family.res)))], FUN=table)
#one family different with age, none with morbidity. 

set.seed(123)
ancom.md2.visit4.diarrhea.output.family <-ancombc2(data=phylobj.visit4.filtered, tax_level="Family", fix_formula="stool_age_days + infant.diarrhea", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit4.diarrhea.output.family.res <- ancom.md2.visit4.diarrhea.output.family$res 
lapply(ancom.md2.visit4.diarrhea.output.family.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit4.diarrhea.output.family.res)))], FUN=table)
#none with diarrhea. 

set.seed(123)
ancom.md2.visit4.fever.output.family <-ancombc2(data=phylobj.visit4.filtered, tax_level="Family", fix_formula="stool_age_days + infant.fever", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit4.fever.output.family.res <- ancom.md2.visit4.fever.output.family$res 
lapply(ancom.md2.visit4.fever.output.family.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit4.fever.output.family.res)))], FUN=table)
#none with fever. 

set.seed(123)
ancom.md2.visit4.vomit.output.family <-ancombc2(data=phylobj.visit4.filtered, tax_level="Family", fix_formula="stool_age_days + infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit4.vomit.output.family.res <- ancom.md2.visit4.vomit.output.family$res 
lapply(ancom.md2.visit4.vomit.output.family.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit4.vomit.output.family.res)))], FUN=table)
#none with vomiting. 

set.seed(123)
ancom.md1.visit4.output.order <-ancombc2(data=phylobj.visit4.filtered, tax_level="Order", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md1.visit4.output.order.res <- ancom.md1.visit4.output.order$res 
lapply(ancom.md1.visit4.output.order.res[,c(grep(pattern="diff_", colnames(ancom.md1.visit4.output.order.res)))], FUN=table)
#two orders different with age, none with morbidity. 

set.seed(123)
ancom.md2.visit4.diarrhea.output.order <-ancombc2(data=phylobj.visit4.filtered, tax_level="Order", fix_formula="stool_age_days + infant.diarrhea", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit4.diarrhea.output.order.res <- ancom.md2.visit4.diarrhea.output.order$res 
lapply(ancom.md2.visit4.diarrhea.output.order.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit4.diarrhea.output.order.res)))], FUN=table)
#none with diarrhea. 

set.seed(123)
ancom.md2.visit4.fever.output.order <-ancombc2(data=phylobj.visit4.filtered, tax_level="Order", fix_formula="stool_age_days + infant.fever", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit4.fever.output.order.res <- ancom.md2.visit4.fever.output.order$res 
lapply(ancom.md2.visit4.fever.output.order.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit4.fever.output.order.res)))], FUN=table)
#none with fever. 

set.seed(123)
ancom.md2.visit4.vomit.output.order <-ancombc2(data=phylobj.visit4.filtered, tax_level="Order", fix_formula="stool_age_days + infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.visit4.vomit.output.order.res <- ancom.md2.visit4.vomit.output.order$res 
lapply(ancom.md2.visit4.vomit.output.order.res[,c(grep(pattern="diff_", colnames(ancom.md2.visit4.vomit.output.order.res)))], FUN=table)
#none with vomiting. 

#Part 4: exploring the hypothesis that the stool microbiota in visit 2 is predictive
#of diarrhea later on in visits 3 and/or 4. Using ANCOM-BC2 to determine if any 
#taxa in visit 2 are associated with morbidity later on in visits 3 or 4.

#Step 1: ASVs in visit 2 associated with morbidity in visit 3, 4. 

#full model for diarrhea
set.seed(123)
ancom.visit2.predicting.diarrhea.md1.output.ASV <-ancombc2(data=phylobj.visit2.filtered, tax_level=NULL, fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit+infant.diarrhea.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.diarrhea.md1.output.ASV.res <- ancom.visit2.predicting.diarrhea.md1.output.ASV$res 
lapply(ancom.visit2.predicting.diarrhea.md1.output.ASV.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.diarrhea.md1.output.ASV.res)))], FUN=table)
#one ASV different with age, one with laxative, one ASV different with mode of birth
#and no ASVs different with diarrhea 

#cut down model.
set.seed(123)
ancom.visit2.predicting.diarrhea.md2.output.ASV <-ancombc2(data=phylobj.visit2.filtered, tax_level=NULL, fix_formula="stool_age_days +infant.laxative+ mode.of.delivery+infant.diarrhea.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.diarrhea.md2.output.ASV.res <- ancom.visit2.predicting.diarrhea.md2.output.ASV$res 
lapply(ancom.visit2.predicting.diarrhea.md2.output.ASV.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.diarrhea.md2.output.ASV.res)))], FUN=table)
#none associated with diarrhea in visit 3, 4. 

#full model for fever
set.seed(123)
ancom.visit2.predicting.fever.md1.output.ASV <-ancombc2(data=phylobj.visit2.filtered, tax_level=NULL, fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit+infant.fever.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.fever.md1.output.ASV.res <- ancom.visit2.predicting.fever.md1.output.ASV$res 
lapply(ancom.visit2.predicting.fever.md1.output.ASV.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.fever.md1.output.ASV.res)))], FUN=table)
#one ASV different with age, one with laxative, one ASV different with mode of birth
#and no ASVs different with fever

#cut-down model for fever
set.seed(123)
ancom.visit2.predicting.fever.md2.output.ASV <-ancombc2(data=phylobj.visit2.filtered, tax_level=NULL, fix_formula="stool_age_days +infant.laxative + mode.of.delivery+infant.fever.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.fever.md2.output.ASV.res <- ancom.visit2.predicting.fever.md2.output.ASV$res 
lapply(ancom.visit2.predicting.fever.md2.output.ASV.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.fever.md2.output.ASV.res)))], FUN=table)
#none associated with fever.

#full model for vomiting
set.seed(123)
ancom.visit2.predicting.vomit.md1.output.ASV <-ancombc2(data=phylobj.visit2.filtered, tax_level=NULL, fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit+infant.vomit.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.vomit.md1.output.ASV.res <- ancom.visit2.predicting.vomit.md1.output.ASV$res 
lapply(ancom.visit2.predicting.vomit.md1.output.ASV.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.vomit.md1.output.ASV.res)))], FUN=table)
#one ASV differed with age, one differed with laxative, and one with mode of birth, 
#and none with vomiting. 

#cut-down model for vomiting
set.seed(123)
ancom.visit2.predicting.vomit.md2.output.ASV <-ancombc2(data=phylobj.visit2.filtered, tax_level=NULL, fix_formula="stool_age_days +infant.laxative + mode.of.delivery+infant.vomit.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.vomit.md2.output.ASV.res <- ancom.visit2.predicting.vomit.md2.output.ASV$res 
lapply(ancom.visit2.predicting.vomit.md2.output.ASV.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.vomit.md2.output.ASV.res)))], FUN=table)
#none associated with vomit in visit 3, 4. 

#Step 2: species in visit 2 associated with morbidity in visit 3, 4.


#full model for diarrhea
set.seed(123)
ancom.visit2.predicting.diarrhea.md1.output.species <-ancombc2(data=phylobj.visit2.filtered, tax_level="Species", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit+infant.diarrhea.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.diarrhea.md1.output.species.res <- ancom.visit2.predicting.diarrhea.md1.output.species$res 
lapply(ancom.visit2.predicting.diarrhea.md1.output.species.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.diarrhea.md1.output.species.res)))], FUN=table)
#one species differed with probiotic, one with laxative, and none with diarrhea. 

#cut down model for diarrhea.
set.seed(123)
ancom.visit2.predicting.diarrhea.md2.output.species <-ancombc2(data=phylobj.visit2.filtered, tax_level="Species", fix_formula="infant.probiotic+infant.laxative+infant.diarrhea.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.diarrhea.md2.output.species.res <- ancom.visit2.predicting.diarrhea.md2.output.species$res 
lapply(ancom.visit2.predicting.diarrhea.md2.output.species.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.diarrhea.md2.output.species.res)))], FUN=table)
#none differed with diarrhea.

#full model for fever
set.seed(123)
ancom.visit2.predicting.fever.md1.output.species <-ancombc2(data=phylobj.visit2.filtered, tax_level="Species", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit+infant.fever.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.fever.md1.output.species.res <- ancom.visit2.predicting.fever.md1.output.species$res 
lapply(ancom.visit2.predicting.fever.md1.output.species.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.fever.md1.output.species.res)))], FUN=table)
#one species differed with probiotic, and one with laxative, and none with fever in visit 3,4. 

#cut down model for fever
set.seed(123)
ancom.visit2.predicting.fever.md2.output.species <-ancombc2(data=phylobj.visit2.filtered, tax_level="Species", fix_formula="infant.probiotic+infant.laxative+infant.fever.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.fever.md2.output.species.res <- ancom.visit2.predicting.fever.md2.output.species$res 
lapply(ancom.visit2.predicting.fever.md2.output.species.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.fever.md2.output.species.res)))], FUN=table)
#none differed with fever in visit 3, 4. 

#full model for vomit
set.seed(123)
ancom.visit2.predicting.vomit.md1.output.species <-ancombc2(data=phylobj.visit2.filtered, tax_level="Species", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit+infant.vomit.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.vomit.md1.output.species.res <- ancom.visit2.predicting.vomit.md1.output.species$res 
lapply(ancom.visit2.predicting.vomit.md1.output.species.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.vomit.md1.output.species.res)))], FUN=table)
#one species differed with laxative, and none with vomiting. 

#cut down model for vomit
set.seed(123)
ancom.visit2.predicting.vomit.md2.output.species <-ancombc2(data=phylobj.visit2.filtered, tax_level="Species", fix_formula="infant.laxative+infant.vomit.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.vomit.md2.output.species.res <- ancom.visit2.predicting.vomit.md2.output.species$res 
lapply(ancom.visit2.predicting.vomit.md2.output.species.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.vomit.md2.output.species.res)))], FUN=table)
#none different with vomiting in visit 3, 4. 

#Step 3: genera in visit 2 associated with morbidity in visit 3, 4.

#full model for diarrhea
set.seed(123)
ancom.visit2.predicting.diarrhea.md1.output.genus <-ancombc2(data=phylobj.visit2.filtered, tax_level="Genus", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit+infant.diarrhea.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.diarrhea.md1.output.genus.res <- ancom.visit2.predicting.diarrhea.md1.output.genus$res 
lapply(ancom.visit2.predicting.diarrhea.md1.output.genus.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.diarrhea.md1.output.genus.res)))], FUN=table)
#4 genera differed with laxative, and none with diarrhea. 

#cut-down model for diarrhea
set.seed(123)
ancom.visit2.predicting.diarrhea.md2.output.genus <-ancombc2(data=phylobj.visit2.filtered, tax_level="Genus", fix_formula="infant.laxative +infant.diarrhea.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.diarrhea.md2.output.genus.res <- ancom.visit2.predicting.diarrhea.md2.output.genus$res 
lapply(ancom.visit2.predicting.diarrhea.md2.output.genus.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.diarrhea.md2.output.genus.res)))], FUN=table)
#none differed with diarrhea.

#full model for fever
set.seed(123)
ancom.visit2.predicting.fever.md1.output.genus <-ancombc2(data=phylobj.visit2.filtered, tax_level="Genus", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit+infant.fever.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.fever.md1.output.genus.res <- ancom.visit2.predicting.fever.md1.output.genus$res 
lapply(ancom.visit2.predicting.fever.md1.output.genus.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.fever.md1.output.genus.res)))], FUN=table)
#4 genera differed with laxative, none with fever. 

#cut-down model for fever
set.seed(123)
ancom.visit2.predicting.fever.md2.output.genus <-ancombc2(data=phylobj.visit2.filtered, tax_level="Genus", fix_formula="infant.laxative +infant.fever.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.fever.md2.output.genus.res <- ancom.visit2.predicting.fever.md2.output.genus$res 
lapply(ancom.visit2.predicting.fever.md2.output.genus.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.fever.md2.output.genus.res)))], FUN=table)
#none differed with fever in visit 3, 4. 

#full model for vomit
set.seed(123)
ancom.visit2.predicting.vomit.md1.output.genus <-ancombc2(data=phylobj.visit2.filtered, tax_level="Genus", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit+infant.vomit.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.vomit.md1.output.genus.res <- ancom.visit2.predicting.vomit.md1.output.genus$res 
lapply(ancom.visit2.predicting.vomit.md1.output.genus.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.vomit.md1.output.genus.res)))], FUN=table)
#4 genera differed with laxative, none with vomit. 

#cut-down model for vomit
set.seed(123)
ancom.visit2.predicting.vomit.md2.output.genus <-ancombc2(data=phylobj.visit2.filtered, tax_level="Genus", fix_formula="infant.laxative +infant.vomit.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.vomit.md2.output.genus.res <- ancom.visit2.predicting.vomit.md2.output.genus$res 
lapply(ancom.visit2.predicting.vomit.md2.output.genus.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.vomit.md2.output.genus.res)))], FUN=table)
#none differed with vomit in visit 3, 4. 

#Step 4: families in visit 2 associated with morbidity in visit 3, 4.

#full model for diarrhea
set.seed(123)
ancom.visit2.predicting.diarrhea.md1.output.family <-ancombc2(data=phylobj.visit2.filtered, tax_level="Family", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit+infant.diarrhea.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.diarrhea.md1.output.family.res <- ancom.visit2.predicting.diarrhea.md1.output.family$res 
lapply(ancom.visit2.predicting.diarrhea.md1.output.family.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.diarrhea.md1.output.family.res)))], FUN=table)
#three families differed with laxative, none with diarrhea. 

#cut down model for diarrhea
set.seed(123)
ancom.visit2.predicting.diarrhea.md2.output.family <-ancombc2(data=phylobj.visit2.filtered, tax_level="Family", fix_formula="infant.laxative+infant.diarrhea.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.diarrhea.md2.output.family.res <- ancom.visit2.predicting.diarrhea.md2.output.family$res 
lapply(ancom.visit2.predicting.diarrhea.md2.output.family.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.diarrhea.md2.output.family.res)))], FUN=table)
#none differed with diarrhea in visit 3, 4. 

#full model for fever
set.seed(123)
ancom.visit2.predicting.fever.md1.output.family <-ancombc2(data=phylobj.visit2.filtered, tax_level="Family", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit+infant.fever.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.fever.md1.output.family.res <- ancom.visit2.predicting.fever.md1.output.family$res 
lapply(ancom.visit2.predicting.fever.md1.output.family.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.fever.md1.output.family.res)))], FUN=table)
#three families differed with laxative, none with fever. 

#cut down model for fever
set.seed(123)
ancom.visit2.predicting.fever.md2.output.family <-ancombc2(data=phylobj.visit2.filtered, tax_level="Family", fix_formula="infant.laxative+infant.fever.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.fever.md2.output.family.res <- ancom.visit2.predicting.fever.md2.output.family$res 
lapply(ancom.visit2.predicting.fever.md2.output.family.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.fever.md2.output.family.res)))], FUN=table)
#none differed with fever. 

#full model for vomit
set.seed(123)
ancom.visit2.predicting.vomit.md1.output.family <-ancombc2(data=phylobj.visit2.filtered, tax_level="Family", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit+infant.vomit.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.vomit.md1.output.family.res <- ancom.visit2.predicting.vomit.md1.output.family$res 
lapply(ancom.visit2.predicting.vomit.md1.output.family.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.vomit.md1.output.family.res)))], FUN=table)
#three families differed with laxative, none differed with vomiting

#cut-down model for vomit
set.seed(123)
ancom.visit2.predicting.vomit.md2.output.family <-ancombc2(data=phylobj.visit2.filtered, tax_level="Family", fix_formula="infant.laxative+infant.vomit.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.vomit.md2.output.family.res <- ancom.visit2.predicting.vomit.md2.output.family$res 
lapply(ancom.visit2.predicting.vomit.md2.output.family.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.vomit.md2.output.family.res)))], FUN=table)
#none differed with vomiting. 

#Step 5: orders in visit 2 associated with morbidity in visit 3, 4.

#full model for diarrhea
set.seed(123)
ancom.visit2.predicting.diarrhea.md1.output.order <-ancombc2(data=phylobj.visit2.filtered, tax_level="Order", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit+infant.diarrhea.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.diarrhea.md1.output.order.res <- ancom.visit2.predicting.diarrhea.md1.output.order$res 
lapply(ancom.visit2.predicting.diarrhea.md1.output.order.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.diarrhea.md1.output.order.res)))], FUN=table)
#two orders were different with laxative, none with diarrhea. 

#cut-down model for diarrhea
set.seed(123)
ancom.visit2.predicting.diarrhea.md2.output.order <-ancombc2(data=phylobj.visit2.filtered, tax_level="Order", fix_formula="infant.laxative+infant.diarrhea.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.diarrhea.md2.output.order.res <- ancom.visit2.predicting.diarrhea.md2.output.order$res 
lapply(ancom.visit2.predicting.diarrhea.md2.output.order.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.diarrhea.md2.output.order.res)))], FUN=table)
#none differed with diarrhea. 

#full model for fever
set.seed(123)
ancom.visit2.predicting.fever.md1.output.order <-ancombc2(data=phylobj.visit2.filtered, tax_level="Order", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit+infant.fever.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.fever.md1.output.order.res <- ancom.visit2.predicting.fever.md1.output.order$res 
lapply(ancom.visit2.predicting.fever.md1.output.order.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.fever.md1.output.order.res)))], FUN=table)
#one order differed with laxative, none with fever. 

#cut down for fever
set.seed(123)
ancom.visit2.predicting.fever.md2.output.order <-ancombc2(data=phylobj.visit2.filtered, tax_level="Order", fix_formula="infant.laxative+infant.fever.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.fever.md2.output.order.res <- ancom.visit2.predicting.fever.md2.output.order$res 
lapply(ancom.visit2.predicting.fever.md2.output.order.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.fever.md2.output.order.res)))], FUN=table)
#none differed with fever. 

#full model for vomit
set.seed(123)
ancom.visit2.predicting.vomit.md1.output.order <-ancombc2(data=phylobj.visit2.filtered, tax_level="Order", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+first.pregnancy+household.size+ExclusiveBreastfeeding.before.4months+infant.diarrhea+infant.fever+infant.vomit+infant.vomit.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.vomit.md1.output.order.res <- ancom.visit2.predicting.vomit.md1.output.order$res 
lapply(ancom.visit2.predicting.vomit.md1.output.order.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.vomit.md1.output.order.res)))], FUN=table)
#one order differed with laxative, none with vomiting in visit 3, 4. 

#cut down model for vomit
set.seed(123)
ancom.visit2.predicting.vomit.md2.output.order <-ancombc2(data=phylobj.visit2.filtered, tax_level="Order", fix_formula="infant.laxative+infant.vomit.in.visit3.or.visit4", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.visit2.predicting.vomit.md2.output.order.res <- ancom.visit2.predicting.vomit.md2.output.order$res 
lapply(ancom.visit2.predicting.vomit.md2.output.order.res[,c(grep(pattern="diff_", colnames(ancom.visit2.predicting.vomit.md2.output.order.res)))], FUN=table)
#none differed with vomit. 

#Part 5: determine if Granulicatella is differentially abundant by diarrhea
#when using rarefied counts and non-parametric testing

feature.table.rarefied <- read_qza("../feature-table-rarefied-6963.qza")
feature.table.rarefied <- feature.table.rarefied$data

taxonomy <- read_qza("../taxonomy-classification.qza")
taxonomy2 <- as.data.frame(taxonomy$data) %>% column_to_rownames("Feature.ID") %>% separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>% as.matrix()
#remove "Confidence" column
taxonomy2 <- taxonomy2[,-8]

#trim taxonomy to include only ASVs in the ASV table
taxonomy3 <- taxonomy2[c(rownames(taxonomy2)%in%rownames(feature.table)),]

phylobj.rarefied <- phyloseq(otu_table(feature.table.rarefied, taxa_are_rows=T),tax_table(taxonomy3), sample_data(metadata.v2 %>% as.data.frame() %>% column_to_rownames("sample.id")))
#get rid of the taxa that are not present at all
phylobj.rarefied<- filter_taxa(phylobj.rarefied, function(x) sum(x)>0, TRUE)
#look at the number of ASVs among the feature table
dim(otu_table(phylobj.rarefied))

#collapse the counts to genus level designations. 
phylobj.rarefied.genus <- tax_glom(phylobj.rarefied, 'Genus', NArm=F)

melted.phylgenus.rarefied <- psmelt(phylobj.rarefied.genus)
unique(melted.phylgenus.rarefied$Genus)
granulicatella.rarefied <- melted.phylgenus.rarefied[c(grep(pattern="Granulica",melted.phylgenus.rarefied$Genus)),]
table(granulicatella.rarefied$Abundance>0) #60 samples had a count above zero. 

wilcox.test(x=(granulicatella.rarefied[granulicatella.rarefied$infant.diarrhea=="yes","Abundance"]/6963*100), y= (granulicatella.rarefied[granulicatella.rarefied$infant.diarrhea=="no","Abundance"]/6963*100), paired=F, alternative="greater")
#significant by wilcox test, with diarrhea group greater than healthy group. 
granulicatella.rarefied$Relative.Abundance <- (granulicatella.rarefied$Abundance/6963)*100

ggplot(aes(x=infant.diarrhea, y=Relative.Abundance, color=as.factor(infant.diarrhea)), data=granulicatella.rarefied)+geom_boxplot()+xlab("diarrhea")+ylab("Granulicatella relative abundance %")+scale_x_discrete(labels=c("no\nn=287", "yes\nn=40"))+scale_color_manual(values=c(no="black", yes="#F8766D"))+theme(legend.position="none")
ggsave("granulicatella.vs.diarrhea.jpeg", dpi=600, width=5, height=5, plot=last_plot())

#also do the graph with non-rarefied counts but with relative abundance
phylobj.genus <- tax_glom(phylobj.filtered, 'Genus', NArm=F)

phylobj.genus.relabund <- transform_sample_counts(phylobj.genus, function(x) x/sum(x) * 100)

melted.phylgenus.relabund <- psmelt(phylobj.genus.relabund)
unique(melted.phylgenus.relabund$Genus)
granulicatella <- melted.phylgenus.relabund[c(grep(pattern="Granulica",melted.phylgenus.relabund$Genus)),]
table(granulicatella$Abundance>0) #65 samples had a count above zero. 

ggplot(aes(x=infant.diarrhea, y=Abundance), data=granulicatella)+geom_boxplot()+xlab("did infant have diarrhea during visit?") +ylab("Granulicatella relative abundance, %")

#since Granulicatella was positively associated with age too, making a scatter plot
#of age vs. Granulicatella abundance with color by diarrhea. 

a <- ggplot(aes(x=stool_age_days, y=Relative.Abundance), data=(granulicatella.rarefied[granulicatella.rarefied$infant.diarrhea.ever.in.study=="no",]))+geom_line(aes(group=mid),color="grey", data=(granulicatella.rarefied[granulicatella.rarefied$infant.diarrhea.ever.in.study=="no",]))+geom_point(aes(color=infant.diarrhea), data=(granulicatella.rarefied[granulicatella.rarefied$infant.diarrhea.ever.in.study=="no",]))+scale_y_continuous("", lim=c(0.0, 1.15))+scale_color_manual(values=c(no="black", yes="#F8766D"))+scale_x_continuous("age (days)",lim=c(0,260))+labs(color="diarrhea")+ggtitle("A\ndiarrhea never reported")+theme(legend.direction="horizontal")
b <- ggplot(aes(x=stool_age_days, y=Relative.Abundance), data=(granulicatella.rarefied[granulicatella.rarefied$infant.diarrhea.ever.in.study=="yes",]))+geom_line(aes(group=mid), color="grey", data=(granulicatella.rarefied[granulicatella.rarefied$infant.diarrhea.ever.in.study=="yes",]))+geom_point(aes(color=infant.diarrhea), data=(granulicatella.rarefied[granulicatella.rarefied$infant.diarrhea.ever.in.study=="yes",]))+scale_y_continuous("Granulicatella relative abundance %",lim=c(0.0, 1.15))+scale_color_manual(values=c(no="black", yes="#F8766D"))+scale_x_continuous("age (days)",lim=c(0,260))+labs(color="diarrhea")+ggtitle("B\ndiarrhea reported at least once")+theme(legend.direction="horizontal")
ggarrange(a,b, legend.grob = get_legend(b), legend=c("top"))
ggsave("age.vs.granulicatella.jpeg", plot=last_plot(), dpi=600, height=6,width=10, units="in")

#Part 6: determine if Lactobacillus rhamnosus is differentially abundant by probiotic use
#when using rarefied counts and non-parametric testing

phylobj.rarefied.species <- tax_glom(phylobj.rarefied, 'Species', NArm=F)

melted.phylspecies.rarefied <- psmelt(phylobj.rarefied.species)
unique(melted.phylspecies.rarefied$Species)
Lrhamnosus.rarefied <- melted.phylspecies.rarefied[c(grep(pattern="rhamnosus",melted.phylspecies.rarefied$Species)),]
#only look at the samples from visit 2 that were analyzed. 
Lrhamnosus.rarefied.vt2 <- Lrhamnosus.rarefied[c(which(Lrhamnosus.rarefied$Sample %in% metadata.v2.visit2$sample.id)),]
table(Lrhamnosus.rarefied.vt2$Abundance>0) #29 samples had a count above zero. 

wilcox.test(x=(Lrhamnosus.rarefied.vt2[Lrhamnosus.rarefied.vt2$infant.probiotic=="yes","Abundance"]/6963*100), y= (Lrhamnosus.rarefied.vt2[Lrhamnosus.rarefied.vt2$infant.probiotic=="no","Abundance"]/6963*100), paired=F, alternative="greater")
#yes, L.rhamnosus was significantly greater in stool from infants that had probiotics 
#than in stool from infants that didn't have probiotic in visit 2. 
table(Lrhamnosus.rarefied.vt2$infant.probiotic) #24 had probiotic and 85 did not
ggplot(aes(x=infant.probiotic, y=((Abundance)/6963 *100)), data=Lrhamnosus.rarefied.vt2)+geom_violin()+xlab("probiotic use")+ylab("L. rhamnosus relative abundance %")+ggtitle("visit 2 stool samples")+scale_x_discrete(labels=c("no\nn=85", "yes\nn=24"))
ggsave("Lactobacillus.rhamnosus.vs.probiotic.visit2.jpeg", plot=last_plot(), dpi=600, height=4, width=5)


#Part 7: among the 24 infants that received probiotic in visit 2 and were healthy,
#look for any reports of the type of probiotic used and determine if L. rhamnosus 
#was in them. And then do the same for the 24 infants in visit 2 (including those 
#with a morbidity) that had probiotic.  

#read in redcap
redcap <- read.csv("../../../../REDCap_downloads/Denmark/MILQMAINSTUDYDENMARK_DATA_2022-10-11_2251.csv")

#select the probiotic information for the infants of interest. 
redcap.probiotic <- redcap[,c(grep(pattern="^mid|f204_supplst_q10_5|f204_brndty_q10_5_1|f204_brndty_q10_5_1a", colnames(redcap)))]
redcap.probiotic <- redcap.probiotic[c(which(redcap.probiotic$mid %in% (Lrhamnosus.rarefied.vt2[Lrhamnosus.rarefied.vt2$infant.probiotic=="yes",][["mid"]]))),]
length(unique(redcap.probiotic$mid)) #24 infants
table(redcap.probiotic$f204_brndty_q10_5_1) #nineteen knew the brand and one did not
redcap.probiotic$f204_brndty_q10_5_1a
#13 infants had "lactocare" which contains Lactobacillus rhamnosus as well as
#Lactobacillus reuteri, five infants had "duolac" which contains bifidobacteria
#(infantis, breve, longum, bifidum), two infants had "semper (gaia)" and one infant
#had "Biogaia" which both contain Lactobacillus reuteri, and one infant had probiotic
#vitamin D due to stomach pain but no product name provided. 

#Part 8: determine if Bacteroidales and Bacteroidota were negatively associated 
#with fever in visit 3 when using rarefied counts.

#step 1: create a phyloseq object with count summation at the Order level. 
phylobj.rarefied.order <- tax_glom(phylobj.rarefied, 'Order', NArm=F)

melted.phylorder.rarefied <- psmelt(phylobj.rarefied.order)

unique(melted.phylorder.rarefied$Order)
bacteroidales.rarefied <- melted.phylorder.rarefied[c(grep(pattern="Bacteroidales",melted.phylorder.rarefied$Order)),]
bacteroidales.rarefied.visit3 <- bacteroidales.rarefied[bacteroidales.rarefied$visit=="3",]
#step 2: change to relative abundance. 
bacteroidales.rarefied.visit3$relative.abundance <- (bacteroidales.rarefied.visit3$Abundance/6963)*100
table(bacteroidales.rarefied.visit3$Abundance>0) #79 samples had a count above zero. 

#step 3: perform wilcoxon-rank sum test. 
wilcox.test(x=(bacteroidales.rarefied.visit3[bacteroidales.rarefied.visit3$infant.fever=="yes","relative.abundance"]), y= (bacteroidales.rarefied.visit3[bacteroidales.rarefied.visit3$infant.fever=="no","relative.abundance"]), paired=F, alternative="less")
#yes, significantly less Bacteroidales counts among infants with fever in visit 3 versus infants without fever in visit 3. 

#how many infants had fever in visit 3? 
table(bacteroidales.rarefied.visit3$infant.fever)

#step 4: graph and save. 
ggplot(aes(x=infant.fever,y=(relative.abundance)), data=bacteroidales.rarefied.visit3)+ geom_violin()
ggplot(aes(x=infant.fever,y=(relative.abundance), color=as.factor(infant.fever)), data=bacteroidales.rarefied.visit3)+ geom_boxplot()+xlab("fever in visit 3")+ylab("Bacteroidales relative abundance %") +scale_x_discrete(labels=c("no\nn=68", "yes\nn=41"))+scale_color_manual(values=c(no="black", yes="#00BA38"))+theme(legend.position = "none")
ggsave("fever.visit3.vs.Bacteroidales.jpeg", plot=last_plot(), dpi=600, width=10, height=10, units="cm")

#step 5: create a phyloseq object with count summation at the Phylum level
phylobj.rarefied.phylum<- tax_glom(phylobj.rarefied, 'Phylum', NArm=F)

melted.phyl.phylum.rarefied <- psmelt(phylobj.rarefied.phylum)

unique(melted.phyl.phylum.rarefied$Phylum)
bacteroidota.rarefied <- melted.phyl.phylum.rarefied[c(grep(pattern="Bacteroidota",melted.phyl.phylum.rarefied$Phylum)),]
bacteroidota.rarefied.visit3 <- bacteroidota.rarefied[bacteroidota.rarefied$visit=="3",]
table(bacteroidota.rarefied.visit3$Abundance>0) #82 samples had a count above zero. 

bacteroidota.rarefied.visit3$relative.abundance <- (bacteroidota.rarefied.visit3$Abundance/6963)*100
#step 6: run wilcoxon rank sum test to determine if relative abundance (of rarefied counts)
#differed with infant fever in visit 3
wilcox.test(x=(bacteroidota.rarefied.visit3[bacteroidota.rarefied.visit3$infant.fever=="yes","relative.abundance"]), y= (bacteroidota.rarefied.visit3[bacteroidota.rarefied.visit3$infant.fever=="no","relative.abundance"]), paired=F, alternative="less")
#yes, bacteroidota relative abundance was greater among infants that didn't have fever in visit 3
#versus infants that did have fever in visit 3. 
ggplot(aes(x=infant.fever,y=(relative.abundance)), data=bacteroidota.rarefied.visit3)+ geom_violin()
ggplot(aes(x=infant.fever,y=(relative.abundance)), data=bacteroidota.rarefied.visit3)+ geom_boxplot() +xlab("fever in visit 3")+ylab("Bacteroidota relative abundance %") + scale_x_discrete(labels=c("no\nn=68", "yes\nn=41"))
ggsave("fever.visit3.vs.Bacteroidota.jpeg", plot=last_plot(), dpi=600, width=10, height=10, units="cm")


