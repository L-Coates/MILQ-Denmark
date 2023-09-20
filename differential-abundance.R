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

#STEP 1: read in the metadata, (unrarefied) feature table, and the taxonomy table
metadata <- read.csv("../morbidity_medication_prevalences/infant.morbidities.medications.csv", header = T)
metadata <- metadata[,-1]

gender.mod.parity <- read.csv("../sociodemographics_and_sample_collection_info/birthmode.gender.parity.csv", header = T)
gender.mod.parity <- gender.mod.parity[,-1]
gender.mod.parity <- gender.mod.parity[,c(grep(pattern="^mid|visit$|parity|infant[.]gender|mode[.]of[.]delivery|number[.]people[.]in", colnames(gender.mod.parity), ignore.case = T))]
gender.mod.parity <- distinct(gender.mod.parity)

metadata.v2 <- merge(x=metadata, y=gender.mod.parity, by="mid", all=F)
str(metadata.v2)

metadata.v2 <- metadata.v2[,c(2:38,1)]

feature.table <- read_qza("../feature-table-no-chloroplast-eukarya-mitochondria.qza")
feature.table <- feature.table$data

taxonomy <- read_qza("../taxonomy-classification.qza")
taxonomy2 <- as.data.frame(taxonomy$data) %>% column_to_rownames("Feature.ID") %>% separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>% as.matrix()
#remove "Confidence" column
taxonomy2 <- taxonomy2[,-8]

#trim taxonomy to include only ASVs in the ASV table
taxonomy3 <- taxonomy2[c(rownames(taxonomy2)%in%rownames(feature.table)),]

#STEP 2: run ANCOM-BC2
phylobj <- phyloseq(otu_table(feature.table, taxa_are_rows=T),tax_table(taxonomy3), sample_data(metadata.v2 %>% as.data.frame() %>% column_to_rownames("sample.id")))
#get rid of the taxa that are not present at all
phylobj.filtered <- filter_taxa(phylobj, function(x) sum(x)>0, TRUE)
#look at the number of ASVs among the feature table
dim(otu_table(phylobj.filtered))
#there are 1782 ASVs
dim(sample_data(phylobj.filtered)) #327 samples. 


#Running ANCOM-BC2 to determine if diarrhea associates with stool microbial taxa
#and including the covariates suspected of impacting stool microbiota. 
set.seed(123)
ancom.md1.output.ASV <-ancombc2(data=phylobj.filtered, tax_level=NULL, fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+parity+number.people.in.household+infant.diarrhea+infant.fever+infant.vomit", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md1.output.ASV.res <- ancom.md1.output.ASV$res 
lapply(ancom.md1.output.ASV.res[,c(grep(pattern="diff_", colnames(ancom.md1.output.ASV.res)))], FUN=table)
#28 ASVs were different by infant age

#removing the co-variates that did not associate with any ASVs, and running
#a similar model. 
set.seed(123)
ancom.md2.output.ASV <-ancombc2(data=phylobj.filtered, tax_level=NULL, fix_formula="stool_age_days+infant.diarrhea", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.output.ASV.res <- ancom.md2.output.ASV$res 
lapply(ancom.md2.output.ASV.res[,c(grep(pattern="diff_", colnames(ancom.md2.output.ASV.res)))], FUN=table)
#25 ASVs were different by infant age. 

#look at the species level with all co-variates
set.seed(123)
ancom.md1.output.species <-ancombc2(data=phylobj.filtered, tax_level="Species", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+parity+number.people.in.household+infant.diarrhea+infant.fever+infant.vomit", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md1.output.species.res <- ancom.md1.output.species$res 
lapply(ancom.md1.output.species.res[,c(grep(pattern="diff_", colnames(ancom.md1.output.species.res)))], FUN=table)
#14 species differed by age and one by diarrhea -- s__uncultured_bacterium_17. 
print(ancom.md1.output.species.res[ancom.md1.output.species.res$diff_infant.diarrheayes=="TRUE",])
print(ancom.md1.output.species.res[ancom.md1.output.species.res$diff_stool_age_days=="TRUE",])

#make a figure of the log fold change of these species in relation to age in days. 
species.diff.by.age <- ancom.md1.output.species.res[ancom.md1.output.species.res$diff_stool_age_days=="TRUE",]
ggplot(aes(x=reorder(taxon,-lfc_stool_age_days), y=lfc_stool_age_days), data=species.diff.by.age)+geom_col()+ylab("infant age at stool collection (days)")+xlab("")+theme(axis.text.x = element_text(angle = 45, hjust=1))


#cut down model at species level
set.seed(123)
ancom.md2.output.species <-ancombc2(data=phylobj.filtered, tax_level="Species", fix_formula="stool_age_days+infant.diarrhea", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
#"Warning message:Model failed to converge with 1 negative eigenvalue: -1.1e-01" 

#look at the genus level
set.seed(123)
ancom.md1.output.genus <-ancombc2(data=phylobj.filtered, tax_level="Genus", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+parity+number.people.in.household+infant.diarrhea+infant.fever+infant.vomit", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md1.output.genus.res <- ancom.md1.output.genus$res 
lapply(ancom.md1.output.genus.res[,c(grep(pattern="diff_", colnames(ancom.md1.output.genus.res)))], FUN=table)
#15 genera were different by age, and one genus -- Granulicatella -- was different with diarrhea.
print(ancom.md1.output.genus.res[ancom.md1.output.genus.res$diff_infant.diarrheayes=="TRUE",])
print(ancom.md1.output.genus.res[ancom.md1.output.genus.res$diff_stool_age_days=="TRUE",])

#make a figure of the log fold change of these species in relation to age in days. 
genus.diff.by.age <- ancom.md1.output.genus.res[ancom.md1.output.genus.res$diff_stool_age_days=="TRUE",]
genus.diff.by.age$taxon <- gsub(pattern="Genus: g__", replacement="", genus.diff.by.age$taxon) 
#pull out the taxonomy for just the 15 genera different by age. 
taxonomy4 <- data.frame(taxonomy3)
rownames(taxonomy4)=NULL
taxonomy4 <- distinct(taxonomy4)
taxonomy4$Genus <- gsub(pattern=" g__", replacement="", taxonomy4$Genus)
genus.diff.by.age.taxonomy <- taxonomy4[c(which(taxonomy4$Genus %in% genus.diff.by.age$taxon)),]
ggplot(aes(x=reorder(taxon,-lfc_stool_age_days), y=lfc_stool_age_days, fill=lfc_stool_age_days), data=genus.diff.by.age)+geom_col()+ylab("log(fold change) of normalized genus counts with infant age (days)")+xlab("")+theme(axis.text.x = element_text(angle = 45, hjust=1))+scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0)+labs(fill="log(fold change)")
ggsave("stool.genus.logFC.by.age.jpeg", dpi=600, plot=last_plot(), width=10, height=6.5)

#cut down model to remove covariates that were not significantly associated with a genus
set.seed(123)
ancom.md2.output.genus <-ancombc2(data=phylobj.filtered, tax_level="Genus", fix_formula="stool_age_days+infant.diarrhea", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
#"Warning message: Model failed to converge with 1 negative eigenvalue: -1.1e-01"
ancom.md2.output.genus.res <- ancom.md2.output.genus$res 
lapply(ancom.md2.output.genus.res[,c(grep(pattern="diff_", colnames(ancom.md2.output.genus.res)))], FUN=table)
#one genus different by diarrhea, 16 genera different by age. 
print(ancom.md2.output.genus.res[ancom.md2.output.genus.res$diff_infant.diarrheayes=="TRUE",])
#Granulicatella was different by diarrhea. 
#checking that this genus was not sensitive to pseudo count for the condition of interest (infant diarrhea)
pseudo_sens = ancom.md2.output.genus$pseudo_sens_tab

#look at the family level
set.seed(123)
ancom.md1.output.family <-ancombc2(data=phylobj.filtered, tax_level="Family", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+parity+number.people.in.household+infant.diarrhea+infant.fever+infant.vomit", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md1.output.family.res <- ancom.md1.output.family$res 
lapply(ancom.md1.output.family.res[,c(grep(pattern="diff_", colnames(ancom.md1.output.family.res)))], FUN=table)
print(ancom.md1.output.family.res[ancom.md1.output.family.res$diff_infant.diarrheayes=="TRUE",])
print(ancom.md1.output.family.res[ancom.md1.output.family.res$diff_stool_age_days=="TRUE",])
#11 families were different by age, and one family -- Carnobacteriaceae -- was different with diarrhea.

#look at the family level with the cut-down model leaving out the co-variates
#that didn't associate with any family-level microbial taxa
set.seed(123)
ancom.md2.output.family <-ancombc2(data=phylobj.filtered, tax_level="Family", fix_formula="stool_age_days + infant.diarrhea", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md2.output.family.res <- ancom.md2.output.family$res 
lapply(ancom.md2.output.family.res[,c(grep(pattern="diff_", colnames(ancom.md2.output.family.res)))], FUN=table)
print(ancom.md2.output.family.res[ancom.md2.output.family.res$diff_infant.diarrheayes=="TRUE",])
print(ancom.md2.output.family.res[ancom.md2.output.family.res$diff_stool_age_days=="TRUE",])
#14 families were different by age, and one family -- Carnobacteriaceae -- was different with diarrhea.

#look at the order level
set.seed(123)
ancom.md1.output.order <-ancombc2(data=phylobj.filtered, tax_level="Order", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+parity+number.people.in.household+infant.diarrhea+infant.fever+infant.vomit", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.md1.output.order.res <- ancom.md1.output.order$res 
lapply(ancom.md1.output.order.res[,c(grep(pattern="diff_", colnames(ancom.md1.output.order.res)))], FUN=table)
print(ancom.md1.output.order.res[ancom.md1.output.order.res$diff_infant.diarrheayes=="TRUE",])
print(ancom.md1.output.order.res[ancom.md1.output.order.res$diff_stool_age_days=="TRUE",])
#9 orders were different by age, and none were different with diarrhea or other co-variates

#look at the order level with the cut-down model leaving out the co-variates
#that didn't associate with any order-level microbial taxa
set.seed(123)
ancom.md2.output.order <-ancombc2(data=phylobj.filtered, tax_level="Order", fix_formula="stool_age_days + infant.diarrhea", rand_formula = "(1|mid)", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
#"warning message: model failed to converge with 1 negative eigenvalue: -1.1e-01"
ancom.md2.output.order.res <- ancom.md2.output.order$res 
lapply(ancom.md2.output.order.res[,c(grep(pattern="diff_", colnames(ancom.md2.output.order.res)))], FUN=table)
print(ancom.md2.output.order.res[ancom.md2.output.order.res$diff_infant.diarrheayes=="TRUE",])
print(ancom.md2.output.order.res[ancom.md2.output.order.res$diff_stool_age_days=="TRUE",])
#9 orders were different with age, and none were different with diarrhea


#STEP 3: testing the hypothesis that the stool microbiota in visit 2 is predictive
#of diarrhea later on in visits 3 and/or 4. 

#select just visit 2 stool samples
visit2.samples <- metadata.v2[metadata.v2$visit=="2",]
#remove infants that had diarrhea in that visit 2
visit2.samples <- visit2.samples[visit2.samples$infant.diarrhea=="no",] 

visit2.samples.list <- visit2.samples$sample.id
phylobj.visit2 <- prune_samples(samples=visit2.samples.list, x=phylobj)

#get rid of the taxa that are not present at all
phylobj.visit2.filtered <- filter_taxa(phylobj.visit2, function(x) sum(x)>0, TRUE)
#look at the number of ASVs among the feature table
dim(otu_table(phylobj.visit2.filtered))
#there are 775 ASVs, 87 samples

#use ANCOM-BC2 to determine if any taxa in visit 2 are associated
#with diarrhea later on in visits 3 or 4.

#full model for ASVs 
set.seed(123)
ancom.vt2.md1.output.ASV <-ancombc2(data=phylobj.visit2.filtered, tax_level=NULL, fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+parity+number.people.in.household+infant.vomit+infant.fever+infant.diarrhea.ever.in.study", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.vt2.md1.output.ASV.res <- ancom.vt2.md1.output.ASV$res 
lapply(ancom.vt2.md1.output.ASV.res[,c(grep(pattern="diff_", colnames(ancom.vt2.md1.output.ASV.res)))], FUN=table)
print(ancom.vt2.md1.output.ASV.res[ancom.vt2.md1.output.ASV.res$diff_mode.of.deliveryvaginal=="TRUE",])
print(ancom.vt2.md1.output.ASV.res[ancom.vt2.md1.output.ASV.res$diff_stool_age_days=="TRUE",])

#one ASV, 7e5340dea51bfd59b29795487da71dd6, was negatively associated with vaginal birth
#and it belongs to an unknown taxa within Enterobacteriaceae. 
#and one ASV, e45861c23b60134373a1cd62ea569b80, was negatively associated with age
#and it belongs to Parabacteroides distasonis 

#cut-down model for ASVs which still includes age, vaginal birth since those had
#significant associations with stool microbes

set.seed(123)
ancom.vt2.md2.output.ASV <-ancombc2(data=phylobj.visit2.filtered, tax_level=NULL, fix_formula="stool_age_days +mode.of.delivery+infant.diarrhea.ever.in.study", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.vt2.md2.output.ASV.res <- ancom.vt2.md2.output.ASV$res 
lapply(ancom.vt2.md2.output.ASV.res[,c(grep(pattern="diff_", colnames(ancom.vt2.md2.output.ASV.res)))], FUN=table)
#still no ASVs in visit 2 were associated with diarrhea later in visits 3,4. 

#full model for species
set.seed(123)
ancom.vt2.md1.output.species <-ancombc2(data=phylobj.visit2.filtered, tax_level="Species", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+parity+number.people.in.household+infant.vomit+infant.fever+infant.diarrhea.ever.in.study", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.vt2.md1.output.species.res <- ancom.vt2.md1.output.species$res 
lapply(ancom.vt2.md1.output.species.res[,c(grep(pattern="diff_", colnames(ancom.vt2.md1.output.species.res)))], FUN=table)
print(ancom.vt2.md1.output.species.res[ancom.vt2.md1.output.species.res$diff_infant.probioticyes=="TRUE",])
print(ancom.vt2.md1.output.species.res[ancom.vt2.md1.output.species.res$diff_infant.laxativeyes=="TRUE",])
#one species, Lactobacillus rhamnosus, was different with probiotic use in visit 2
#and it was a positive association.
#one species, Peptostreptococcus anaerobius, was positively associated with laxative use, 
#but there was only one infant that received laxative in visit 2. 
pseudo_sens = ancom.vt2.md1.output.species$pseudo_sens_tab
#not sensitive to pseudo count

#cut-down model for species (keeping only the covariates: probiotic use and laxative use)
set.seed(123)
ancom.vt2.md2.output.species <-ancombc2(data=phylobj.visit2.filtered, tax_level="Species", fix_formula="infant.probiotic+infant.laxative+infant.diarrhea.ever.in.study", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.vt2.md2.output.species.res <- ancom.vt2.md2.output.species$res 
lapply(ancom.vt2.md2.output.species.res[,c(grep(pattern="diff_", colnames(ancom.vt2.md2.output.species.res)))], FUN=table)
#still no differentially abundant species with diarrhea. 

#full model for genus
set.seed(123)
ancom.vt2.md1.output.genus <-ancombc2(data=phylobj.visit2.filtered, tax_level="Genus", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+parity+number.people.in.household+infant.vomit+infant.fever+infant.diarrhea.ever.in.study", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.vt2.md1.output.genus.res <- ancom.vt2.md1.output.genus$res 
lapply(ancom.vt2.md1.output.genus.res[,c(grep(pattern="diff_", colnames(ancom.vt2.md1.output.genus.res)))], FUN=table)
print(ancom.vt2.md1.output.genus.res[ancom.vt2.md1.output.genus.res$diff_infant.laxativeyes=="TRUE",])
#four genera were different with laxative but only one infant had laxative in visit 2. 


#cut-down model for genus level taxa, and only including laxative use as covariate
set.seed(123)
ancom.vt2.md2.output.genus <-ancombc2(data=phylobj.visit2.filtered, tax_level="Genus", fix_formula="infant.laxative+infant.diarrhea.ever.in.study", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.vt2.md2.output.genus.res <- ancom.vt2.md2.output.genus$res 
lapply(ancom.vt2.md2.output.genus.res[,c(grep(pattern="diff_", colnames(ancom.vt2.md2.output.genus.res)))], FUN=table)
#still no association between diarrhea and stool microbial genera.

#look at the family level now 
set.seed(123)
ancom.vt2.md1.output.family <-ancombc2(data=phylobj.visit2.filtered, tax_level="Family", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+parity+number.people.in.household+infant.vomit+infant.fever+infant.diarrhea.ever.in.study", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.vt2.md1.output.family.res <- ancom.vt2.md1.output.family$res 
lapply(ancom.vt2.md1.output.family.res[,c(grep(pattern="diff_", colnames(ancom.vt2.md1.output.family.res)))], FUN=table)
#two family level taxa were differentially abundant with laxative use, but none 
#were associated with diarrhea. 

#cut down model at family level with laxative use retained as covariate.  
set.seed(123)
ancom.vt2.md2.output.family <-ancombc2(data=phylobj.visit2.filtered, tax_level="Family", fix_formula="infant.laxative+infant.diarrhea.ever.in.study", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.vt2.md2.output.family.res <- ancom.vt2.md2.output.family$res 
lapply(ancom.vt2.md2.output.family.res[,c(grep(pattern="diff_", colnames(ancom.vt2.md2.output.family.res)))], FUN=table)
#still no families differentially abundant with diarrhea. 


#look at the order level now 
set.seed(123)
ancom.vt2.md1.output.order <-ancombc2(data=phylobj.visit2.filtered, tax_level="Order", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+parity+number.people.in.household+infant.vomit+infant.fever+infant.diarrhea.ever.in.study", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.vt2.md1.output.order.res <- ancom.vt2.md1.output.order$res 
lapply(ancom.vt2.md1.output.order.res[,c(grep(pattern="diff_", colnames(ancom.vt2.md1.output.order.res)))], FUN=table)
#one order level taxon was differentially abundant with laxative use, but none 
#were associated with diarrhea. 

#cut-down model at the order level, and only include laxative use as covariates 
set.seed(123)
ancom.vt2.md2.output.order <-ancombc2(data=phylobj.visit2.filtered, tax_level="Order", fix_formula="infant.laxative+infant.diarrhea.ever.in.study", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.vt2.md2.output.order.res <- ancom.vt2.md2.output.order$res 
lapply(ancom.vt2.md2.output.order.res[,c(grep(pattern="diff_", colnames(ancom.vt2.md2.output.order.res)))], FUN=table)
#still no orders were associated with diarrhea. 


#STEP 4: determine if Granulicatella is differentially abundant by diarrhea
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

ggplot(aes(x=infant.diarrhea, y=Relative.Abundance), data=granulicatella.rarefied)+geom_boxplot()+xlab("diarrhea")+ylab("Granulicatella relative abundance %")+scale_x_discrete(labels=c("no\nn=287", "yes\nn=40"))
ggsave("granulicatella.vs.diarrhea.jpeg", dpi=600, width=5, height=5, plot=last_plot())

#also do the graph with non-rarefied counts but with relative abundance
phylobj.genus <- tax_glom(phylobj.filtered, 'Genus', NArm=F)

#Create a list of unique taxa names that can be assigned to the features within the phyloseq object
taxa.names <- tax_table(phylobj.genus)[,6]

#replace empty names with family level identities
taxa.names2 <- tax_table(phylobj.genus)[,5] 
taxa.names[c(which(is.na(taxa.names)), grep("g__$", taxa.names))] <-taxa.names2[c(which(is.na(taxa.names)),grep("g__$", taxa.names))]

#replace empty names with order level identities
taxa.names3 <-tax_table(phylobj.genus)[,4] 
taxa.names[c(which(is.na(taxa.names2)), grep( "g__$", taxa.names))] <- taxa.names3[c(which(is.na(taxa.names2)), grep("g__$", taxa.names))]

#replace empty names with class level identities
taxa.names4 <- tax_table(phylobj.genus)[,3]
taxa.names[c(which(is.na(taxa.names3)), grep("g__$", taxa.names))] <-taxa.names4[c(which(is.na(taxa.names3)),grep("g__$", taxa.names))]

#replace empty names with phylum level identities
taxa.names5 <- tax_table(phylobj.genus)[,2] 
taxa.names[c(which(is.na(taxa.names4)),grep("g__$", taxa.names))] <-taxa.names5[c(which(is.na(taxa.names4)),grep("g__$", taxa.names))]

#replace empty names with kingdom level identities
taxa.names6 <- tax_table(phylobj.genus)[,1] 
taxa.names[c(which(is.na(taxa.names5)),grep("g__$", taxa.names))] <-taxa.names6[c(which(is.na(taxa.names5)),grep("g__$", taxa.names))]

#remove leading white space
taxa.names <- gsub(" ","", taxa.names)

#change the taxa designations to names compatible with R
taxa.names <- make.names(taxa.names, unique = TRUE)

#Make all of the feature names within the phyloseq object correspond to the genus level designation
taxa_names(phylobj.genus) <- taxa.names

phylobj.genus.relabund <- transform_sample_counts(phylobj.genus, function(x) x/sum(x) * 100)

melted.phylgenus.relabund <- psmelt(phylobj.genus.relabund)
unique(melted.phylgenus.relabund$Genus)
granulicatella <- melted.phylgenus.relabund[c(grep(pattern="Granulica",melted.phylgenus.relabund$Genus)),]
table(granulicatella$Abundance>0) #65 samples had a count above zero. 

ggplot(aes(x=infant.diarrhea, y=Abundance), data=granulicatella)+geom_boxplot()+xlab("did infant have diarrhea during visit?") +ylab("Granulicatella relative abundance, %")

#since Granulicatella was positively associated with age too, making a scatter plot
#of age vs. Granulicatella abundance with color by diarrhea. 

b <- ggplot(aes(x=stool_age_days, y=Relative.Abundance), data=(granulicatella.rarefied[granulicatella.rarefied$infant.diarrhea.ever.in.study=="no",]))+geom_line(aes(group=mid), data=(granulicatella.rarefied[granulicatella.rarefied$infant.diarrhea.ever.in.study=="no",]))+geom_point(aes(color=infant.diarrhea), data=(granulicatella.rarefied[granulicatella.rarefied$infant.diarrhea.ever.in.study=="no",]))+scale_y_continuous("", lim=c(0.0, 1.15))+scale_x_continuous("age (days)",lim=c(0,260))+labs(color="diarrhea")+ggtitle("B\ndiarrhea never reported")
a <- ggplot(aes(x=stool_age_days, y=Relative.Abundance), data=(granulicatella.rarefied[granulicatella.rarefied$infant.diarrhea.ever.in.study=="yes",]))+geom_line(aes(group=mid), data=(granulicatella.rarefied[granulicatella.rarefied$infant.diarrhea.ever.in.study=="yes",]))+geom_point(aes(color=infant.diarrhea), data=(granulicatella.rarefied[granulicatella.rarefied$infant.diarrhea.ever.in.study=="yes",]))+scale_y_continuous("Granulicatella relative abundance %",lim=c(0.0, 1.15))+scale_x_continuous("age (days)",lim=c(0,260))+labs(color="diarrhea")+ggtitle("A\ndiarrhea reported at least once")
ggarrange(a,b, common.legend = TRUE)
ggsave("age.vs.granulicatella.jpeg", plot=last_plot(), dpi=600, height=6,width=10, units="in")

#STEP 5: determine if Lactobacillus rhamnosus is differentially abundant by probiotic use
#when using rarefied counts and non-parametric testing

phylobj.rarefied.species <- tax_glom(phylobj.rarefied, 'Species', NArm=F)

melted.phylspecies.rarefied <- psmelt(phylobj.rarefied.species)
unique(melted.phylspecies.rarefied$Species)
Lrhamnosus.rarefied <- melted.phylspecies.rarefied[c(grep(pattern="rhamnosus",melted.phylspecies.rarefied$Species)),]
#only look at the samples from visit 2 that were analyzed. 
Lrhamnosus.rarefied.vt2 <- Lrhamnosus.rarefied[c(which(Lrhamnosus.rarefied$Sample %in% visit2.samples.list)),]
table(Lrhamnosus.rarefied.vt2$Abundance>0) #26 samples had a count above zero. 

wilcox.test(x=(Lrhamnosus.rarefied.vt2[Lrhamnosus.rarefied.vt2$infant.probiotic=="yes","Abundance"]/6963*100), y= (Lrhamnosus.rarefied.vt2[Lrhamnosus.rarefied.vt2$infant.probiotic=="no","Abundance"]/6963*100), paired=F, alternative="greater")
#yes, L.rhamnosus was significantly greater in stool from infants that had probiotics 
#than in stool from infants that didn't have probiotic in visit 2. 
table(Lrhamnosus.rarefied.vt2$infant.probiotic) #20 had probiotic and 81 did not
ggplot(aes(x=infant.probiotic, y=((Abundance)/6963 *100)), data=Lrhamnosus.rarefied.vt2)+geom_boxplot()+xlab("probiotic use")+ylab("L. rhamnosus relative abundance %")+ggtitle("visit 2 stool samples")+scale_x_discrete(labels=c("no\nn=81", "yes\nn=20"))
ggsave("Lactobacillus.rhamnosus.vs.probiotic.visit2.jpeg", plot=last_plot(), dpi=600, height=4, width=5)

Lrhamnosus.rarefied.vt2.all <- Lrhamnosus.rarefied[Lrhamnosus.rarefied$visit=="2",]
table(Lrhamnosus.rarefied.vt2.all$infant.probiotic) #24 had probiotic and 85 did not
table(Lrhamnosus.rarefied.vt2.all$Abundance>0) #29 samples had a count above zero. 

wilcox.test(x=(Lrhamnosus.rarefied.vt2.all[Lrhamnosus.rarefied.vt2.all$infant.probiotic=="yes","Abundance"]/6963*100), y= (Lrhamnosus.rarefied.vt2.all[Lrhamnosus.rarefied.vt2.all$infant.probiotic=="no","Abundance"]/6963*100), paired=F, alternative="greater")
#yes, when we add back in the infants that had diarrhea in visit 2, we still
#find that L.rhamnosus was significantly greater in stool from infants that had 
#probiotics than in stool from infants that didn't have probiotic in visit 2. 

#seeing if ANCOM-BC2 produces the same result
#select all visit 2 stool samples this time
phylobj.visit2.all <- prune_samples(samples=(metadata.v2[metadata.v2$visit=="2",][["sample.id"]]), x=phylobj)

#get rid of the taxa that are not present at all
phylobj.visit2.all.filtered <- filter_taxa(phylobj.visit2.all, function(x) sum(x)>0, TRUE)
#look at the number of ASVs among the feature table
dim(otu_table(phylobj.visit2.all.filtered))
#there are 794 ASVs, 109 samples

set.seed(123)
ancom.vt2.all.species <-ancombc2(data=phylobj.visit2.all.filtered, tax_level="Species", fix_formula="stool_age_days + infant.antibiotic+ infant.probiotic+infant.laxative+infant.gender + mode.of.delivery+parity+number.people.in.household+infant.fever+infant.diarrhea+infant.vomit", pseudo=0, pseudo_sens=TRUE, p_adj_method="BH")
ancom.vt2.all.species.res <- ancom.vt2.all.species$res 
lapply(ancom.vt2.all.species.res[,c(grep(pattern="diff_", colnames(ancom.vt2.all.species.res)))], FUN=table)
#by ANCOM-BC2, the relationship between L.rhamnosus and probiotic use among all 109 infants in visit 2 was no longer significant. 

table(Lrhamnosus.rarefied.vt2.all$infant.probiotic) #24 had probiotic and 85 did not. 
ggplot(aes(x=infant.probiotic, y=((Abundance)/6963 *100)), data=Lrhamnosus.rarefied.vt2.all)+geom_boxplot()+xlab("probiotic use")+ylab("L. rhamnosus relative abundance %")+ggtitle("visit 2 stool samples")+scale_x_discrete(labels=c("no\nn=85", "yes\nn=24"))
ggsave("Lactobacillus.rhamnosus.vs.probiotic.visit2.allsamples.jpeg", plot=last_plot(), dpi=600, height=4, width=5)

#STEP 6: among the 16 infants that received probiotic in visit 2 and were healthy,
#look for any reports of the type of probiotic used and determine if L. rhamnosus 
#was in them. And then do the same for the 24 infants in visit 2 (including those 
#with a morbidity) that had probiotic.  

#read in redcap
redcap <- read.csv("../../../../REDCap_downloads/Denmark/MILQMAINSTUDYDENMARK_DATA_2022-10-11_2251.csv")

#select the probiotic information for the infants of interest. 
redcap.probiotic <- redcap[,c(grep(pattern="^mid|f204_supplst_q10_5|f204_brndty_q10_5_1|f204_brndty_q10_5_1a", colnames(redcap)))]
redcap.probiotic <- redcap.probiotic[c(which(redcap.probiotic$mid %in% (Lrhamnosus.rarefied.vt2[Lrhamnosus.rarefied.vt2$infant.probiotic=="yes",][["mid"]]))),]
length(unique(redcap.probiotic$mid)) #20 infants
table(redcap.probiotic$f204_brndty_q10_5_1) #nineteen knew the brand and one did not
redcap.probiotic$f204_brndty_q10_5_1a
#11 infants had "lactocare" which contains Lactobacillus rhamnosus as well as
#Lactobacillus reuteri, five infants had "duolac" which contains bifidobacteria
#(infantis, breve, longum, bifidum), one infant had "semper gaia" and one infant
#had "Biogaia" which both contain Lactobacillus reuteri, and one infant had probiotic
#vitamin D due to stomach pain but no product name provided. 

#now look among all 24 infants in visit 2 that received probiotic. 
redcap.probiotic.all <- redcap[,c(grep(pattern="^mid|f204_supplst_q10_5|f204_brndty_q10_5_1|f204_brndty_q10_5_1a", colnames(redcap)))]
redcap.probiotic.all <- redcap.probiotic.all[c(which(redcap.probiotic.all$mid %in% (Lrhamnosus.rarefied.vt2.all[Lrhamnosus.rarefied.vt2.all$infant.probiotic=="yes",][["mid"]]))),]
length(unique(redcap.probiotic.all$mid)) #24 infants
table(redcap.probiotic.all$f204_brndty_q10_5_1) #23 knew the brand and one did not
table(redcap.probiotic.all$f204_brndty_q10_5_1a)
#13 infants receive lactocare, four received semper biogaia, five received duolac, 
#and two infants had some unknown probiotic. 
