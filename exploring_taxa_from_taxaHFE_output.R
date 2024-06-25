#this file contains the code used to explore the microbial taxa from visit 2 that were identified
#with "taxaHFE 2.0" as potentially associated with a morbidity in visits 3, 4. 

#loading libraries
library(qiime2R)
library(tidyverse)
library(phyloseq)
library(caret)
library(MASS)
library(ranger)
library(fastshap)
library(shapviz)
library(ROCR)
library(vip)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(cowplot)
library(ggforce)



#Part 1: create phyloseq objects for each taxonomic level that appeared as potentially
#associated with a morbidity (by taxaHFE 2.0). And also see if age at stool collection
#in visit 2 differed by morbidity development in visits 3 & 4. 

#Step 1: read in the taxa HFE tables for the taxa potentially associated with a morbidity.
taxaHFE.output.diarrhea.visit3and4 <- read.csv("taxaHFE.output.diarrhea.visit3visit4.csv")
colnames(taxaHFE.output.diarrhea.visit3and4)

taxaHFE.output.fever.visit3and4 <- read.csv("taxaHFE.output.fever.visit3visit4.csv")
colnames(taxaHFE.output.fever.visit3and4)

taxaHFE.output.vomit.visit3and4 <- read.csv("taxaHFE.output.vomit.visit3visit4.csv")
colnames(taxaHFE.output.vomit.visit3and4)

#Step 2: reading in the rarefied feature file, taxonomy, and metadata file and
#creating a phyloseq object (for all 327 samples).
feature.table.rarefied <- read_qza("../feature-table-rarefied-6963.qza")
feature.table.rarefied <- feature.table.rarefied$data

taxonomy <- read_qza("../taxonomy-classification.qza")
taxonomy2 <- as.data.frame(taxonomy$data) %>% column_to_rownames("Feature.ID") %>% separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>% as.matrix()
#remove "Confidence" column
taxonomy2 <- taxonomy2[,-8]

#trim taxonomy to include only ASVs in the ASV table
taxonomy3 <- taxonomy2[c(rownames(taxonomy2)%in%rownames(feature.table.rarefied)),]

metadata.part1 <- read.csv("../morbidity_medication_prevalences/infant.morbidities.medications.csv")
metadata.part1 <- metadata.part1[,-1]

metadata.part2 <- read.csv("../sociodemographics_and_sample_collection_info/birthmode.gender.parity.csv")

metadata <- merge(x=metadata.part1, y=metadata.part2[,c(1,20:25)], by="sample.id", all=F)

#step 3: determine if there is a statistically significant difference in age at 
#stool collection in visit 2 with morbidity development in visits 3, 4. And similarly,
#if there is a statistically significant difference in days at room temperature
#for stool in visit 2 with morbidity development in visits 3, 4. 
metadata.visit2 <- metadata[metadata$visit==2,]
ggplot(aes(x=infant.diarrhea.in.visit3.or.visit4, y=stool_age_days), data=metadata.visit2)+geom_boxplot() #the groups look very similar
ggplot(aes(x=infant.fever.in.visit3.or.visit4, y=stool_age_days), data=metadata.visit2)+geom_boxplot()#the groups look very similar
ggplot(aes(x=infant.vomit.in.visit3.or.visit4, y=stool_age_days), data=metadata.visit2)+geom_boxplot()#the groups look very similar
shapiro.test(metadata.visit2$stool_age_days) #not normal distribution
wilcox.test(x=metadata.visit2[metadata.visit2$infant.diarrhea.in.visit3.or.visit4=="yes","stool_age_days"], y=metadata.visit2[metadata.visit2$infant.diarrhea.in.visit3.or.visit4=="no", "stool_age_days"]) #not different
wilcox.test(x=metadata.visit2[metadata.visit2$infant.fever.in.visit3.or.visit4=="yes","stool_age_days"], y=metadata.visit2[metadata.visit2$infant.fever.in.visit3.or.visit4=="no", "stool_age_days"]) #not different
wilcox.test(x=metadata.visit2[metadata.visit2$infant.vomit.in.visit3.or.visit4=="yes","stool_age_days"], y=metadata.visit2[metadata.visit2$infant.vomit.in.visit3.or.visit4=="no", "stool_age_days"]) #not different
#can we normalize the distribution and use the t-test?
bestNormalize::bestNormalize(metadata.visit2$stool_age_days)
metadata.visit2$stool_age_days_ordNorm_transformed <-bestNormalize::bestNormalize(metadata.visit2$stool_age_days)$x.t
shapiro.test(metadata.visit2$stool_age_days_ordNorm_transformed)
t.test(x=metadata.visit2[metadata.visit2$infant.diarrhea.in.visit3.or.visit4=="yes","stool_age_days_ordNorm_transformed"], y=metadata.visit2[metadata.visit2$infant.diarrhea.in.visit3.or.visit4=="no", "stool_age_days_ordNorm_transformed"]) #not different
t.test(x=metadata.visit2[metadata.visit2$infant.fever.in.visit3.or.visit4=="yes","stool_age_days_ordNorm_transformed"], y=metadata.visit2[metadata.visit2$infant.fever.in.visit3.or.visit4=="no", "stool_age_days_ordNorm_transformed"]) #not different
t.test(x=metadata.visit2[metadata.visit2$infant.vomit.in.visit3.or.visit4=="yes","stool_age_days_ordNorm_transformed"], y=metadata.visit2[metadata.visit2$infant.vomit.in.visit3.or.visit4=="no", "stool_age_days_ordNorm_transformed"]) #not different

ggplot(aes(x=infant.diarrhea.in.visit3.or.visit4, y=stool_days_at_room_temperature), data=metadata.visit2)+geom_boxplot() #the groups look very similar
ggplot(aes(x=infant.fever.in.visit3.or.visit4, y=stool_days_at_room_temperature), data=metadata.visit2)+geom_boxplot() #the groups look very similar
ggplot(aes(x=infant.vomit.in.visit3.or.visit4, y=stool_days_at_room_temperature), data=metadata.visit2)+geom_boxplot() #the groups look very similar
shapiro.test(metadata.visit2$stool_days_at_room_temperature) #not normal distribution
wilcox.test(x=metadata.visit2[metadata.visit2$infant.diarrhea.in.visit3.or.visit4=="yes","stool_days_at_room_temperature"], y=metadata.visit2[metadata.visit2$infant.diarrhea.in.visit3.or.visit4=="no", "stool_days_at_room_temperature"]) #not different
wilcox.test(x=metadata.visit2[metadata.visit2$infant.fever.in.visit3.or.visit4=="yes","stool_days_at_room_temperature"], y=metadata.visit2[metadata.visit2$infant.fever.in.visit3.or.visit4=="no", "stool_days_at_room_temperature"]) #not different
wilcox.test(x=metadata.visit2[metadata.visit2$infant.vomit.in.visit3.or.visit4=="yes","stool_days_at_room_temperature"], y=metadata.visit2[metadata.visit2$infant.vomit.in.visit3.or.visit4=="no", "stool_days_at_room_temperature"]) #not different
bestNormalize::bestNormalize(metadata.visit2$stool_days_at_room_temperature)
metadata.visit2$stool_days_at_room_temperature_log_transformed <-bestNormalize::bestNormalize(metadata.visit2$stool_days_at_room_temperature)$x.t
shapiro.test(metadata.visit2$stool_days_at_room_temperature_log_transformed) #still non-normal distribution

#make phyloseq object now
phylobj <- phyloseq(otu_table(feature.table.rarefied, taxa_are_rows=T),tax_table(taxonomy3), sample_data(metadata %>% as.data.frame() %>% column_to_rownames("sample.id")))

#Step 4: select only visit 2 samples and remove taxa that aren't present in any of the samples
phylobj.visit2 <- prune_samples(x=phylobj, samples=c(metadata[metadata$visit=="2","sample.id"]))
phylobj.visit2 <- filter_taxa(phylobj.visit2, function(x) sum(x)>0, TRUE)

#Step 5: gather counts at species level and transform to relative abundance
phylobj.species <- tax_glom(phylobj.visit2, "Species", NArm=F)
phylobj.species.relabund <- transform_sample_counts(phylobj.species, function(x) x/sum(x) * 100)
melted.species <- psmelt(phylobj.species.relabund)

#Step 6: gather counts at genus level and transform to relative abundance
phylobj.genus <- tax_glom(phylobj.visit2, "Genus", NArm=F)
phylobj.genus.relabund <- transform_sample_counts(phylobj.genus, function(x) x/sum(x) * 100)
melted.genus <- psmelt(phylobj.genus.relabund)

#Step 7: gather counts at family level and transform to relative abundance
phylobj.family <- tax_glom(phylobj.visit2, "Family", NArm=F)
phylobj.family.relabund <- transform_sample_counts(phylobj.family, function(x) x/sum(x) * 100)
melted.family <- psmelt(phylobj.family.relabund)

#Step 8: gather counts at order level and transform to relative abundance
phylobj.order <- tax_glom(phylobj.visit2, "Order", NArm=F)
phylobj.order.relabund <- transform_sample_counts(phylobj.order, function(x) x/sum(x) * 100)
melted.order <- psmelt(phylobj.order.relabund)

#Step 9: gather counts at class level and transform to relative abundance
phylobj.class <- tax_glom(phylobj.visit2, "Class", NArm=F)
phylobj.class.relabund <- transform_sample_counts(phylobj.class, function(x) x/sum(x) * 100)
melted.class <- psmelt(phylobj.class.relabund)

#Step 10: gather counts at family level and transform to relative abundance
phylobj.phylum <- tax_glom(phylobj.visit2, "Phylum", NArm=F)
phylobj.phylum.relabund <- transform_sample_counts(phylobj.phylum, function(x) x/sum(x) * 100)
melted.phylum <- psmelt(phylobj.phylum.relabund)

#Step 11: combine the taxonomic abundance data frames. 

#The species names will be a concatenation of the other taxonomic levels since 
#some species names are not unique in the data set and since the species names
#from the taxaHFE output are the concatenation of the taxonomy. 
melted.species$Species <- apply(X=melted.species[,c(grep(pattern="^Kingdom|^Phylum|^Class|^Order|^Family|^Genus|^Species", colnames(melted.species)))],MARGIN=1, function(x) paste(x, collapse="_"))
#confirm that all species names are unique now that they have been combined with 
#genus names. 
table(melted.species$Species)[table(melted.species$Species)>109]

melted.species$taxon <- melted.species$Species
melted.species <- melted.species[,-c(grep(pattern="^OTU|^Kingdom|^Phylum|^Class|^Order|^Family|^Genus|^Species", colnames(melted.species)))]

melted.genus$taxon <- melted.genus$Genus
melted.genus <- melted.genus[,-c(grep(pattern="^OTU|^Kingdom|^Phylum|^Class|^Order|^Family|^Genus", colnames(melted.genus)))]

melted.family$taxon <- melted.family$Family
melted.family <- melted.family[,-c(grep(pattern="^OTU|^Kingdom|^Phylum|^Class|^Order|^Family", colnames(melted.family)))]

melted.order$taxon <- melted.order$Order
melted.order <- melted.order[,-c(grep(pattern="^OTU|^Kingdom|^Phylum|^Class|^Order", colnames(melted.order)))]

melted.class$taxon <- melted.class$Class
melted.class <- melted.class[,-c(grep(pattern="^OTU|^Kingdom|^Phylum|^Class", colnames(melted.class)))]

melted.phylum$taxon <- melted.phylum$Phylum
melted.phylum <- melted.phylum[,-c(grep(pattern="^OTU|^Kingdom|^Phylum", colnames(melted.phylum)))]

melted.alltaxa <- rbind(melted.species, melted.genus, melted.family, melted.order, melted.class, melted.phylum)
melted.alltaxa <- melted.alltaxa[!is.na(melted.alltaxa$taxon),]
melted.alltaxa <- distinct(melted.alltaxa)

#Part 2: Among the taxa flagged as potentially associated with diarrhea in visit 3,4,
#look for associations using random forest.  

#step 1: select taxa flagged by taxaHFE 2.0 for diarrhea.
colnames(taxaHFE.output.diarrhea.visit3and4)
taxaHFE.output.diarrhea.taxa <- colnames(taxaHFE.output.diarrhea.visit3and4[,3:11])
print(taxaHFE.output.diarrhea.taxa)
#all of these species names are unique 
taxaHFE.output.diarrhea.taxa <- "staphylococcales$|bacteroides_vulgatus$|bacteroides_uniformis$|prevotella_timonensis$|flavobacteriales$|corynebacteriales$|sutterella_wadsworthensis$|rhodoferax$|haemophilus$"
diarrhea.df <- melted.alltaxa[c(grep(pattern=taxaHFE.output.diarrhea.taxa, melted.alltaxa$taxon, ignore.case = T)),]
unique(diarrhea.df$taxon)

#step 2: fit random forest model with these taxa. 

#change data frame into wide format
diarrhea.df <- pivot_wider(data=diarrhea.df, names_from="taxon", values_from="Abundance")

#remove leading space at beginning of taxon name
colnames(diarrhea.df) <- gsub(pattern=" ", replacement="", colnames(diarrhea.df))

#check the class type of feature
str(diarrhea.df) #all taxa were numeric and all categorical were character

#change the values for "infant.diarrhea.in.visit3.or.visit4" from "yes" and "no"
#to 1 and 0, respectively. 
diarrhea.df$infant.diarrhea.in.visit3.or.visit4 <- gsub(pattern="no", replacement="0", ignore.case = T, diarrhea.df$infant.diarrhea.in.visit3.or.visit4)
diarrhea.df$infant.diarrhea.in.visit3.or.visit4 <- gsub(pattern="yes", replacement="1", ignore.case = T, diarrhea.df$infant.diarrhea.in.visit3.or.visit4)
diarrhea.df$infant.diarrhea.in.visit3.or.visit4 <- as.factor(diarrhea.df$infant.diarrhea.in.visit3.or.visit4)

#determine if any of the features (aka predictors) are highly correlated. This matters because if 
#any of the features are highly correlated then the feature importance values/shapley values
#will not be accurate for the highly correlated features. 
ggplot(aes(x=infant.diarrhea, y=household.size), data=diarrhea.df)+geom_jitter()
#change the features to numeric values to determine if they are
#correlated with each other
diarrhea.df$infant.diarrhea.numeric <- diarrhea.df$infant.diarrhea
diarrhea.df[diarrhea.df$infant.diarrhea.numeric=="yes", "infant.diarrhea.numeric"] <- "1"
diarrhea.df[diarrhea.df$infant.diarrhea.numeric=="no", "infant.diarrhea.numeric"] <- "0"
diarrhea.df$infant.diarrhea.numeric <- as.numeric(diarrhea.df$infant.diarrhea.numeric)

diarrhea.df$household.size.numeric <- diarrhea.df$household.size
diarrhea.df[diarrhea.df$household.size.numeric=="> 3", "household.size.numeric"] <- "1"
diarrhea.df[diarrhea.df$household.size.numeric=="<= 3", "household.size.numeric"] <- "0"
diarrhea.df$household.size.numeric <- as.numeric(diarrhea.df$household.size.numeric)

corr.matrix.diarrhea.df <- cor(x=diarrhea.df[,c(grep(pattern="household.size.numeric|infant.diarrhea.numeric$|^[sgfocpd]_", colnames(diarrhea.df)))])
length(findCorrelation(corr.matrix.diarrhea.df, cutoff=0.9))
#no features met this correlation cutoff, so it's ok to move ahead with training
#the model with all these features.

#split the data into train and test sets
set.seed(12)
training.diarrhea <- createDataPartition(y=diarrhea.df$infant.diarrhea.in.visit3.or.visit4,p = 0.7, list = FALSE)
diarrhea.df.train <- diarrhea.df[training.diarrhea, ]
diarrhea.df.test <- diarrhea.df[-training.diarrhea, ]
table(diarrhea.df.train$infant.diarrhea.in.visit3.or.visit4)
table(diarrhea.df.test$infant.diarrhea.in.visit3.or.visit4) 

diarrhea.formula <- paste(collapse="","infant.diarrhea.in.visit3.or.visit4 ~", paste(collapse=" + ",c(grep(pattern="^[a-z]__|infant.diarrhea$|household.size$", colnames(diarrhea.df), value=T))))
print(diarrhea.formula)

diarrhea.RF.classification.model <- ranger(formula = diarrhea.formula, data=diarrhea.df.train, write.forest=TRUE, oob.error=TRUE, classification = TRUE,seed=123, importance = "permutation")

#step 3: assess feature importance in the model
importance(diarrhea.RF.classification.model)
#this value is the increase in the model's prediction error if the values for that feature were randomly scrambled.
#so a positive and higher value means the feature was important to the model's prediction power. 
diarrhea.RF.feature.importance <- data.frame(importance(diarrhea.RF.classification.model))
diarrhea.RF.feature.importance$feature <- rownames(diarrhea.RF.feature.importance)
colnames(diarrhea.RF.feature.importance) <- c("feature.importance", "feature")

diarrhea.RF.feature.importance %>%
    mutate(feature = fct_reorder(feature, feature.importance)) %>%
    ggplot(aes(y=feature, x=feature.importance))+geom_col()


diarrhea.RF.feature.importance %>%
    mutate(feature = fct_reorder(feature, feature.importance)) %>%
    ggplot(aes(y=feature, x=feature.importance))+geom_col(aes(fill=feature.importance), color="black")+ scale_fill_gradient(low = "black", high="yellow")+scale_y_discrete(labels=c("Prevotella timonensis", "Flavobacteriales", "Rhodoferax", "Sutterella wadsworthensis", "infant diarrhea (in visit 2)", "household size (<= 3 or > 3)", "Haemophilus","Bacteroides uniformis","Corynebacteriales", "Staphylococcales", "Bacteroides vulgatus"))+xlab("permutation feature importance")+labs(fill="feature importance")
ggsave(filename="diarrhea.RF.model.feature.importance.jpeg", plot=last_plot(), width=8, height=6, units=c("in"), dpi=600)


#step 4: make predictions on test data and assess model performance
diarrhea.RF.classification.model.prediction <- predict(object = diarrhea.RF.classification.model, data=diarrhea.df.test, type = "response", seed=123)

confmatrix.diarrhea <-confusionMatrix(data=as.factor(diarrhea.RF.classification.model.prediction$predictions), reference=as.factor(diarrhea.df.test$infant.diarrhea.in.visit3.or.visit4), positive="1")

#determine accuracy, balanced accuracy, and F1-score 
confmatrix.diarrhea$byClass #balanced accuracy 0.646, F1-score 0.462
confmatrix.diarrhea$overall #accuracy 0.781

#determine ROC AUC
diarrhea.RF.probability.model <- ranger(formula = diarrhea.formula, data=diarrhea.df.train, write.forest=TRUE, oob.error=TRUE, probability = TRUE,seed=123, importance = "permutation")
diarrhea.RF.probabiliy.model.probabilities <- predict(object = diarrhea.RF.probability.model, data=diarrhea.df.test, type = "response", seed=123)
ROCpred.diarrhea.RF <- prediction(diarrhea.RF.probabiliy.model.probabilities$predictions[,2], diarrhea.df.test$infant.diarrhea.in.visit3.or.visit4)
ROCper.diarrhea.RF <- performance(ROCpred.diarrhea.RF, measure = "auc")
print(ROCper.diarrhea.RF@y.values) #auc of 0.615, which is a little better than chance.

#step 5: calculate (and visualize the Shapley values) to determine feature contribution to model predictions.
pfun <- function(object, newdata){
    predict(object, newdata)$predictions[,2]
}

shap.output.diarrhea.RF <-fastshap::explain(object=diarrhea.RF.probability.model, pred_wrapper = pfun, X=diarrhea.df.train[,c(grep(pattern="^[a-z]__|infant.diarrhea$|household.size$", colnames(diarrhea.df.train)))],nsim=100, shap_only=F, adjust=T)

shap.viz.diarrhea.RF <- shapviz(shap.output.diarrhea.RF, X=diarrhea.df.train[,c(grep(pattern="^[a-z]__|infant.diarrhea$|household.size$", colnames(diarrhea.df.train)))])
shap.viz.diarrhea.RF.svimportance <-sv_importance(shap.viz.diarrhea.RF, max_display = Inf, kind="beeswarm")
shap.viz.diarrhea.RF.svimportance
shap.viz.diarrhea.RF.svimportance+scale_y_discrete(labels=c("Rhodoferax","Flavobacteriales","Prevotella timonensis", "infant diarrhea (in visit 2)", "Sutterella wadsworthensis", "household size (<= 3 or > 3)","Corynebacteriales","Bacteroides vulgatus","Bacteroides uniformis","Staphylococcales", "Haemophilus"))
ggsave(filename="diarrhea.RF.model.SHAPvalue.jpeg", plot=last_plot(), width=8, height=6, units=c("in"), dpi=600)

#step 6: look at the abundances of Bacteroides vulgatus, Staphylococcales, and Haemophilus
#since they were either the top two important features to the model, or the top
#two important features in predicting diarrhea cases in the test data set (i.e. high on the SHAP value list).  
melted.alltaxa.Bvulgatus <- melted.alltaxa[c(grep(pattern="Bacteroides_vulgatus", melted.alltaxa$taxon)),]
melted.alltaxa.Staphylococcales <- melted.alltaxa[melted.alltaxa$taxon==" o__Staphylococcales",]
melted.alltaxa.Haemophilus <- melted.alltaxa[melted.alltaxa$taxon==" g__Haemophilus",]

ggplot(aes(x=infant.diarrhea.in.visit3.or.visit4, y=Abundance, color=as.factor(infant.diarrhea.in.visit3.or.visit4)),data=melted.alltaxa.Bvulgatus)+geom_boxplot()+scale_x_discrete(labels=c("no\nn=80", "yes\nn=29"))+scale_color_manual(values=c("black", "#F8766D"))+theme(legend.position="none")+xlab("diarrhea in visits 3 or 4")+ylab("B. vulgatus relative abundance (%)")+facet_zoom(ylim=c(0,5))
ggsave(filename="diarrheaVisit3or4.vs.Bvulgatus.jpeg", plot=last_plot(), dpi=600, width=4, height=4)
ggplot(aes(x=infant.diarrhea.in.visit3.or.visit4, y=Abundance, color=as.factor(infant.diarrhea.in.visit3.or.visit4)),data=melted.alltaxa.Staphylococcales)+geom_boxplot()+scale_x_discrete(labels=c("no\nn=80", "yes\nn=29"))+scale_color_manual(values=c("black", "#F8766D"))+theme(legend.position="none")+xlab("diarrhea in visits 3 or 4")+ylab("Staphylococcales relative abundance (%)")+facet_zoom(ylim=c(0,5))
ggsave(filename="diarrheaVisit3or4.vs.Staphylococcales.jpeg", plot=last_plot(), dpi=600, width=4, height=4)
ggplot(aes(x=infant.diarrhea.in.visit3.or.visit4, y=Abundance, color=as.factor(infant.diarrhea.in.visit3.or.visit4)),data=melted.alltaxa.Haemophilus)+geom_boxplot()+scale_x_discrete(labels=c("no\nn=80", "yes\nn=29"))+scale_color_manual(values=c("black", "#F8766D"))+theme(legend.position="none")+xlab("diarrhea in visits 3 or 4")+ylab("Haemophilus relative abundance (%)")+facet_zoom(ylim=c(0,2.5))
ggsave(filename="diarrheaVisit3or4.vs.Haemophilus.jpeg", plot=last_plot(), dpi=600, width=4, height=4)

wilcox.test(x=melted.alltaxa.Bvulgatus[melted.alltaxa.Bvulgatus$infant.diarrhea.in.visit3.or.visit4=="no","Abundance"], y=melted.alltaxa.Bvulgatus[melted.alltaxa.Bvulgatus$infant.diarrhea.in.visit3.or.visit4=="yes","Abundance"], alternative="less") #not significant
wilcox.test(x=melted.alltaxa.Staphylococcales[melted.alltaxa.Staphylococcales$infant.diarrhea.in.visit3.or.visit4=="no","Abundance"], y=melted.alltaxa.Staphylococcales[melted.alltaxa.Staphylococcales$infant.diarrhea.in.visit3.or.visit4=="yes","Abundance"], alternative="less") #significant, p-value = 0.001443
wilcox.test(x=melted.alltaxa.Haemophilus[melted.alltaxa.Haemophilus$infant.diarrhea.in.visit3.or.visit4=="no","Abundance"], y=melted.alltaxa.Haemophilus[melted.alltaxa.Haemophilus$infant.diarrhea.in.visit3.or.visit4=="yes","Abundance"], alternative="less") #significant, p-value = 0.01206
#p-value adjustment for these three tests
p.adjust(p=c(0.07642, 0.001443, 0.01206), method="fdr")
#resulting p-values are: 0.076420 (B. valgatus), 0.004329 (Staphylococcales), 0.018090 (Haemophilus)

#Staphylococcales and Haemophilus were significantly higher in visit 2 for infants that had diarrhea
#in visits 3 or 4 versus those that didn't have diarrhea in visits 3 or 4. 

#step 7: make a combined plot (A, B, C) of the feature importance and Staphylococcales and Haemophilus graphs
diarrhea.RF.feature.importance.plot <- diarrhea.RF.feature.importance %>%
    mutate(feature = fct_reorder(feature, feature.importance)) %>%
    ggplot(aes(y=feature, x=feature.importance))+geom_col(aes(fill=feature.importance), color="black")+ scale_fill_gradient(low = "black", high="yellow")+scale_y_discrete(labels=c("Prevotella timonensis", "Flavobacteriales", "Rhodoferax", "Sutterella wadsworthensis", "infant diarrhea in visit 2", "household size\n(<= 3 or > 3)", "Haemophilus","Bacteroides uniformis","Corynebacteriales", "Staphylococcales", "Bacteroides vulgatus"))+xlab("permutation feature importance")+labs(fill="feature importance")

shap.viz.diarrhea.RF.svimportance.plot<- shap.viz.diarrhea.RF.svimportance+scale_y_discrete(labels=c("Rhodoferax","Flavobacteriales","Prevotella timonensis", "infant diarrhea in visit 2", "Sutterella wadsworthensis", "household size\n(<= 3 or > 3)","Corynebacteriales","Bacteroides vulgatus","Bacteroides uniformis","Staphylococcales", "Haemophilus"))
diarrhea.staphylococcales <-ggplot(aes(x=infant.diarrhea.in.visit3.or.visit4, y=Abundance, color=as.factor(infant.diarrhea.in.visit3.or.visit4)),data=melted.alltaxa.Staphylococcales)+geom_boxplot()+scale_x_discrete(labels=c("no\nn=80", "yes\nn=29"))+scale_color_manual(values=c("black", "#F8766D"))+theme(legend.position="none")+xlab("diarrhea in visits 3 or 4")+ylab("Staphylococcales relative abundance (%)")+facet_zoom(ylim=c(0,5))
diarrhea.haemophilus <-ggplot(aes(x=infant.diarrhea.in.visit3.or.visit4, y=Abundance, color=as.factor(infant.diarrhea.in.visit3.or.visit4)),data=melted.alltaxa.Haemophilus)+geom_boxplot()+scale_x_discrete(labels=c("no\nn=80", "yes\nn=29"))+scale_color_manual(values=c("black", "#F8766D"))+theme(legend.position="none")+xlab("diarrhea in visits 3 or 4")+ylab("Haemophilus relative abundance (%)")+facet_zoom(ylim=c(0,2.5))

plot_grid(diarrhea.RF.feature.importance.plot,shap.viz.diarrhea.RF.svimportance.plot, diarrhea.staphylococcales, diarrhea.haemophilus, ncol=2, nrow=2, labels=LETTERS[1:4])
ggsave("diarrhea.RF.model.featureimport.SHAP.Staphylococcales.Haemophilus.jpeg", plot=last_plot(), width=14, height=10, units=c("in"), dpi=600)

#Part 3: Among the taxa flagged as potentially associated with fever in visit 3,4,
#look for associations using random forest.  

#step 1: select taxa flagged by taxaHFE 2.0 for fever.
colnames(taxaHFE.output.fever.visit3and4)
taxaHFE.output.fever.taxa <- colnames(taxaHFE.output.fever.visit3and4[,3:7])
print(taxaHFE.output.fever.taxa)
#s_human_gut is not a unique species name, but we have added the full taxonomy
#to the name so it's now unique. 
#Need to look at how this taxon is named in qiime feature table. 
melted.alltaxa[c(grep(pattern="human_gut", melted.alltaxa$taxon)),"taxon"] %>% unique()
taxaHFE.output.fever.taxa <- "negativicutes$|clostridium_sensu_stricto_1_ s__human_gut$|bacilli$|bacteroides_fragilis$|actinobacteriota$"
fever.df <- melted.alltaxa[c(grep(pattern=taxaHFE.output.fever.taxa, melted.alltaxa$taxon, ignore.case = T)),]
table(fever.df$taxon)

#step 2: run random forest with these taxa. 

#change data frame into wide format
fever.df <- pivot_wider(data=fever.df, names_from="taxon", values_from="Abundance")

#remove leading space at beginning of taxon name
colnames(fever.df) <- gsub(pattern=" ", replacement="", colnames(fever.df))

#check the class type of feature
str(fever.df) #all taxa were numeric and all categorical were character

#change the values for "infant.fever.in.visit3.or.visit4" from "yes" and "no"
#to 1 and 0, respectively. 
fever.df$infant.fever.in.visit3.or.visit4 <- gsub(pattern="no", replacement="0", ignore.case = T, fever.df$infant.fever.in.visit3.or.visit4)
fever.df$infant.fever.in.visit3.or.visit4 <- gsub(pattern="yes", replacement="1", ignore.case = T, fever.df$infant.fever.in.visit3.or.visit4)
fever.df$infant.fever.in.visit3.or.visit4 <- as.factor(fever.df$infant.fever.in.visit3.or.visit4)

#determine if any of the features (aka predictors) are highly correlated. This matters because if 
#any of the features are highly correlated then the feature importance values/shapley values
#will not be accurate for the highly correlated features. 
ggplot(aes(x=infant.fever, y=household.size), data=fever.df)+geom_jitter()
#change the features to numeric values to determine if they are
#correlated with each other
fever.df$infant.fever.numeric <- fever.df$infant.fever
fever.df[fever.df$infant.fever.numeric=="yes", "infant.fever.numeric"] <- "1"
fever.df[fever.df$infant.fever.numeric=="no", "infant.fever.numeric"] <- "0"
fever.df$infant.fever.numeric <- as.numeric(fever.df$infant.fever.numeric)

fever.df$household.size.numeric <- fever.df$household.size
fever.df[fever.df$household.size.numeric=="> 3", "household.size.numeric"] <- "1"
fever.df[fever.df$household.size.numeric=="<= 3", "household.size.numeric"] <- "0"
fever.df$household.size.numeric <- as.numeric(fever.df$household.size.numeric)

corr.matrix.fever.df <- cor(x=fever.df[,c(grep(pattern="household.size.numeric|infant.fever.numeric$|^[sgfocpd]_", colnames(fever.df)))])
length(findCorrelation(corr.matrix.fever.df, cutoff=0.9))
#no features met this correlation cutoff, so it's ok to move ahead with training
#the model with all these features.


#split the data into train and test sets
set.seed(123)
training.fever <- createDataPartition(y=fever.df$infant.fever.in.visit3.or.visit4,p = 0.7, list = FALSE)
fever.df.train <- fever.df[training.fever, ]
fever.df.test <- fever.df[-training.fever, ]
table(fever.df.train$infant.fever.in.visit3.or.visit4)
table(fever.df.test$infant.fever.in.visit3.or.visit4) 

fever.formula <- paste(collapse="","infant.fever.in.visit3.or.visit4 ~", paste(collapse=" + ",c(grep(pattern="^[a-z]__|infant.fever$|household.size$", colnames(fever.df), value=T))))
print(fever.formula)

fever.RF.classification.model <- ranger(formula = fever.formula, data=fever.df.train, write.forest=TRUE, oob.error=TRUE, classification = TRUE, seed=123, importance = "permutation")

#step 3: assess feature importance in the model
importance(fever.RF.classification.model)
#this value is the increase in the model's prediction error if the values for that feature were randomly scrambled.
#so a positive and higher value means the feature was important to the model's prediction power. 
fever.RF.feature.importance <- data.frame(importance(fever.RF.classification.model))
fever.RF.feature.importance$feature <- rownames(fever.RF.feature.importance)
colnames(fever.RF.feature.importance) <- c("feature.importance", "feature")

fever.RF.feature.importance %>%
    mutate(feature = fct_reorder(feature, feature.importance)) %>%
    ggplot(aes(y=feature, x=feature.importance))+geom_col()

fever.RF.feature.importance %>%
    mutate(feature = fct_reorder(feature, feature.importance)) %>%
    ggplot(aes(y=feature, x=feature.importance))+geom_col(aes(fill=feature.importance), color="black")+ scale_fill_gradient(low = "black", high="yellow")+scale_y_discrete(labels=c("Bacilli", "Clostridium sensu stricto 1, species: human gut", "infant fever (in visit 2)", "Negativicutes", "household size (<= 3 or > 3)","Bacteroides fragilis","Actinobacteriota"))+xlab("permutation feature importance")+labs(fill="feature importance")

ggsave(filename="fever.RF.model.feature.importance.jpeg", plot=last_plot(), width=9, height=4, units=c("in"), dpi=600)

#step 4: make predictions on test data and assess model performance
fever.RF.classification.model.prediction <- predict(object = fever.RF.classification.model, data=fever.df.test, type = "response", seed=123)

confmatrix.fever <-confusionMatrix(data=as.factor(fever.RF.classification.model.prediction$predictions), reference=as.factor(fever.df.test$infant.fever.in.visit3.or.visit4), positive="1")

#determine accuracy, balanced accuracy, and F1-score 
confmatrix.fever$byClass #balanced accuracy 0.603, F1-score 0.667
confmatrix.fever$overall #accuracy 0.613

#determine ROC AUC
fever.RF.probability.model <- ranger(formula = fever.formula, data=fever.df.train, write.forest=TRUE, oob.error=TRUE, probability = TRUE,seed=123, importance = "permutation")
fever.RF.probabiliy.model.probabilities <- predict(object = fever.RF.probability.model, data=fever.df.test, type = "response", seed=123)
ROCpred.fever.RF <- prediction(fever.RF.probabiliy.model.probabilities$predictions[,2], fever.df.test$infant.fever.in.visit3.or.visit4)
ROCper.fever.RF <- performance(ROCpred.fever.RF, measure = "auc")
print(ROCper.fever.RF@y.values) #auc of 0.620, which is a little better than chance.

#step 5: calculate (and visualize the Shapley values) to determine feature contribution to model predictions.
pfun <- function(object, newdata){
    predict(object, newdata)$predictions[,2]
}

shap.output.fever.RF <-fastshap::explain(object=fever.RF.probability.model, pred_wrapper = pfun, X=fever.df.train[,c(grep(pattern="^[a-z]__|infant.fever$|household.size$", colnames(fever.df.train)))],nsim=100, shap_only=F, adjust=T)

shap.viz.fever.RF <- shapviz(shap.output.fever.RF, X=fever.df.train[,c(grep(pattern="^[a-z]__|infant.fever$|household.size$", colnames(fever.df.train)))])
shap.viz.fever.RF.svimportance <-sv_importance(shap.viz.fever.RF, max_display = Inf, kind="beeswarm")
shap.viz.fever.RF.svimportance
shap.viz.fever.RF.svimportance+scale_y_discrete(labels=c("household size (<= 3 or > 3)","infant fever (in visit 2)", "Clostridium sensu stricto 1, species: human gut", "Bacilli","Bacteroides fragilis","Negativicutes", "Actinobacteriota"))
ggsave(filename="fever.RF.model.SHAPvalue.jpeg", plot=last_plot(), width=9, height=5, units=c("in"), dpi=600)

#step 6: look at the abundance of Actinobacteriota since it had by far the greatest
#feature importance and SHAP value. 
ggplot(aes(x=infant.fever.in.visit3.or.visit4, y=Abundance),data=melted.alltaxa[melted.alltaxa$taxon==" p__Actinobacteriota",])+geom_boxplot()+scale_x_discrete(labels=c("no\nn=46", "yes\nn=63"))+scale_color_manual(values=c("black", "#00BA38"))+theme(legend.position = "none")+xlab("infant fever in visits 3 or 4")+ylab("Actinobacteriota relative abundance (%)")
#using wilcoxon rank sum test to test for greater actinobacteriota with infant fever in visit 3 or 4. 
wilcox.test(x=melted.alltaxa[melted.alltaxa$taxon==" p__Actinobacteriota" & melted.alltaxa$infant.fever.in.visit3.or.visit4=="no","Abundance"], y=melted.alltaxa[melted.alltaxa$taxon==" p__Actinobacteriota" & melted.alltaxa$infant.fever.in.visit3.or.visit4=="yes","Abundance"], alternative="less")
#actinobacteriota was significantly higher in visit 2 for infants that had fever in visits 3 or 4
#versus those that didn't get fever in visits 3 or 4. 

#step 7: make a combined plot (A, B, C) of the feature importances and the Actinobacteriota
fever.RF.feature.importance.plot <- fever.RF.feature.importance %>%
    mutate(feature = fct_reorder(feature, feature.importance)) %>%
    ggplot(aes(y=feature, x=feature.importance))+geom_col(aes(fill=feature.importance), color="black")+ scale_fill_gradient(low = "black", high="yellow")+scale_y_discrete(labels=c("Bacilli", "Clostridium sensu stricto 1,\nspecies: human gut", "infant fever in visit 2", "Negativicutes", "household size\n(<= 3 or > 3)","Bacteroides fragilis","Actinobacteriota"))+xlab("permutation feature importance")+labs(fill="feature importance")

shap.viz.fever.RF.svimportance.plot <- shap.viz.fever.RF.svimportance+scale_y_discrete(labels=c("household size\n(<= 3 or > 3)","infant fever in visit 2", "Clostridium sensu stricto 1,\nspecies: human gut", "Bacilli","Bacteroides fragilis","Negativicutes", "Actinobacteriota"))
fever.actinobacteriota <- ggplot(aes(x=infant.fever.in.visit3.or.visit4, color=as.factor(infant.fever.in.visit3.or.visit4), y=Abundance),data=melted.alltaxa[melted.alltaxa$taxon==" p__Actinobacteriota",])+geom_boxplot()+scale_x_discrete(labels=c("no\nn=46", "yes\nn=63"))+scale_color_manual(values=c("black", "#00BA38"))+theme(legend.position = "none")+xlab("infant fever in visits 3 or 4")+ylab("Actinobacteriota relative abundance (%)")
plot_grid(fever.RF.feature.importance.plot, shap.viz.fever.RF.svimportance.plot, fever.actinobacteriota, ncol=3, nrow=1, labels=LETTERS[1:3], rel_widths = c(2,2,0.8))
ggsave("fever.RF.model.featureimport.SHAP.actinobacteriota.jpeg", plot=last_plot(), width=16, height=5, units=c("in"), dpi=600)


#Part 4: Among the taxa flagged as potentially associated with vomit in visit 3,4,
#look for associations using random forest.  

#step 1: select taxa flagged by taxaHFE 2.0 for vomit.
colnames(taxaHFE.output.vomit.visit3and4)
taxaHFE.output.vomit.taxa <- colnames(taxaHFE.output.vomit.visit3and4[,3:20])
print(taxaHFE.output.vomit.taxa)
#s_uncultured_bacteroidetes (from g_sediminibacterium) is a non-unique species name 
#so checking the taxonomies of other species named s_uncultured_bacteroidetes
melted.alltaxa[c(grep(pattern="uncultured_bacteroidetes", melted.alltaxa$taxon, ignore.case=T)), "taxon"] %>% unique()
#there is one other species with this name but they differ at the genus level so 
#will include genus name to distinguish the species. 
taxaHFE.output.vomit.taxa <- "veillonella_ratti$|clostridium_butyricum$|erysipelatoclostridium_ramosum$|streptococcus_agalactiae$|lactobacillus_reuteri$|gemellaceae$|bacteroides_caccae$|bacteroides_thetaiotaomicron$|parabacteroides_distasonis$|flavobacteriales$|sediminibacterium_ s__uncultured_bacteroidetes$|actinomyces$|corynebacterium_pseudodiphtheriticum$|bifidobacterium_dentium$|bifidobacterium_bifidum$|frankiales$|haemophilus$|raoultella$"
vomit.df <- melted.alltaxa[c(grep(pattern=taxaHFE.output.vomit.taxa, melted.alltaxa$taxon, ignore.case = T)),]
table(vomit.df$taxon) 
length(unique(vomit.df$taxon)) #18 taxa

#step 2: run random forest with these taxa. 

#change data frame into wide format
vomit.df <- pivot_wider(data=vomit.df, names_from="taxon", values_from="Abundance")

#remove leading space at beginning of taxon name, and replace "-" with "_"
colnames(vomit.df) <- gsub(pattern=" ", replacement="", colnames(vomit.df))
colnames(vomit.df) <- gsub(pattern="-", replacement="_", colnames(vomit.df))

#check the class type of feature
str(vomit.df) #all taxa were numeric and all categorical were character

#change the values for "infant.vomit.in.visit3.or.visit4" from "yes" and "no"
#to 1 and 0, respectively. 
vomit.df$infant.vomit.in.visit3.or.visit4 <- gsub(pattern="no", replacement="0", ignore.case = T, vomit.df$infant.vomit.in.visit3.or.visit4)
vomit.df$infant.vomit.in.visit3.or.visit4 <- gsub(pattern="yes", replacement="1", ignore.case = T, vomit.df$infant.vomit.in.visit3.or.visit4)
vomit.df$infant.vomit.in.visit3.or.visit4 <- as.factor(vomit.df$infant.vomit.in.visit3.or.visit4)

#determine if any of the features (aka predictors) are highly correlated. This matters because if 
#any of the features are highly correlated then the feature importance values/shapley values
#will not be accurate for the highly correlated features. 
ggplot(aes(x=infant.vomit, y=household.size), data=vomit.df)+geom_jitter()
#change the features to numeric values to determine if they are
#correlated with each other
vomit.df$infant.vomit.numeric <- vomit.df$infant.vomit
vomit.df[vomit.df$infant.vomit.numeric=="yes", "infant.vomit.numeric"] <- "1"
vomit.df[vomit.df$infant.vomit.numeric=="no", "infant.vomit.numeric"] <- "0"
vomit.df$infant.vomit.numeric <- as.numeric(vomit.df$infant.vomit.numeric)

vomit.df$household.size.numeric <- vomit.df$household.size
vomit.df[vomit.df$household.size.numeric=="> 3", "household.size.numeric"] <- "1"
vomit.df[vomit.df$household.size.numeric=="<= 3", "household.size.numeric"] <- "0"
vomit.df$household.size.numeric <- as.numeric(vomit.df$household.size.numeric)

corr.matrix.vomit.df <- cor(x=vomit.df[,c(grep(pattern="household.size.numeric|infant.vomit.numeric$|^[sgfocpd]_", colnames(vomit.df)))])
length(findCorrelation(corr.matrix.vomit.df, cutoff=0.9))
#two correlations met this correlation cutoff, so need to determine which feature(s)
#must be removed before moving forward. 
findCorrelation(corr.matrix.vomit.df, cutoff=0.9)
#columns 12 and 18
corr.vomit.df <- as.data.frame(corr.matrix.vomit.df)
vomit.df.feature.to.remove <-corr.vomit.df[,c(12,18)]
vomit.df.feature.to.remove <-vomit.df.feature.to.remove[(vomit.df.feature.to.remove[1]>=0.9)|(vomit.df.feature.to.remove[2]>=0.9),]
print(vomit.df.feature.to.remove)
#Frankiales Order is highly correlated with uncultured Bacteroidetes species and Flavobacteriales
#so removing this feature. 
dim(vomit.df)
vomit.df <- vomit.df[,-c(grep(pattern="o__Frankiales", colnames(vomit.df)))]
dim(vomit.df)
colnames(vomit.df)

#split the data into train and test sets
set.seed(123)
training.vomit <- createDataPartition(y=vomit.df$infant.vomit.in.visit3.or.visit4,p = 0.7, list = FALSE)
vomit.df.train <- vomit.df[training.vomit, ]
vomit.df.test <- vomit.df[-training.vomit, ]
table(vomit.df.train$infant.vomit.in.visit3.or.visit4)
table(vomit.df.test$infant.vomit.in.visit3.or.visit4) 

vomit.formula <- paste(collapse="","infant.vomit.in.visit3.or.visit4 ~", paste(collapse=" + ",c(grep(pattern="^[a-z]__|infant.vomit$|household.size$", colnames(vomit.df), value=T))))
print(vomit.formula)

vomit.RF.classification.model <- ranger(formula = vomit.formula, data=vomit.df.train, write.forest=TRUE, oob.error=TRUE, classification = TRUE, seed=123, importance = "permutation")

#step 3: assess feature importance in the model
importance(vomit.RF.classification.model)
#this value is the increase in the model's prediction error if the values for that feature were randomly scrambled.
#so a positive and higher value means the feature was important to the model's prediction power. 
vomit.RF.feature.importance <- data.frame(importance(vomit.RF.classification.model))
vomit.RF.feature.importance$feature <- rownames(vomit.RF.feature.importance)
colnames(vomit.RF.feature.importance) <- c("feature.importance", "feature")

vomit.RF.feature.importance %>%
    mutate(feature = fct_reorder(feature, feature.importance)) %>%
    ggplot(aes(y=feature, x=feature.importance))+geom_col()

vomit.RF.feature.importance %>%
    mutate(feature = fct_reorder(feature, feature.importance)) %>%
    ggplot(aes(y=feature, x=feature.importance))+geom_col(aes(fill=feature.importance), color="black")+ scale_fill_gradient(low = "black", high="yellow")+scale_y_discrete(labels=c("Haemophilus", "Erysipelatoclostridium ramosum", "Bacteroides thetaiotaomicron", "infant vomit (in visit 2)", "Veillonella ratti", "Flavobacteriales", "Gemellaceae", "Actinomyces", "Clostridium butyricum", "household size (<= 3 or > 3)", "Raoultella", "Lactobacillus reuteri", "Bifidobacterium dentium", "Sediminibacterium, species: uncultured Bacteroidetes", "Streptococcus agalactiae", "Bifidobacterium bifidum", "Parabacteroides distasonis", "Bacteroides caccae", "Corynebacterium pseudodiphtheriticum"))+xlab("permutation feature importance")+labs(fill="feature importance")

ggsave(filename="vomit.RF.model.feature.importance.jpeg", plot=last_plot(), width=9, height=5, units=c("in"), dpi=600)

#step 4: make predictions on test data and assess model performance
vomit.RF.classification.model.prediction <- predict(object = vomit.RF.classification.model, data=vomit.df.test, type = "response", seed=123)

confmatrix.vomit <-confusionMatrix(data=as.factor(vomit.RF.classification.model.prediction$predictions), reference=as.factor(vomit.df.test$infant.vomit.in.visit3.or.visit4), positive="1")

#determine accuracy, balanced accuracy, and F1-score 
confmatrix.vomit$byClass #balanced accuracy 0.480, F1-score NaN (which means no true positives were predicted)
confmatrix.vomit$overall #accuracy 0.750

#determine ROC AUC
vomit.RF.probability.model <- ranger(formula = vomit.formula, data=vomit.df.train, write.forest=TRUE, oob.error=TRUE, probability = TRUE,seed=123, importance = "permutation")
vomit.RF.probabiliy.model.probabilities <- predict(object = vomit.RF.probability.model, data=vomit.df.test, type = "response", seed=123)
ROCpred.vomit.RF <- prediction(vomit.RF.probabiliy.model.probabilities$predictions[,2], vomit.df.test$infant.vomit.in.visit3.or.visit4)
ROCper.vomit.RF <- performance(ROCpred.vomit.RF, measure = "auc")
print(ROCper.vomit.RF@y.values) #auc of 0.857, which is better than chance.

#step 5: calculate (and visualize the Shapley values) to determine feature contribution to model predictions.
pfun <- function(object, newdata){
    predict(object, newdata)$predictions[,2]
}

shap.output.vomit.RF <-fastshap::explain(object=vomit.RF.probability.model, pred_wrapper = pfun, X=vomit.df.train[,c(grep(pattern="^[a-z]__|infant.vomit$|household.size$", colnames(vomit.df.train)))],nsim=100, shap_only=F, adjust=T)

shap.viz.vomit.RF <- shapviz(shap.output.vomit.RF, X=vomit.df.train[,c(grep(pattern="^[a-z]__|infant.vomit$|household.size$", colnames(vomit.df.train)))])
shap.viz.vomit.RF.svimportance <-sv_importance(shap.viz.vomit.RF, max_display = Inf, kind="beeswarm")
shap.viz.vomit.RF.svimportance
shap.viz.vomit.RF.svimportance +scale_y_discrete(labels=c("Raoultella", "household size (<= 3 or > 3)", "Lactobacillus reuteri", "Veillonella ratti","Flavobacteriales", "Gemellaceae", "infant vomit (in visit 2)", "Sediminibacterium, species: uncultured Bacteroidetes", "Erysipelatoclostridium ramosum", "Streptococcus agalactiae", "Actinomyces", "Bifidobacterium dentium", "Clostridium butyricum", "Parabacteroides distasonis", "Bacteroides caccae", "Bacteroides thetaiotaomicron", "Corynebacterium pseudodiphtheriticum", "Haemophilus", "Bifidobacterium bifidum"))
ggsave(filename="vomit.RF.model.SHAPvalue.jpeg", plot=last_plot(), width=9, height=5, units=c("in"), dpi=600)

#make a combined plot (A, B) of the feature importances
vomit.RF.feature.importance.plot <- vomit.RF.feature.importance %>%
    mutate(feature = fct_reorder(feature, feature.importance)) %>%
    ggplot(aes(y=feature, x=feature.importance))+geom_col(aes(fill=feature.importance), color="black")+ scale_fill_gradient(low = "black", high="yellow")+scale_y_discrete(labels=c("Haemophilus", "Erysipelatoclostridium ramosum", "Bacteroides thetaiotaomicron", "infant vomit (in visit 2)", "Veillonella ratti", "Flavobacteriales", "Gemellaceae", "Actinomyces", "Clostridium butyricum", "household size (<= 3 or > 3)", "Raoultella", "Lactobacillus reuteri", "Bifidobacterium dentium", "Sediminibacterium, species: uncultured Bacteroidetes", "Streptococcus agalactiae", "Bifidobacterium bifidum", "Parabacteroides distasonis", "Bacteroides caccae", "Corynebacterium pseudodiphtheriticum"))+xlab("permutation feature importance")+labs(fill="feature importance")

shap.viz.vomit.RF.svimportance.plot<- shap.viz.vomit.RF.svimportance +scale_y_discrete(labels=c("Raoultella", "household size (<= 3 or > 3)", "Lactobacillus reuteri", "Veillonella ratti","Flavobacteriales", "Gemellaceae", "infant vomit (in visit 2)", "Sediminibacterium, species: uncultured Bacteroidetes", "Erysipelatoclostridium ramosum", "Streptococcus agalactiae", "Actinomyces", "Bifidobacterium dentium", "Clostridium butyricum", "Parabacteroides distasonis", "Bacteroides caccae", "Bacteroides thetaiotaomicron", "Corynebacterium pseudodiphtheriticum", "Haemophilus", "Bifidobacterium bifidum"))
plot_grid(vomit.RF.feature.importance.plot, shap.viz.vomit.RF.svimportance.plot, ncol=2, nrow=1, labels=LETTERS[1:2])
ggsave("vomit.RF.model.featureimport.SHAP.jpeg", plot=last_plot(), width=16, height=6, units=c("in"), dpi=600)

#step 6: look at the abundance of Corynebacterium pseudodiphtheriticum and B. bifidum
#since they had the greatest feature importance and SHAP value, respectively. 
ggplot(aes(x=infant.vomit.in.visit3.or.visit4, y=Abundance),data=melted.alltaxa[c(grep(pattern="s__Corynebacterium_pseudodiphtheriticum",melted.alltaxa$taxon)),])+geom_boxplot()
ggplot(aes(x=infant.vomit.in.visit3.or.visit4, y=Abundance),data=melted.alltaxa[c(grep(pattern="s__Bifidobacterium_bifidum",melted.alltaxa$taxon)),])+geom_boxplot()
#visually, there doesn't seem to be a difference in abundance for either taxon between the vomit outcome groups. 


