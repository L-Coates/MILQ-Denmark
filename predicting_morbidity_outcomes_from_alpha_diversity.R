#This script contains the code used to determine if alpha diversity is predictive 
#of morbidity outcomes. 

library(qiime2R)
library(tidyverse)
library(caret)
library(MASS)
library(ROCR)
library(ggpubr)

#PART 1: read in datasets containing metadata (with morbidity outcomes and covariates),
#and alpha diversity measurements.

##STEP 1: read in the metadata
metadata.part1 <- read.csv("../morbidity_medication_prevalences/infant.morbidities.medications.csv")
metadata.part1 <- metadata.part1[,-1]

metadata.part2 <- read.csv("../sociodemographics_and_sample_collection_info/birthmode.gender.parity.csv")

metadata.part3 <- read.csv("../complementary_feeding_vs_Granulicatella/Exclusive.breastfeeding.before.4months.csv")

metadata.v2 <- merge(x=metadata.part2, y=metadata.part1[,-c(2:20)], all=F, by="sample.id")
metadata.v2 <- merge(x=metadata.v2, y=metadata.part3, all=F, by="mid")
#put sample.id column back in first position
metadata.v2 <- metadata.v2[,c(2:55, 1)]


##STEP 2: read in the alpha-diversity measurements
shannon <- read_qza("../diversity/diversity_core_metrics_samples_327/shannon_vector.qza")
shannon <- as.data.frame(shannon$data)
shannon$sample.id <- rownames(shannon)

faith <- read_qza("../diversity/diversity_core_metrics_samples_327/faith_pd_vector.qza")
faith <- as.data.frame(faith$data)
colnames(faith) <- c("sample.id", "faith.diversity")

obsvOTU <- read_qza("../diversity/diversity_core_metrics_samples_327/observed_features_vector.qza")
obsvOTU <- obsvOTU$data
obsvOTU$sample.id <- rownames(obsvOTU)

evenness <- read_qza("../diversity/diversity_core_metrics_samples_327/evenness_vector.qza")
evenness <- evenness$data
evenness$sample.id <- rownames(evenness)


##STEP 3: make a table combining the metadata with the alpha diversity
alpha.diversity <- merge(x=evenness, y=faith, by="sample.id", all=F)
alpha.diversity <- merge(x=alpha.diversity, y=obsvOTU, by="sample.id", all=F)
alpha.diversity <- merge(x=alpha.diversity, y=shannon, by="sample.id", all=F)
dataset <- merge(x=metadata.v2, y=alpha.diversity, by="sample.id", all=F)

#PART 2: Run logistic regression to determine alpha-diversity
#predictors of diarrhea, fever, vomit in visit 3 and 4. 

##STEP 1: select visit 2 samples, and determine if evenness in stool microbiota 
#in visit 2 is a predictor of diarrhea, fever, or vomiting in visits 3 and 4.
dataset.visit2 <- dataset[dataset$visit=="2",]

#change the values for "infant.diarrhea.in.visit3.or.visit4" from "yes" and "no"
#to 1 and 0, respectively. 
dataset.visit2$infant.diarrhea.in.visit3.or.visit4 <- gsub(pattern="no", replacement="0", ignore.case = T, dataset.visit2$infant.diarrhea.in.visit3.or.visit4)
dataset.visit2$infant.diarrhea.in.visit3.or.visit4 <- gsub(pattern="yes", replacement="1", ignore.case = T, dataset.visit2$infant.diarrhea.in.visit3.or.visit4)
dataset.visit2$infant.diarrhea.in.visit3.or.visit4 <- as.factor(dataset.visit2$infant.diarrhea.in.visit3.or.visit4)
str(dataset.visit2$infant.diarrhea.in.visit3.or.visit4)
dataset.visit2$infant.diarrhea.in.visit3.or.visit4

dataset.visit2$infant.fever.in.visit3.or.visit4 <- gsub(pattern="no", replacement="0", ignore.case = T, dataset.visit2$infant.fever.in.visit3.or.visit4)
dataset.visit2$infant.fever.in.visit3.or.visit4 <- gsub(pattern="yes", replacement="1", ignore.case = T, dataset.visit2$infant.fever.in.visit3.or.visit4)
dataset.visit2$infant.fever.in.visit3.or.visit4 <- as.factor(dataset.visit2$infant.fever.in.visit3.or.visit4)
str(dataset.visit2$infant.fever.in.visit3.or.visit4)
dataset.visit2$infant.fever.in.visit3.or.visit4

dataset.visit2$infant.vomit.in.visit3.or.visit4 <- gsub(pattern="no", replacement="0", ignore.case = T, dataset.visit2$infant.vomit.in.visit3.or.visit4)
dataset.visit2$infant.vomit.in.visit3.or.visit4 <- gsub(pattern="yes", replacement="1", ignore.case = T, dataset.visit2$infant.vomit.in.visit3.or.visit4)
dataset.visit2$infant.vomit.in.visit3.or.visit4 <- as.factor(dataset.visit2$infant.vomit.in.visit3.or.visit4)
str(dataset.visit2$infant.vomit.in.visit3.or.visit4)
dataset.visit2$infant.vomit.in.visit3.or.visit4


#Find the total number of infants with each morbidity in visit 3 or 4, 
#because this will be used to determine how many variables can be included in the model
#to keep the number of "events per variable" to around 10. 
table(dataset.visit2$infant.diarrhea.in.visit3.or.visit4)
#29 infants had diarrhea in visit 3 or 4, and 80 did not. So about three variables could
#be included in the model. 
table(dataset.visit2$infant.fever.in.visit3.or.visit4)
#63 infants had fever in visit 3 or 4, and 46 did not. So about 4 variables could
#be included in the model.
table(dataset.visit2$infant.vomit.in.visit3.or.visit4)
#25 infants had vomit in visit 3 or 4, and 84 did not. So about 2 variables could
#be included in the model.

#Since 2-4 variables is close, we will choose three variables to include, especially
#because there are three variables that make the most biological sense to include:
#prior morbidity (i.e. visit 2), household size (i.e. infant's potential level of 
#exposure to infectious agents causing the morbidities of interest), and
#alpha diversity (our variable of interest).


#split the data into train and test sets
set.seed(123)
training.diarrhea <- createDataPartition(y=dataset.visit2$infant.diarrhea.in.visit3.or.visit4,p = 0.7, list = FALSE)
dataset.visit2.train.diarrhea  <- dataset.visit2[training.diarrhea, ]
dataset.visit2.test.diarrhea <- dataset.visit2[-training.diarrhea, ]
table(dataset.visit2.train.diarrhea$infant.diarrhea.in.visit3.or.visit4)
table(dataset.visit2.test.diarrhea$infant.diarrhea.in.visit3.or.visit4)

set.seed(123)
training.fever <- createDataPartition(y=dataset.visit2$infant.fever.in.visit3.or.visit4,p = 0.7, list = FALSE)
dataset.visit2.train.fever  <- dataset.visit2[training.fever, ]
dataset.visit2.test.fever <- dataset.visit2[-training.fever, ]
table(dataset.visit2.train.fever$infant.fever.in.visit3.or.visit4)
table(dataset.visit2.test.fever$infant.fever.in.visit3.or.visit4)

set.seed(123)
training.vomit <- createDataPartition(y=dataset.visit2$infant.vomit.in.visit3.or.visit4,p = 0.7, list = FALSE)
dataset.visit2.train.vomit  <- dataset.visit2[training.vomit, ]
dataset.visit2.test.vomit <- dataset.visit2[-training.vomit, ]
table(dataset.visit2.train.vomit$infant.vomit.in.visit3.or.visit4)
table(dataset.visit2.test.vomit$infant.vomit.in.visit3.or.visit4)


#fit the model (cannot include laxative in the model because there is only one 
#sample/infant that had it in visit 2 and the sample isn't in the training set)

set.seed(123)
dataset.visit2.model.diarrhea.evenness <- glm(infant.diarrhea.in.visit3.or.visit4 ~household.size+infant.diarrhea+pielou_evenness, data = dataset.visit2.train.diarrhea, family = binomial)

#summarize the final selected model
summary(dataset.visit2.model.diarrhea.evenness)
#None of these predictors had significant associations with diarrhea outcome  

#model for Pielou's evenness in visit 2 and later fever outcomes
set.seed(123)
dataset.visit2.model.fever.evenness <- glm(infant.fever.in.visit3.or.visit4 ~ household.size+infant.fever+pielou_evenness, data = dataset.visit2.train.fever, family = binomial)

#summarize the final selected model
summary(dataset.visit2.model.fever.evenness)
#Pielou's evenness had a significant negative association
#with fever outcome. 

exp(coef(dataset.visit2.model.fever.evenness))

#make predictions
dataset.visit2.probabilities.fever.evenness <- predict(dataset.visit2.model.fever.evenness, dataset.visit2.test.fever, type = "response")
dataset.visit2.predicted.classes.fever.evenness<- ifelse(dataset.visit2.probabilities.fever.evenness > 0.5, 1, 0)

confmatrix.fever.evenness <-confusionMatrix(as.factor(dataset.visit2.predicted.classes.fever.evenness), as.factor(dataset.visit2.test.fever$infant.fever.in.visit3.or.visit4), positive="1")

#determine accuracy, balanced accuracy, and F1-score 
confmatrix.fever.evenness$byClass #balanced accuracy 0.459, F1-score 0.579
confmatrix.fever.evenness$overall #accuracy 0.484

#determine ROC AUC
ROCPred.model.fever.evenness <- prediction(dataset.visit2.probabilities.fever.evenness, dataset.visit2.test.fever$infant.fever.in.visit3.or.visit4)
ROCPer.model.fever.evenness <- performance(ROCPred.model.fever.evenness, measure = "auc")
auc.model.fever.evenness <- ROCPer.model.fever.evenness@y.values[[1]]
print(auc.model.fever.evenness) #auc of 0.526, which is a little better than chance.

#plot evenness by fever
table(dataset.visit2$infant.fever.in.visit3.or.visit4) #46 no, 63 yes
ggplot(aes(x=as.factor(infant.fever.in.visit3.or.visit4), y=pielou_evenness, color=as.factor(infant.fever.in.visit3.or.visit4)), data=dataset.visit2)+geom_boxplot()+xlab("fever in visits 3 or 4")+ylab("Pielou's evenness")+scale_x_discrete(labels=c("no\nn=46", "yes\nn=63"))+scale_color_manual(values=c("black", "#00BA38"))+theme(legend.position = "none")
ggsave("evenness.visit2.vs.fever.visit3.4.jpeg",dpi=600, width=5, height=5)
hist(dataset.visit2$pielou_evenness)
shapiro.test(dataset.visit2$pielou_evenness)
#not normal distribution so finding best normalization before applying parametric test. 
bestNormalize::bestNormalize(dataset.visit2$pielou_evenness) #the best normalization technique was determined to be log_b(x+a) with b = 10
dataset.visit2$pielou_evenness_bestNormalized <- bestNormalize::bestNormalize(dataset.visit2$pielou_evenness)$x.t
shapiro.test(dataset.visit2$pielou_evenness_bestNormalized) #still not normal, so will have to use non-parametric test
wilcox.test(x=dataset.visit2[dataset.visit2$infant.fever.in.visit3.or.visit4=="0","pielou_evenness"], y= dataset.visit2[dataset.visit2$infant.fever.in.visit3.or.visit4=="1","pielou_evenness"], alternative="greater")
#significant difference in Pielou's evenness during visit 2 between infants that
#developed fever in visit 3 or 4 versus those that didn't develop fever. 


#evenness and vomit
set.seed(123)
dataset.visit2.model.vomit.evenness <- glm(infant.vomit.in.visit3.or.visit4 ~household.size+infant.vomit+pielou_evenness, data = dataset.visit2.train.vomit, family = binomial)
#summarize the final selected model
summary(dataset.visit2.model.vomit.evenness)
#no significant predictors. 

##STEP 2: among visit 2 samples, and determine if Faith's phylogenetic diversity
#in stool microbiota in visit 2 is a predictor of diarrhea, fever, or vomit in visits 3 and 4.

#fit the model (cannot include laxative in the model because there is only one 
#sample/infant that had it in visit 2 and the sample isn't in the training set)

#faith's PD in visit 2 as a predictor of diarrhea in visit 3, 4. 
set.seed(123)
dataset.visit2.model.diarrhea.faith <- glm(infant.diarrhea.in.visit3.or.visit4 ~household.size+infant.diarrhea+faith.diversity, data = dataset.visit2.train.diarrhea, family = binomial)
#summarize the final selected model
summary(dataset.visit2.model.diarrhea.faith)
#faith's PD was positively associated with diarrhea
exp(coef(dataset.visit2.model.diarrhea.faith))

dataset.visit2.probabilities.diarrhea.faith <- predict(dataset.visit2.model.diarrhea.faith, dataset.visit2.test.diarrhea, type = "response")
dataset.visit2.predicted.classes.diarrhea.faith<- ifelse(dataset.visit2.probabilities.diarrhea.faith > 0.5, 1, 0)

confmatrix.diarrhea.faith <-confusionMatrix(as.factor(dataset.visit2.predicted.classes.diarrhea.faith), as.factor(dataset.visit2.test.diarrhea$infant.diarrhea.in.visit3.or.visit4), positive="1")

#determine accuracy, balanced accuracy, and F1-score 
confmatrix.diarrhea.faith$byClass #balanced accuracy 0.500, F1-score 0.167. 
confmatrix.diarrhea.faith$overall #accuracy 0.688


#determine ROC AUC
ROCPred.model.diarrhea.faith <- prediction(dataset.visit2.probabilities.diarrhea.faith, dataset.visit2.test.diarrhea$infant.diarrhea.in.visit3.or.visit4)
ROCPer.model.diarrhea.faith <- performance(ROCPred.model.diarrhea.faith, measure = "auc")
auc.model.diarrhea.faith <- ROCPer.model.diarrhea.faith@y.values[[1]]
print(auc.model.diarrhea.faith) #auc of 0.526, which is a little better than chance.

#plot faith PD by diarrhea
ggplot(aes(x=as.factor(infant.diarrhea.in.visit3.or.visit4), y=faith.diversity, color=as.factor(infant.diarrhea.in.visit3.or.visit4)), data=dataset.visit2)+geom_boxplot()+xlab("diarrhea in visits 3 or 4")+ylab("faith's PD")+scale_x_discrete(labels=c("no\nn=80", "yes\nn=29"))+scale_color_manual(values=c("black", "#F8766D"))+theme(legend.position="none")
ggsave(plot=last_plot(), filename="faithPDvisit2.vs.diarrheaVisits34.jpeg")
wilcox.test(dataset.visit2[dataset.visit2$infant.diarrhea.in.visit3.or.visit4=="0","faith.diversity"], y= dataset.visit2[dataset.visit2$infant.diarrhea.in.visit3.or.visit4=="1","faith.diversity"],alternative="less")
#significant difference
hist(dataset.visit2$faith.diversity)
shapiro.test(dataset.visit2$faith.diversity)
#not normal distribution so finding best normalization before applying parametric test. 
bestNormalize::bestNormalize(dataset.visit2$faith.diversity) #the best normalization technique was determined to be box cox transformation
dataset.visit2$faith_diversity_bestNormalized <- bestNormalize::bestNormalize(dataset.visit2$faith.diversity)$x.t
shapiro.test(dataset.visit2$faith_diversity_bestNormalized) #now normal, so will use parametric test
t.test(x=dataset.visit2[dataset.visit2$infant.diarrhea.in.visit3.or.visit4=="0","faith_diversity_bestNormalized"], y= dataset.visit2[dataset.visit2$infant.diarrhea.in.visit3.or.visit4=="1","faith_diversity_bestNormalized"], alternative="less")
#significant difference


#faith's PD in visit 2 as predictor of fever in visit 3, 4.  
dataset.visit2.model.fever.faith <- glm(infant.fever.in.visit3.or.visit4 ~household.size+infant.fever+faith.diversity, data = dataset.visit2.train.fever, family = binomial)
#summarize the final selected model
summary(dataset.visit2.model.fever.faith)
#faith's PD was not significantly associated with fever outcome. 


#faith's pd in visit 2 to predict vomiting in visit 3, 4. 
set.seed(123)
dataset.visit2.model.vomit.faith <- glm(infant.vomit.in.visit3.or.visit4 ~household.size+infant.vomit+faith.diversity, data = dataset.visit2.train.vomit, family = binomial)
#summarize the final selected model
summary(dataset.visit2.model.vomit.faith)
#faith PD was not associated with vomit.  

##STEP 3: among visit 2 samples, determine if total observed ASVs
#in stool microbiota in visit 2 is a predictor of diarrhea, fever, or vomit in visits 3 and 4.
#fit the model (cannot include laxative in the model because there is only one 
#sample/infant that had it in visit 2 and the sample isn't in the training set)
set.seed(123)
dataset.visit2.model.diarrhea.obsvOTU <- glm(infant.diarrhea.in.visit3.or.visit4 ~household.size+infant.diarrhea+observed_features, data = dataset.visit2.train.diarrhea, family = binomial)
#summarize the final selected model
summary(dataset.visit2.model.diarrhea.obsvOTU)
#observed OTUs is not associated with diarrhea.

#total OTUs in visit 2 predicting fever in visit 3, 4. 
set.seed(123)
dataset.visit2.model.fever.obsvOTU <- glm(infant.fever.in.visit3.or.visit4 ~household.size+infant.fever+observed_features, data = dataset.visit2.train.fever, family = binomial)
#summarize the final selected model
summary(dataset.visit2.model.fever.obsvOTU)
#total OTUs didn't associate with fever. 

#total OTUs in visit 2 predicting vomit in visit 3, 4. 
set.seed(123)
dataset.visit2.model.vomit.obsvOTU <- glm(infant.vomit.in.visit3.or.visit4 ~household.size+infant.vomit+observed_features, data = dataset.visit2.train.vomit, family = binomial)
#summarize the final selected model
summary(dataset.visit2.model.vomit.obsvOTU)
#only total OTUs were not associated with vomit.  

##STEP 4: among visit 2 samples, determine if shannon entropy
#in stool microbiota in visit 2 is a predictor of diarrhea, fever or vomit in visits 3 and 4.
#fit the model (cannot include laxative in the model because there is only one 
#sample/infant that had it in visit 2 and the sample isn't in the training set)

#shannon entropy in visit 2 predicting diarrhea in visit 3, 4. 
set.seed(123)
dataset.visit2.model.diarrhea.shannon <- glm(infant.diarrhea.in.visit3.or.visit4 ~household.size+infant.diarrhea+shannon_entropy, data = dataset.visit2.train.diarrhea, family = binomial)
#summarize the final selected model
summary(dataset.visit2.model.diarrhea.shannon)
#shannon entropy was not associated with diarrhea

#shannon entropy in visit 2 predicting fever in visit 3, 4. 
set.seed(123)
dataset.visit2.model.fever.shannon <- glm(infant.fever.in.visit3.or.visit4 ~household.size+infant.fever+shannon_entropy, data = dataset.visit2.train.fever, family = binomial)
#summarize the final selected model
summary(dataset.visit2.model.fever.shannon)
#shannon entropy was negatively associated with fever. 

exp(coef(dataset.visit2.model.fever.shannon))

dataset.visit2.probabilities.fever.shannon <- predict(dataset.visit2.model.fever.shannon, dataset.visit2.test.fever, type = "response")
dataset.visit2.predicted.classes.fever.shannon<- ifelse(dataset.visit2.probabilities.fever.shannon > 0.5, 1, 0)

confmatrix.fever.shannon <-confusionMatrix(as.factor(dataset.visit2.predicted.classes.fever.shannon), as.factor(dataset.visit2.test.fever$infant.fever.in.visit3.or.visit4), positive="1")

#determine accuracy, balanced accuracy, and F1-score 
confmatrix.fever.shannon$byClass #balanced accuracy 0.459, F1-score 0.579. 
confmatrix.fever.shannon$overall #accuracy 0.484


#determine ROC AUC
ROCPred.model.fever.shannon <- prediction(dataset.visit2.probabilities.fever.shannon, dataset.visit2.test.fever$infant.fever.in.visit3.or.visit4)
ROCPer.model.fever.shannon <- performance(ROCPred.model.fever.shannon, measure = "auc")
auc.model.fever.shannon <- ROCPer.model.fever.shannon@y.values[[1]]
print(auc.model.fever.shannon) #auc of 0.496, which is a little better than chance.

#plot shannon entropy by fever in visit 3 or 4. 
ggplot(aes(x=as.factor(infant.fever.in.visit3.or.visit4), y=shannon_entropy, color=as.factor(infant.fever.in.visit3.or.visit4)), data=dataset.visit2)+geom_boxplot()+xlab("fever in visits 3 or 4")+ylab("Shannon entropy")+scale_x_discrete(labels=c("no\nn=46", "yes\nn=63"))+scale_color_manual(values=c("black", "#00BA38"))+theme(legend.position = "none")
ggsave(plot=last_plot(), filename="ShannonEntropyvisit2.vs.FeverVisits34.jpeg")
wilcox.test(dataset.visit2[dataset.visit2$infant.diarrhea.in.visit3.or.visit4=="0","shannon_entropy"], y= dataset.visit2[dataset.visit2$infant.diarrhea.in.visit3.or.visit4=="1","shannon_entropy"])
#significant difference
hist(dataset.visit2$shannon_entropy)
shapiro.test(dataset.visit2$shannon_entropy)
#normal distribution
t.test(x=dataset.visit2[dataset.visit2$infant.fever.in.visit3.or.visit4=="0","shannon_entropy"], y= dataset.visit2[dataset.visit2$infant.fever.in.visit3.or.visit4=="1","shannon_entropy"], "greater")
#significant difference

#shannon entropy in visit 2 predicting vomit in visit 3, 4. 
set.seed(123)
dataset.visit2.model.vomit.shannon <- glm(infant.vomit.in.visit3.or.visit4 ~household.size+infant.vomit+shannon_entropy, data = dataset.visit2.train.vomit, family = binomial)
#summarize the final selected model
summary(dataset.visit2.model.vomit.shannon)
#shannon entropy was not associated with vomit. 


#make a joint figure for alpha diversity metrics that associated with morbidity
ggplot.faithPD.diarrhea <- ggplot(aes(x=as.factor(infant.diarrhea.in.visit3.or.visit4), y=faith.diversity, color=as.factor(infant.diarrhea.in.visit3.or.visit4)), data=dataset.visit2)+
    geom_boxplot()+
    theme_bw()+
    xlab("diarrhea in\nvisits 3 or 4")+
    ylab("Faith's phylogenetic diversity")+
    scale_x_discrete(labels=c("no\nn=80", "yes\nn=29"))+
    scale_color_manual(values=c("black", "#F8766D"))+
    theme(legend.position="none", axis.title=element_text(size=9))
ggplot.evenness.fever <- ggplot(aes(x=as.factor(infant.fever.in.visit3.or.visit4), y=pielou_evenness, color=as.factor(infant.fever.in.visit3.or.visit4)), data=dataset.visit2)+
    geom_boxplot()+
    theme_bw()+
    xlab("fever in\nvisits 3 or 4")+
    ylab("Pielou's evenness")+
    scale_x_discrete(labels=c("no\nn=46", "yes\nn=63"))+
    scale_color_manual(values=c("black", "#00BA38"))+
    theme(legend.position = "none", axis.title=element_text(size=9))
ggplot.shannon.fever <- ggplot(aes(x=as.factor(infant.fever.in.visit3.or.visit4), y=shannon_entropy, color=as.factor(infant.fever.in.visit3.or.visit4)), data=dataset.visit2)+
    geom_boxplot()+
    theme_bw()+
    xlab("fever in\nvisits 3 or 4")+
    ylab("Shannon entropy")+
    scale_x_discrete(labels=c("no\nn=46", "yes\nn=63"))+
    scale_color_manual(values=c("black", "#00BA38"))+
    theme(legend.position = "none", axis.title=element_text(size=9))
ggarrange(ggplot.faithPD.diarrhea, ggplot.evenness.fever, ggplot.shannon.fever, ncol=3, nrow=1, labels="AUTO")
ggsave(plot=last_plot(), filename="morbidityVisits34.vs.alphaDiversityVisit2.tiff", device="tiff", height = 9, width = 14, units="cm")
