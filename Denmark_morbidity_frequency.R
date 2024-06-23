#This script describes the steps taken to determine the frequency of morbidity
#and medications in the Danish cohort.  

library(tidyverse)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(ggvenn)

#STEP 1: Read in data from REDCap and the metadata file for the 327 infant stool
#samples that are being used for analyses. 
redcap <- read.csv('../../../../REDCap_downloads/Denmark/MILQMAINSTUDYDENMARK_DATA_2022-10-11_2251.csv')
#we already know that there are two maternal IDs (MIDs) in the REDCap dataset that have 
#a lowercase check letter instead of the standard uppercase check letter. Correcting
#these MIDs now.
redcap[redcap$mid=="MD346m","mid"] <- "MD346M"
redcap[redcap$mid=="MD378x","mid"] <- "MD378X"

#Read in Kable lab data on Denmark samples processed and sent for sequencing
metadata <- read.csv("../sociodemographics_and_sample_collection_info/metadata-327-samples-with-select-stool-information.csv")

#cut down redcap data to the maternal IDs in the metadata file
redcap <- redcap[c(which(redcap$mid %in% metadata$mid)),]
    
#STEP 2: determine if all mothers completed morbidity questionnaires for themselves
#and their infants at each visit. 
unique(redcap$redcap_event_name)
table(redcap[redcap$redcap_event_name=="10__349_months_arm_1","f204_maternal_infant_morbidity_1_to_349m_complete"])
#11 were incomplete, 98 were complete
table(redcap[redcap$redcap_event_name=="35__599_months_arm_1","f304_maternal_infant_morbidity_35_to_599m_complete"])
#1 was incomplete, 108 were complete
table(redcap[redcap$redcap_event_name=="60__849_months_arm_1","f404_maternal_infant_morbidity_35_to_599m_complete"])
#5 were incomplete, and 104 were complete. 

#STEP 3:  Change data into wide format so that each child only has one row of data.  
#selecting the columns of interest
redcap.visit1 <- redcap[redcap$redcap_event_name=="2__3_days_postpart_arm_1",]
redcap.visit1 <- redcap.visit1[,c(1,grep(pattern="f101|f003", colnames(redcap.visit1)))]
redcap.visit2 <- redcap[redcap$redcap_event_name=="10__349_months_arm_1",]
redcap.visit2 <- redcap.visit2[,c(1,grep(pattern="f2", colnames(redcap.visit2)))]
redcap.visit3 <- redcap[redcap$redcap_event_name=="35__599_months_arm_1",]
redcap.visit3 <- redcap.visit3[,c(1,grep(pattern="f3", colnames(redcap.visit3)))]
redcap.visit4 <- redcap[redcap$redcap_event_name=="60__849_months_arm_1",]
redcap.visit4 <- redcap.visit4[,c(1,grep(pattern="f4", colnames(redcap.visit4)))]

redcap.wide <- merge(redcap.visit1, redcap.visit2, by="mid", all=F)
redcap.wide <- merge(redcap.wide, redcap.visit3, by="mid", all=F)
redcap.wide <- merge(redcap.wide, redcap.visit4, by="mid", all=F)

#STEP 4: look at the number of infants with an illness in visit 2. 

#diarrhea in the past week
table(redcap.wide$f204_diarrhoea_q5_1) #two infants had diarrhea

#diarrhea since birth, not including the past week. 
table(redcap.wide$f204_diarrhoea_q6_1) #seven infants had diarrhea


#vomiting in the past week
table(redcap.wide$f204_vomitn_q5_2) #two infants had vomiting

#vomiting since birth, not including the past week.
table(redcap.wide$f204_vomitn_q6_2) #four infants had vomiting

#fever in the past week
table(redcap.wide$f204_fever_q5_5) #five infants had fever

#fever since birth, not including the past week. 
table(redcap.wide$f204_fever_q6_5) #three infants had fever.

#other illnesses in the past week
table(redcap.wide$f204_other_q5_8) #four infants with an illness other than
#what was already listed. 
#looking specifically at these other illnesses
unique(redcap.wide$f204_specify_q5_8a)
#yeast infections, eye inflammation, stomachache which we won't include in 
#morbidity analyses. 

#other illnesses since birth, but not including the past week. 
table(redcap.wide$f204_other_q6_8) #fourteen infants were sick
unique(redcap.wide$f204_other_q6_8a)
#irritated lacrimal, jaundice, yeast infection, liquid in the ear canal, hip dysplasia,
#jaundice, abdominal cramps from week 2-7, newborn jaundice, eye infection, fungus in the skin,
#urticaria, constipation, and eye infection were the illnesses but won't be
#included in the analyses. 

table(redcap.wide$f204_infsck_q8) #three infants sick
unique(redcap.wide$f204_spillness1_q8_1_1)
#umbilical hernia, infection, jaundice were the illnesses listed
table(redcap.wide$f204_illness2_q8_1_2)#no additional illnesses.

#find the number of infants that had diarrhea in visit 2
visit2.diarrhea <- data.frame("mid"=redcap.wide$mid, "visit2.diarrhea"=rowSums(redcap.wide[,c(grep(pattern="f204_diarrhoea_q5_1$|f204_diarrhoea_q6_1$", colnames(redcap.wide)))], na.rm=T))
table(visit2.diarrhea$visit2.diarrhea>0)#8 infants total had diarrhea during visit 2. 

#find the number of infants that had vomiting in visit 2
visit2.vomit <- data.frame("mid"=redcap.wide$mid, "visit2.vomit"=rowSums(redcap.wide[,c(grep(pattern="f204_vomitn_q5_2$|f204_vomitn_q6_2$", colnames(redcap.wide)))], na.rm=T))
table(visit2.vomit$visit2.vomit>0)#six infants total vomited during visit 2. 

#find the number of infants that had fever in visit 2
visit2.fever <- data.frame("mid"=redcap.wide$mid, "visit2.fever"=rowSums(redcap.wide[,c(grep(pattern="f204_fever_q5_5$|f204_fever_q6_5$", colnames(redcap.wide)))], na.rm=T))
table(visit2.fever$visit2.fever>0)#eight infants total had fever during visit 2. 

#merge the diarrhea, vomit, and fever morbidities
visit2.morbidities <- merge(visit2.diarrhea, visit2.fever, by="mid", all=F)
visit2.morbidities <- merge(visit2.morbidities, visit2.vomit, by="mid", all=F)

#STEP 5: look at the number of infants with an illness in visit 3

#diarrhea in the past week
table(redcap.wide$f304_diarrhoea_q5_1) #three infants had diarrhea

#diarrhea since birth, not including the past week. 
table(redcap.wide$f304_diarrhoea_q6_1) #12 infants had diarrhea

#vomiting in the past week
table(redcap.wide$f304_vomitn_q5_2) #four infants had vomiting

#vomiting since birth, not including the past week.
table(redcap.wide$f304_vomitn_q6_2) #five infants had vomiting

#fever in the past week
table(redcap.wide$f304_fever_q5_5) #nine infants had fever

#fever since birth, not including the past week. 
table(redcap.wide$f304_fever_q6_5) #thirty-five infants had fever. 

#other illnesses in the past week
table(redcap.wide$f304_other_q5_8) #four infants with an illness other than
#what was already listed. 
#looking specifically at these other illnesses
unique(redcap.wide$f304_specify_q5_8a)
#thrush, stomachache when having feces, and eye infection
#which we won't include in morbidity analyses. 

#other illnesses since birth, but not including the past week. 
table(redcap.wide$f304_other_q6_8) #ten infants were sick
unique(redcap.wide$f304_other_q6_8a)
#yeast infection, rash, thrush in mouth, eye infection, stomachache when having feces,
#umbilical granuloma, eczema/fungus, inflammation of the toe, thrush were the illnesses 
#but won't be included in the analyses. 

table(redcap.wide$f304_infsck_q8) #one infant sick
unique(redcap.wide$f304_spillness1_q8_1_1)
#parechovirus was the illness listed and it can cause diarrhea and vomiting so 
#looking at the infant that had parechovirus and will see later if this same
#infant also had diarrhea or vomiting or fever during this visit 3. 
redcap.wide[c(grep(pattern="parecco", redcap.wide$f304_spillness1_q8_1_1)),"mid"] #MD070Y
table(redcap.wide$f304_illness2_q8_1_2)#no additional illnesses.

#find the number of infants that had diarrhea in visit 3
visit3.diarrhea <- data.frame("mid"=redcap.wide$mid, "visit3.diarrhea"=rowSums(redcap.wide[,c(grep(pattern="f304_diarrhoea_q5_1$|f304_diarrhoea_q6_1$", colnames(redcap.wide)))], na.rm=T))
table(visit3.diarrhea$visit3.diarrhea>0)#fifteen infants total had diarrhea during visit 3. 

#find the number of infants that had vomiting in visit 3
visit3.vomit <- data.frame("mid"=redcap.wide$mid, "visit3.vomit"=rowSums(redcap.wide[,c(grep(pattern="f304_vomitn_q5_2$|f304_vomitn_q6_2$", colnames(redcap.wide)))], na.rm=T))
table(visit3.vomit$visit3.vomit>0)#seven infants total vomited during visit 3. 

#find the number of infants that had fever in visit 3
visit3.fever <- data.frame("mid"=redcap.wide$mid, "visit3.fever"=rowSums(redcap.wide[,c(grep(pattern="f304_fever_q5_5$|f304_fever_q6_5$", colnames(redcap.wide)))], na.rm=T))
table(visit3.fever$visit3.fever>0)#fourty-one infants total had fever during visit 3. 

#merge the diarrhea, vomit, and fever morbidities
visit3.morbidities <- merge(visit3.diarrhea, visit3.fever, by="mid", all=F)
visit3.morbidities <- merge(visit3.morbidities, visit3.vomit, by="mid", all=F)

#now look to see if infant MD070Y which had parecheco virus also had diarrhea, fever or vomiting
visit3.morbidities[visit3.morbidities$mid=="MD070Y",] #the infant had fever. 

#STEP 6: look at the number of infants with an illness in visit 4

#diarrhea in the past week
table(redcap.wide$f404_diarrhoea_q5_1) #four infants had diarrhea

#diarrhea since birth, not including the past week. 
table(redcap.wide$f404_diarrhoea_q6_1) #13 infants had diarrhea

#vomiting in the past week
table(redcap.wide$f404_vomitn_q5_2) #three infants had vomiting

#vomiting since birth, not including the past week.
table(redcap.wide$f404_vomitn_q6_2) #eighteen infants had vomiting

#fever in the past week
table(redcap.wide$f404_fever_q5_5) #fourteen infants had fever

#fever since birth, not including the past week. 
table(redcap.wide$f404_fever_q6_5) #thirty-one infants had fever. 

#other illnesses in the past week
table(redcap.wide$f404_other_q5_8) #six infants with an illness other than
#what was already listed. 
#looking specifically at these other illnesses
unique(redcap.wide$f404_specify_q5_8a)
#constipation, rash all over the body-fever rash, thrush, blood in feces, virus rash.
#The blood in feces is notable however this could be from a variety of things.
#looking to see the infant/maternal ID in order to determine if this infant also
#had diarrhea, vomiting and/or fever. 
redcap.wide[redcap.wide$f404_specify_q5_8a=="blood in feces", "mid"] #MD312P

#other illnesses since birth, but not including the past week. 
table(redcap.wide$f404_other_q6_8) #two infants were sick
unique(redcap.wide$f404_other_q6_8a)
#constipated, constipation. But these won't be included in the analyses. 

table(redcap.wide$f404_infsck_q8) #two infants sick
unique(redcap.wide$f404_spillness1_q8_1_1)
# "See description", "during investigation -- blood spatters on the legs" 
#It isn't clear what description is being referred to, so taking note of this infant
#and looking to see if the infant was already noted as having some other specific 
#illness. 
redcap.wide[redcap.wide$f404_spillness1_q8_1_1=="See description","mid"] #MD280R

table(redcap.wide$f404_illness2_q8_1_2)#no additional illnesses.

#find the number of infants that had diarrhea in visit 4
visit4.diarrhea <- data.frame("mid"=redcap.wide$mid, "visit4.diarrhea"=rowSums(redcap.wide[,c(grep(pattern="f404_diarrhoea_q5_1$|f404_diarrhoea_q6_1$", colnames(redcap.wide)))], na.rm=T))
table(visit4.diarrhea$visit4.diarrhea>0)#seventeen infants total had diarrhea during visit 4. 

#find the number of infants that had vomiting in visit 4
visit4.vomit <- data.frame("mid"=redcap.wide$mid, "visit4.vomit"=rowSums(redcap.wide[,c(grep(pattern="f404_vomitn_q5_2$|f404_vomitn_q6_2$", colnames(redcap.wide)))], na.rm=T))
table(visit4.vomit$visit4.vomit>0)#twenty infants total vomited during visit 4. 

#find the number of infants that had fever in visit 4
visit4.fever <- data.frame("mid"=redcap.wide$mid, "visit4.fever"=rowSums(redcap.wide[,c(grep(pattern="f404_fever_q5_5$|f404_fever_q6_5$", colnames(redcap.wide)))], na.rm=T))
table(visit4.fever$visit4.fever>0)#fourty-two infants total had fever during visit 4. 

#merge the diarrhea, vomit, and fever morbidities
visit4.morbidities <- merge(visit4.diarrhea, visit4.fever, by="mid", all=F)
visit4.morbidities <- merge(visit4.morbidities, visit4.vomit, by="mid", all=F)

#now looking to see if infant MD280R (which had "See description" as an entry for
#specify illness) had diarrhea, fever or vomiting. 
visit4.morbidities[visit4.morbidities$mid=="MD280R",] #the infant had diarrhea, fever, and vomiting. 

#now looking to see if infant MD312P (which had "blood in feces" as entry for 
#illness specification) had diarrhea, fever, or vomiting. 
visit4.morbidities[visit4.morbidities$mid=="MD312P",] #the infant had fever and vomiting in this visit. 



#STEP 7: look at the medications infants received during the study, including supplements
#such as probiotic, laxative, and iron. 

#first, visit 2:
#antibiotic during visit 2
table(redcap.wide$f204_med_q9_2) 
#five infants received antibiotic (designated by "1"). Twelve infants received some 
#other type of medicine. Looking to see what the "other" medicine was..
unique(redcap.wide$f204_other_q9_2a)
#"brentan creme (miconazol) for fungus infection in diaper area", "pain killers",
#"Mycostatin (nystatin) against yeast/fungal infection (thrush)", "fungicide",
#"Painkillers", "Brantan", "mycostatin", "Fungal medication", "Steroid ointment, Brentacort"
#"Laxative", "Brentan", "Chloramphenicol ointment".Brentan is antifungal cream. Brentacort
#is ointment for eczema. None of these, except the laxative are notable for our analyses. 

#The laxative is notable and will be included in a figure of medications. 
redcap.wide[redcap.wide$f204_other_q9_2a=="Laxative","mid"] #MD293N

#other medications in visit 2
table(redcap.wide$f204_medtype_q9_3) #no more medications. 

#Other listings of medications during visit 2
unique(redcap.wide$f204_sptreatmt_q5_1_4) #none
unique(redcap.wide$f204_sptreatmt_q5_2_4) #none
unique(redcap.wide$f204_sptreatmt_q5_3_4) #none
unique(redcap.wide$f204_sptreatmt_q5_4_4) #none
unique(redcap.wide$f204_sptreatmt_q5_5_4) #none
unique(redcap.wide$f204_sptreatmt_q5_6_4) #none
unique(redcap.wide$f204_sptreatmt_q5_7_4) #none
unique(redcap.wide$f204_sptreatmt_q5_8_4) #mycostatin (sweet liquid 4x daily)
unique(redcap.wide$f204_sptreatmt_q6_1_1d) #none
unique(redcap.wide$f204_sptreatmt_q6_1_2d) #none
unique(redcap.wide$f204_sptreatmt_q6_1_3d) #none
unique(redcap.wide$f204_sptreatmt_q6_1_4d) #none
unique(redcap.wide$f204_sptreatmt_q6_2_1d) #none
unique(redcap.wide$f204_sptreatmt_q6_2_2d) #none
unique(redcap.wide$f204_sptreatmt_q6_2_3d) #none
unique(redcap.wide$f204_sptreatmt_q6_2_4d) #none
unique(redcap.wide$f204_sptreatmt_q6_3_1d) #none
unique(redcap.wide$f204_sptreatmt_q6_3_2d) #none
unique(redcap.wide$f204_sptreatmt_q6_3_3d) #none
unique(redcap.wide$f204_sptreatmt_q6_3_4d) #none
unique(redcap.wide$f204_sptreatmt_q6_4_1d) # extra oxygen due to infection
unique(redcap.wide$f204_sptreatmt_q6_4_2d) #none
unique(redcap.wide$f204_sptreatmt_q6_4_3d) #none
unique(redcap.wide$f204_sptreatmt_q6_4_4d) #none
unique(redcap.wide$f204_sptreatmt_q6_5_1d) #antibiotics due to infection
#look at the MID to see if this mother already had marked antibiotic treatment
#for her infant. 
redcap.wide[redcap.wide$f204_sptreatmt_q6_5_1d=="Antibiotics due to infection","f204_med_q9_2"] #Yes, the mother already answered yes to the antibiotic medication question.
unique(redcap.wide$f204_sptreatmt_q6_5_2d) #none
unique(redcap.wide$f204_sptreatmt_q6_5_3d) #none
unique(redcap.wide$f204_sptreatmt_q6_5_4d) #none
unique(redcap.wide$f204_sptreatmt_q6_6_1d) #none
unique(redcap.wide$f204_sptreatmt_q6_6_2d) #none
unique(redcap.wide$f204_sptreatmt_q6_6_3d) #none
unique(redcap.wide$f204_sptreatmt_q6_6_4d) #none
unique(redcap.wide$f204_sptreatmt_q6_7_1d) #none
unique(redcap.wide$f204_sptreatmt_q6_7_2d) #none
unique(redcap.wide$f204_sptreatmt_q6_7_3d) #none
unique(redcap.wide$f204_sptreatmt_q6_7_4d) #none
unique(redcap.wide$f204_sptreatmt_q6_8_1d) #mycostatin, treatment for 6 weeks, light bed for three days, eye drops (kloramfenikol), double light (matte + influence), chloramphenicol ointment. 
#looking for the type of treatment and/or illness for the infant that was recorded to have "Treatment for 6 weeks"
redcap.wide[redcap.wide$f204_sptreatmt_q6_8_1d=="Treatment for 6 weeks","f204_med_q9_2"] #no medication listed. 
redcap.wide[redcap.wide$f204_sptreatmt_q6_8_1d=="Treatment for 6 weeks","f204_other_q6_8a"] #hip dysplasia was the diagnosis

unique(redcap.wide$f204_sptreatmt_q6_8_2d) #none
unique(redcap.wide$f204_sptreatmt_q6_8_3d) #none
unique(redcap.wide$f204_sptreatmt_q6_8_4d) #none
unique(redcap.wide$f204_trtmnt_q8_1_1e) #"no treatment", "unknown infection and antibiotic from 11-2-2019 to 17-2-2019 and 3 days hospitalization", "hospitalized for a week and a half"  
#one of the infants had antibiotic, looking to see which infant and if this infant
#already had antibiotic treatment listed previously
redcap.wide[c(grep(pattern="Unknown infection. Antibiotic from",redcap.wide$f204_trtmnt_q8_1_1e)), "f204_med_q9_2"]
#yes, this mother already listed antibiotic treatment for the infant.  

unique(redcap.wide$f204_trtmnt_q8_1_2f) #none
unique(redcap.wide$f204_trtmnt_q8_1_3f) #none

#look at supplements such as probiotic, laxative, and iron.
table(redcap.wide$f204_supp_q10) #all mothers responded yes to a supplement given.
#looking specifically at multivitamin
table(redcap.wide$f204_supplst_q10_2) #none
#looking specifically at iron
table(redcap.wide$f204_supplst_q10_3) #none
#looking specifically at vitamin D
table(redcap.wide$f204_supplst_q10_4) #all infants received vitamin D. 
#looking specifically at probiotics
table(redcap.wide$f204_supplst_q10_5) #24 infants received probiotic.
#looking at other supplements not listed
table(redcap.wide$f204_supplst_q10_6) #one infant received some other supplement.
#looking at the supplement given. 
unique(redcap.wide$f204_other_q10_6a) #"vitamin K as injection"

#gathering the information on antibiotic use among infants in visit2. 
visit2.antibiotic <- data.frame("mid"=redcap.wide$mid, "visit2.antibiotic"=redcap.wide$f204_med_q9_2)
visit2.antibiotic[is.na(visit2.antibiotic$visit2.antibiotic),"visit2.antibiotic"] <- 0
visit2.antibiotic[visit2.antibiotic$visit2.antibiotic==3,"visit2.antibiotic"] <- 0
table(visit2.antibiotic$visit2.antibiotic>0)#five infants total had antibiotic

#gathering the information on laxative use among infants in visit 2, which only
#includes one infant, MD293N
visit2.laxative <- data.frame("mid"=redcap.wide$mid, "visit2.laxative"=0)
visit2.laxative[visit2.laxative$mid=="MD293N", "visit2.laxative"] <- 1
table(visit2.laxative$visit2.laxative>0) #one infant with laxative

#gathering the information on probiotic use among infants in visit 2
visit2.probiotic <- data.frame("mid"=redcap.wide$mid, "visit2.probiotic"=redcap.wide$f204_supplst_q10_5)
table(visit2.probiotic$visit2.probiotic>0) #24 infants had probiotic

#combine all the medication information from visit 2 in a dataframe
visit2.medication <- merge(visit2.antibiotic, visit2.laxative, by="mid", all=F)
visit2.medication <- merge(visit2.medication, visit2.probiotic, by="mid", all=F)


#now, visit 3:
#antibiotic during visit 3
table(redcap.wide$f304_med_q9_2) 
#two infants received antibiotic (designated by "1"). Eighteen infants received some 
#other type of medicine. Looking to see what the "other" medicine was..
unique(redcap.wide$f304_other_q9_2a)
#"brentan creme (miconazol) for fungus infection in diaper area", "Mycostatin", 
#"hormone cream", "Pain killers", "Diproderm creme", "Pain reliever", "pain killers",
#"paracetamol", "brentacort", "Panodil", "lactulose", "Paracetamol", "Lactulose",
#"For constipation, a few days", "Thrush", "panodil"

#paracetamol and panodil are pain relievers. 
#lactulose was a laxative used for several infants and this will be noted in our analyses
#so identifying which infants had laxative
redcap.wide[c(grep(pattern="lactulose", redcap.wide$f304_other_q9_2a, ignore.case = TRUE)),"mid"] #MD223Y and MD261C

#other medications in visit 3
table(redcap.wide$f304_medtype_q9_3) #four more medications. 
#looking at what the additional medications were
table(redcap.wide$f304_med_q9_4) #all were "other" (i.e. not antibiotic or asthma/allergy medicine)
unique(redcap.wide$f304_other_q9_4a) 
#"Brentan oral cavity gel", "pain killers", "pain killers (max dose) and bonyl (1 time as a shot also pain killers to reduce fever)",
#"anti-emetica"   
#none of these will be included in the analyses. 

table(redcap.wide$f304_medtype_q9_5) #no more medications

#Other listings of medications during visit 3
unique(redcap.wide$f304_sptreatmt_q5_1_4) #none
unique(redcap.wide$f304_sptreatmt_q5_2_4) #none
unique(redcap.wide$f304_sptreatmt_q5_3_4) #none
unique(redcap.wide$f304_sptreatmt_q5_4_4) #none
unique(redcap.wide$f304_sptreatmt_q5_5_4) #none
unique(redcap.wide$f304_sptreatmt_q5_6_4) #none
unique(redcap.wide$f304_sptreatmt_q5_7_4) #none
unique(redcap.wide$f304_sptreatmt_q5_8_4) #"fucithalmic 1% (fucidinsyre), 2 drops", "Mycostatin"  
#fucithalmic is an antibacterial eye drop
unique(redcap.wide$f304_sptreatmt_q6_1_1d) #none
unique(redcap.wide$f304_sptreatmt_q6_1_2d) #none
unique(redcap.wide$f304_sptreatmt_q6_1_3d) #none
unique(redcap.wide$f304_sptreatmt_q6_1_4d) #none
unique(redcap.wide$f304_sptreatmt_q6_2_1d) #"anti-emetica + painkillers"
unique(redcap.wide$f304_sptreatmt_q6_2_2d) #none
unique(redcap.wide$f304_sptreatmt_q6_2_3d) #none
unique(redcap.wide$f304_sptreatmt_q6_2_4d) #none
unique(redcap.wide$f304_sptreatmt_q6_3_1d) #none
unique(redcap.wide$f304_sptreatmt_q6_3_2d) #none
unique(redcap.wide$f304_sptreatmt_q6_3_3d) #none
unique(redcap.wide$f304_sptreatmt_q6_3_4d) #none
unique(redcap.wide$f304_sptreatmt_q6_4_1d) #"antibiotics and painkillers due to a virus (parecco) at 5 weeks old"
#the antibiotic treatment is notable, and it was listed under treatment for rapid breathing. 
#looking to see if the mother indicated specifically under antibiotic
#medication that her infant received antibiotics during visit 3. 
redcap.wide[c(grep(pattern="antibiotic", redcap.wide$f304_sptreatmt_q6_4_1d)), "f304_med_q9_2"]
#yes, the mother already indicated in the medication question that the infant received antibiotics
unique(redcap.wide$f304_sptreatmt_q6_4_2d) #none
unique(redcap.wide$f304_sptreatmt_q6_4_3d) #none
unique(redcap.wide$f304_sptreatmt_q6_4_4d) #none
unique(redcap.wide$f304_sptreatmt_q6_5_1d) #"no treatment, no diagnose", "same as 6.4.1d", "paracetamol", "Panodil"
#these treatments were listed under fever. Looking to see what response is meant by
#the answer "same as 6.4.1d". 
redcap.wide[c(grep(pattern="same as 6[.]4", redcap.wide$f304_sptreatmt_q6_5_1d)),"f304_sptreatmt_q6_4_1d"]
#this matches  to the "antibiotics and painkillers due to a virus (parecco) at 5 weeks old". 
unique(redcap.wide$f304_sptreatmt_q6_5_2d) #none
unique(redcap.wide$f304_sptreatmt_q6_5_3d) #none
unique(redcap.wide$f304_sptreatmt_q6_5_4d) #none
unique(redcap.wide$f304_sptreatmt_q6_6_1d) #none
unique(redcap.wide$f304_sptreatmt_q6_6_2d) #none
unique(redcap.wide$f304_sptreatmt_q6_6_3d) #none
unique(redcap.wide$f304_sptreatmt_q6_6_4d) #none
unique(redcap.wide$f304_sptreatmt_q6_7_1d) #"Penicillin 7 days"
#looking to see if this mother already answered the earlier medication question
#for the infant having antibiotic
redcap.wide[c(grep(pattern="Penicillin", redcap.wide$f304_sptreatmt_q6_7_1d)),"f304_med_q9_2"]
#yes, she answered that the infant had antibiotics
unique(redcap.wide$f304_sptreatmt_q6_7_2d) #none
unique(redcap.wide$f304_sptreatmt_q6_7_3d) #none
unique(redcap.wide$f304_sptreatmt_q6_7_4d) #none
unique(redcap.wide$f304_sptreatmt_q6_8_1d) #"Mycostatin", "Diproderm, 0.05 % creme", "brentacort (creme twice a day)", 
#"cream (don't know what was in it - happened in Thailand)", "mother received Brentan"       
unique(redcap.wide$f304_sptreatmt_q6_8_2d) #none
unique(redcap.wide$f304_sptreatmt_q6_8_3d) #none
unique(redcap.wide$f304_sptreatmt_q6_8_4d) #none
unique(redcap.wide$f304_trtmnt_q8_1_1e) #"antibiotics (one time at the hospital) and pain killers. the full episode was 6-7 days but hospitalization was 3 days."  
#one of the infants had antibiotic, looking to see which infant and if this infant
#already had antibiotic treatment listed previously
redcap.wide[c(grep(pattern="antibiotics",redcap.wide$f304_trtmnt_q8_1_1e)), "f304_med_q9_2"]
#yes, the mother already indicated that the infant had antibiotic

unique(redcap.wide$f304_trtmnt_q8_1_2f) #none
unique(redcap.wide$f304_trtmnt_q8_1_3f) #none

#look at supplements such as probiotic, laxative, and iron.
table(redcap.wide$f304_supp_q10) #all mothers responded yes to a supplement given.
#looking specifically at multivitamin
table(redcap.wide$f304_supplst_q10_2) #none
#looking specifically at iron
table(redcap.wide$f304_supplst_q10_3) #none
#looking specifically at vitamin D
table(redcap.wide$f304_supplst_q10_4) #all infants received vitamin D. 
#looking specifically at probiotics
table(redcap.wide$f304_supplst_q10_5) #18 infants received probiotic.
#looking at other supplements not listed
table(redcap.wide$f304_supplst_q10_6) #no infants received some other supplement.

#gathering the information on antibiotic use among infants in visit 3. 
visit3.antibiotic <- data.frame("mid"=redcap.wide$mid, "visit3.antibiotic"=redcap.wide$f304_med_q9_2)
visit3.antibiotic[is.na(visit3.antibiotic$visit3.antibiotic),"visit3.antibiotic"] <- 0
visit3.antibiotic[visit3.antibiotic$visit3.antibiotic==3,"visit3.antibiotic"] <- 0
table(visit3.antibiotic$visit3.antibiotic>0)#two infants total had antibiotic
visit3.antibiotic[visit3.antibiotic$visit3.antibiotic==1,"mid"] #MD070Y and MD281K are listed here and these are the only two infants that had antibiotics recorded. 

#gathering the information on laxative use among infants in visit 3, which
#includes two infants, MD223Y and MD261C
visit3.laxative <- data.frame("mid"=redcap.wide$mid, "visit3.laxative"=0)
visit3.laxative[visit3.laxative$mid=="MD223Y"|visit3.laxative$mid=="MD261C", "visit3.laxative"] <- 1
table(visit3.laxative$visit3.laxative>0) #two infants with laxative

#gathering the information on probiotic use among infants in visit 3
visit3.probiotic <- data.frame("mid"=redcap.wide$mid, "visit3.probiotic"=redcap.wide$f304_supplst_q10_5)
table(visit3.probiotic$visit3.probiotic>0) #18 infants had probiotic

#combine all the medication information from visit 3 in a dataframe
visit3.medication <- merge(visit3.antibiotic, visit3.laxative, by="mid", all=F)
visit3.medication <- merge(visit3.medication, visit3.probiotic, by="mid", all=F)


#now, visit 4:
#antibiotic during visit 4
table(redcap.wide$f404_med_q9_2) 
#three infants received antibiotic (designated by "1"). Twenty-nine infants received some 
#other type of medicine. Looking to see what the "other" medicine was..
unique(redcap.wide$f404_other_q9_2a)
#"Lactulose", "Fungus infection in the mouth", "Children's panodil", "fine", 
#"Painkiller", "paracetamol", "Antipyretic", "Pain killers", "Prophylactic panodil"
#"klyx", "laktulose", "panodol", "Painkillers", "lactulose", "Panodil", "movicol", "lactolose"
#"Panodil Junior", "children's panodil", "panodil", "Junior Panodil", "penodil",
#"panodil 2 times", "Fungal cream for tongue", "Lactolose for constipation" 

#Movicol is a laxative. The laxative is notable and will be included in a figure of medications. 
redcap.wide[c(grep(pattern="movicol|la[ck]t[uo]lose", ignore.case = TRUE,redcap.wide$f404_other_q9_2a)),"mid"] 
#"MD007M","MD061A", "MD146K", "MD223Y", "MD261C", "MD273F", "MD378X"

#other medications in visit 4
table(redcap.wide$f404_medtype_q9_3) #eight more medications
table(redcap.wide$f404_med_q9_4) #all eight were specified as "other"
unique(redcap.wide$f404_other_q9_4a)
#"Painkiller", "paracetamol", "Paracetamol", "Lactulose", 
#"cream against viral eczema", "Cream for skin rash, Inotyol"

#the use of lactulose is notable so taking note of the infants
redcap.wide[c(grep(pattern="lactulose", ignore.case=T, redcap.wide$f404_other_q9_4a)), "mid"]
# "MD261C"

#other medications in visit 4
table(redcap.wide$f404_medtype_q9_5) #none

#Other listings of medications during visit 4
unique(redcap.wide$f404_sptreatmt_q5_1_4) #none
unique(redcap.wide$f404_sptreatmt_q5_2_4) #none
unique(redcap.wide$f404_sptreatmt_q5_3_4) #"Panodil"
unique(redcap.wide$f404_sptreatmt_q5_4_4) #"Antibiotics", "panodil"
#the antibiotic use is notable so taking note of the infant
redcap.wide[c(grep(pattern="biotics", ignore.case = T, redcap.wide$f404_sptreatmt_q5_4_4)), "mid"]
#MD188M
unique(redcap.wide$f404_sptreatmt_q5_5_4) #"Antibiotics", "Panodil"
#the antibiotic use is notable so taking note of the infant
redcap.wide[c(grep(pattern="Antibiotics", ignore.case = T, redcap.wide$f404_sptreatmt_q5_5_4)), "mid"]
#MD188M
unique(redcap.wide$f404_sptreatmt_q5_6_4) #none
unique(redcap.wide$f404_sptreatmt_q5_7_4) #"Antibiotics"
#the antibiotic use is notable so taking note of the infant
redcap.wide[c(grep(pattern="Antibiotics", ignore.case = T, redcap.wide$f404_sptreatmt_q5_7_4)), "mid"]
#MD188M
unique(redcap.wide$f404_sptreatmt_q5_8_4) #"Lactulose", "Brentan", "creme"
#the lactulose (laxative) use is notable so taking note of the infant
redcap.wide[c(grep(pattern="Lactulose", ignore.case = T, redcap.wide$f404_sptreatmt_q5_8_4)), "mid"]
#MD007M
unique(redcap.wide$f404_sptreatmt_q6_1_1d) #none
unique(redcap.wide$f404_sptreatmt_q6_1_2d) #none
unique(redcap.wide$f404_sptreatmt_q6_1_3d) #none
unique(redcap.wide$f404_sptreatmt_q6_1_4d) #none
unique(redcap.wide$f404_sptreatmt_q6_2_1d) #"klyx"
unique(redcap.wide$f404_sptreatmt_q6_2_2d) #none
unique(redcap.wide$f404_sptreatmt_q6_2_3d) #none
unique(redcap.wide$f404_sptreatmt_q6_2_4d) #none
unique(redcap.wide$f404_sptreatmt_q6_3_1d) #none
unique(redcap.wide$f404_sptreatmt_q6_3_2d) #none
unique(redcap.wide$f404_sptreatmt_q6_3_3d) #none
unique(redcap.wide$f404_sptreatmt_q6_3_4d) #none
unique(redcap.wide$f404_sptreatmt_q6_4_1d) #none
unique(redcap.wide$f404_sptreatmt_q6_4_2d) #none
unique(redcap.wide$f404_sptreatmt_q6_4_3d) #none
unique(redcap.wide$f404_sptreatmt_q6_4_4d) #none
unique(redcap.wide$f404_sptreatmt_q6_5_1d) #"panodil"
unique(redcap.wide$f404_sptreatmt_q6_5_2d) #none
unique(redcap.wide$f404_sptreatmt_q6_5_3d) #none
unique(redcap.wide$f404_sptreatmt_q6_5_4d) #none
unique(redcap.wide$f404_sptreatmt_q6_6_1d) #none
unique(redcap.wide$f404_sptreatmt_q6_6_2d) #none
unique(redcap.wide$f404_sptreatmt_q6_6_3d) #none
unique(redcap.wide$f404_sptreatmt_q6_6_4d) #none
unique(redcap.wide$f404_sptreatmt_q6_7_1d) #none
unique(redcap.wide$f404_sptreatmt_q6_7_2d) #none
unique(redcap.wide$f404_sptreatmt_q6_7_3d) #none
unique(redcap.wide$f404_sptreatmt_q6_7_4d) #none
unique(redcap.wide$f404_sptreatmt_q6_8_1d) #"Lactulose", "Lactofarm junior"
#these are both laxatives and will be included in our analyses
redcap.wide[c(grep(pattern="lactulose|lactofarm", ignore.case=T, redcap.wide$f404_sptreatmt_q6_8_1d)), "mid"] #MD050Q, MD090G
unique(redcap.wide$f404_sptreatmt_q6_8_2d) #none
unique(redcap.wide$f404_sptreatmt_q6_8_3d) #none
unique(redcap.wide$f404_sptreatmt_q6_8_4d) #none
unique(redcap.wide$f404_trtmnt_q8_1_1e) 
#"Q6.1.4: infant had fever and was seen by general practitioner who sent the infant to the hospital; different doctors took blood and urine samples and assumed it was either throat infection, kidney infection and/or diarrhea. antibiotics was prescribed, but was not given to the infant and the illness passed."
#"Åben indlæggelse" is another response which means "open admission"  

unique(redcap.wide$f404_trtmnt_q8_1_2f) #none
unique(redcap.wide$f404_trtmnt_q8_1_3f) #none

#look at supplements such as probiotic, laxative, and iron.
table(redcap.wide$f404_supp_q10) #all mothers responded yes to a supplement given.
#looking specifically at multivitamin
table(redcap.wide$f404_supplst_q10_2) #none
#looking specifically at iron
table(redcap.wide$f404_supplst_q10_3) #none
#looking specifically at vitamin D
table(redcap.wide$f404_supplst_q10_4) #all infants received vitamin D. 
#looking specifically at probiotics
table(redcap.wide$f404_supplst_q10_5) #11 infants received probiotic.
#looking at other supplements not listed
table(redcap.wide$f404_supplst_q10_6) #one infant received some other supplement.
#looking at the supplement given. 
unique(redcap.wide$f404_other_q10_6a) #"Biogaia" 
#Biogaia is a probiotic which is included in our analyses so seeing if this mother
#already recorded this infant as having probiotic in the earlier question. 
redcap.wide[c(grep(pattern="Biogaia", ignore.case=T, redcap.wide$f404_other_q10_6a)),"f404_supplst_q10_5"]
#yes this mother already answered that the infant had probiotic.

#gathering the information on antibiotic use among infants in visit4. 
visit4.antibiotic <- data.frame("mid"=redcap.wide$mid, "visit4.antibiotic"=redcap.wide$f404_med_q9_2)
visit4.antibiotic[is.na(visit4.antibiotic$visit4.antibiotic),"visit4.antibiotic"] <- 0
visit4.antibiotic[visit4.antibiotic$visit4.antibiotic==3,"visit4.antibiotic"] <- 0
table(visit4.antibiotic$visit4.antibiotic>0)#three infants total had antibiotic

visit4.antibiotic[visit4.antibiotic$visit4.antibiotic>0,"mid"] #MD080P, MD188M, MD210C had antibiotic marked already as in treatment list. 

#gathering the information on laxative use among infants in visit 4, which includes
#"MD007M","MD061A", "MD146K", "MD223Y", "MD261C", "MD273F", "MD349R", "MD378X", "MD050Q", "MD090G"

visit4.laxative <- data.frame("mid"=redcap.wide$mid, "visit4.laxative"=0)
visit4.laxative[c(grep(pattern="MD007M|MD061A|MD146K|MD223Y|MD261C|MD273F|MD349R|MD378X|MD050Q|MD090G",visit4.laxative$mid)), "visit4.laxative"] <- 1
table(visit4.laxative$visit4.laxative>0) #ten infants with laxative

#gathering the information on probiotic use among infants in visit 4
visit4.probiotic <- data.frame("mid"=redcap.wide$mid, "visit4.probiotic"=redcap.wide$f404_supplst_q10_5)
table(visit4.probiotic$visit4.probiotic>0) #11 infants had probiotic

#combine all the medication information in a data frame
visit4.medication <- merge(visit4.antibiotic, visit4.laxative, by="mid", all=F)
visit4.medication <- merge(visit4.medication, visit4.probiotic, by="mid", all=F)

#STEP 8: combine the morbidity and medication data together in a table to make a figure. 

df.illness.medx <- merge(visit2.morbidities, visit3.morbidities, by="mid", all=F)
df.illness.medx <- merge(df.illness.medx, visit4.morbidities, by="mid", all=F)
df.illness.medx <- merge(df.illness.medx, visit2.medication, by="mid", all=F)
df.illness.medx <- merge(df.illness.medx, visit3.medication, by="mid", all=F)
df.illness.medx <- merge(df.illness.medx, visit4.medication, by="mid", all=F)

df.illness.medx[df.illness.medx$visit2.antibiotic>0,"visit2.antibiotic"] <- "A"
df.illness.medx[df.illness.medx$visit3.antibiotic>0,"visit3.antibiotic"] <- "A"
df.illness.medx[df.illness.medx$visit4.antibiotic>0,"visit4.antibiotic"] <- "A"
df.illness.medx[df.illness.medx$visit2.probiotic>0,"visit2.probiotic"] <- "P"
df.illness.medx[df.illness.medx$visit3.probiotic>0,"visit3.probiotic"] <- "P"
df.illness.medx[df.illness.medx$visit4.probiotic>0,"visit4.probiotic"] <- "P"
df.illness.medx[df.illness.medx$visit2.laxative>0,"visit2.laxative"] <- "L"
df.illness.medx[df.illness.medx$visit3.laxative>0,"visit3.laxative"] <- "L"
df.illness.medx[df.illness.medx$visit4.laxative>0,"visit4.laxative"] <- "L"

df.illness.medx$visit2.medication <- paste(df.illness.medx$visit2.antibiotic, df.illness.medx$visit2.probiotic, df.illness.medx$visit2.laxative, sep=",")
df.illness.medx$visit2.medication <- gsub(pattern="0,|,0$|^0|0", replacement="", df.illness.medx$visit2.medication)
df.illness.medx$visit2.medication <- gsub(pattern="A,$", replacement="A", df.illness.medx$visit2.medication)

df.illness.medx$visit3.medication <- paste(df.illness.medx$visit3.antibiotic, df.illness.medx$visit3.probiotic, df.illness.medx$visit3.laxative, sep=",")
df.illness.medx$visit3.medication <- gsub(pattern="0,|,0$|^0|0", replacement="", df.illness.medx$visit3.medication)
df.illness.medx$visit3.medication <- gsub(pattern="A,$", replacement="A", df.illness.medx$visit3.medication)

df.illness.medx$visit4.medication <- paste(df.illness.medx$visit4.antibiotic, df.illness.medx$visit4.probiotic, df.illness.medx$visit4.laxative, sep=",")
df.illness.medx$visit4.medication <- gsub(pattern="0,|,0$|^0|0", replacement="", df.illness.medx$visit4.medication)
df.illness.medx$visit4.medication <- gsub(pattern="A,$", replacement="A", df.illness.medx$visit4.medication)

df.medx <- df.illness.medx[,c(1,20:22)]
df.medx <- pivot_longer(df.medx, cols=2:4)
colnames(df.medx) <- c("mid", "visit", "drug")
df.medx$visit <- gsub(pattern="visit", replacement="", df.medx$visit)
df.medx$visit <- gsub(pattern="[.]medication", replacement="", df.medx$visit)

df.illness.medx.v2 <- pivot_longer(df.illness.medx, cols=2:10)

df.illness.medx.v2[is.na(df.illness.medx.v2$value) ,"value"] <- 0
df.illness.medx.v2[df.illness.medx.v2$value>0,"value"] <- 1

df.illness.medx.v2$visit <- gsub(pattern="visit", replacement="", df.illness.medx.v2$name)
df.illness.medx.v2$visit <- gsub(pattern="[.]diarrhea$|[.]fever$|[.]vomit$", replacement="", df.illness.medx.v2$visit)

df.illness.medx.v2$morbidity <- gsub(pattern="visit[2-4][.]", replacement="", df.illness.medx.v2$name)
df.illness.medx.v2$cid <- gsub(pattern="^M", replacement="C", df.illness.medx.v2$mid)

df.illness.medx.v3 <- merge(df.illness.medx.v2[,c(1,14:18)], df.medx, by=c("mid", "visit"), all=F)
total.illness.per.mid <- data.frame(tapply(df.illness.medx.v3$value, INDEX=df.illness.medx.v3$mid, FUN=sum))
total.illness.per.mid$mid <- rownames(total.illness.per.mid)
colnames(total.illness.per.mid) <- c("total.illness", "mid")

df.illness.medx.v4 <- merge(df.illness.medx.v3, total.illness.per.mid, by="mid", all=F)

#change the column "name" from visit # + morbidity to morbidity + morbidity rate (%)
morbidity.rate.per.visit <- data.frame(tapply(X=df.illness.medx.v4$value, INDEX=df.illness.medx.v4$name, FUN=sum))
morbidity.rate.per.visit$name <- rownames(morbidity.rate.per.visit)
colnames(morbidity.rate.per.visit) <- c("total", "name")
morbidity.rate.per.visit$percentage <- round(((morbidity.rate.per.visit$total/109)*100),digits=1)
morbidity.rate.per.visit$percentage <- paste(morbidity.rate.per.visit$percentage, " %", sep="")
morbidity.rate.per.visit$visit.morbidity.rate <- paste(morbidity.rate.per.visit$name, morbidity.rate.per.visit$percentage, sep=", ")
morbidity.rate.per.visit$visit.morbidity.rate <- gsub(pattern="visit[2-4][.]", replacement="", morbidity.rate.per.visit$visit.morbidity.rate)

df.illness.medx.v5 <- merge(x=df.illness.medx.v4, y=morbidity.rate.per.visit[,c(2,4)], by="name", all=F)
library(forcats)
ggplot(df.illness.medx.v5, aes(x=visit.morbidity.rate, y=fct_reorder(cid,total.illness), fill=morbidity, alpha=value)) + geom_tile()+scale_alpha(range=c(0,1))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(text=element_text(size=8)) +geom_text(alpha=1,data=df.illness.medx.v5,aes(label=drug), color="black", size=1.5)+facet_wrap(~visit, scales="free_x")+ylab("infant ID")+xlab("")+guides(alpha="none")
ggsave("morbidity.medication.per.infant.jpeg", plot=last_plot(), dpi=600, height=10, width=4, units="in")

#we cannot display the infant IDs so need to remove them and in place have the 
#birth mode and infant gender listed. 
df.illness.medx.v6 <- merge(x=df.illness.medx.v5, y=redcap.wide[,c(grep(pattern="^mid|f201_infgend_q10$|f201_mod_q11$", colnames(redcap.wide)))], by="mid", all=F)
df.illness.medx.v6$f201_infgend_q10 <- gsub(pattern="1", replacement="boy", df.illness.medx.v6$f201_infgend_q10)
df.illness.medx.v6$f201_infgend_q10 <- gsub(pattern="2", replacement="girl", df.illness.medx.v6$f201_infgend_q10)
df.illness.medx.v6$f201_mod_q11 <- gsub(pattern="1", replacement="vaginal", df.illness.medx.v6$f201_mod_q11)
df.illness.medx.v6$f201_mod_q11 <- gsub(pattern="2|3", replacement="c-section", df.illness.medx.v6$f201_mod_q11)
df.illness.medx.v6$infant.gender.mod <- paste(df.illness.medx.v6$f201_infgend_q10, df.illness.medx.v6$f201_mod_q11, sep=", ")


ggplot(df.illness.medx.v6, aes(x=visit.morbidity.rate, y=fct_reorder(cid,total.illness), fill=morbidity, alpha=value)) + geom_tile()+scale_alpha(range=c(0,1))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(text=element_text(size=8)) +geom_text(alpha=1,data=df.illness.medx.v6,aes(label=drug), color="black", size=1.5)+facet_wrap(~visit, scales="free_x")+ylab("infant gender, mode of delivery")+xlab("")+guides(alpha="none")+scale_y_discrete(breaks=c(df.illness.medx.v6$cid),labels=c(df.illness.medx.v6$infant.gender.mod))
ggsave("morbidity.medication.per.infant.noCID.jpeg", plot=last_plot(), dpi=600, height=10, width=4, units="in")

#STEP 9: adjusting the format of the data frame to make into table that can be 
#merged
visit2.info <- merge(visit2.morbidities, visit2.medication, by="mid", all=F) 
visit2.info$visit <- 2
colnames(visit2.info)
colnames(visit2.info) <- c("mid", "infant.diarrhea", "infant.fever", "infant.vomit", "infant.antibiotic", "infant.laxative", "infant.probiotic", "visit")
colnames(visit2.info)
visit2.info[visit2.info$infant.diarrhea > 0, "infant.diarrhea"] <- "yes"
visit2.info[visit2.info$infant.diarrhea=="0", "infant.diarrhea"] <- "no"
table(visit2.info$infant.diarrhea)

visit2.info[visit2.info$infant.fever > 0, "infant.fever"] <- "yes"
visit2.info[visit2.info$infant.fever=="0", "infant.fever"] <- "no"
table(visit2.info$infant.fever)

visit2.info[visit2.info$infant.vomit > 0, "infant.vomit"] <- "yes"
visit2.info[visit2.info$infant.vomit=="0", "infant.vomit"] <- "no"
table(visit2.info$infant.vomit)

visit2.info[visit2.info$infant.antibiotic > 0, "infant.antibiotic"] <- "yes"
visit2.info[visit2.info$infant.antibiotic=="0", "infant.antibiotic"] <- "no"
table(visit2.info$infant.antibiotic)

visit2.info[visit2.info$infant.laxative > 0, "infant.laxative"] <- "yes"
visit2.info[visit2.info$infant.laxative=="0", "infant.laxative"] <- "no"
table(visit2.info$infant.laxative)

visit2.info[visit2.info$infant.probiotic > 0, "infant.probiotic"] <- "yes"
visit2.info[visit2.info$infant.probiotic=="0", "infant.probiotic"] <- "no"
table(visit2.info$infant.probiotic)

visit3.info <- merge(visit3.morbidities, visit3.medication, by="mid", all=F) 
visit3.info$visit <- 3
colnames(visit3.info)
colnames(visit3.info) <- c("mid", "infant.diarrhea", "infant.fever", "infant.vomit", "infant.antibiotic", "infant.laxative", "infant.probiotic", "visit")
colnames(visit3.info)
visit3.info[visit3.info$infant.diarrhea > 0, "infant.diarrhea"] <- "yes"
visit3.info[visit3.info$infant.diarrhea=="0", "infant.diarrhea"] <- "no"
table(visit3.info$infant.diarrhea)

visit3.info[visit3.info$infant.fever > 0, "infant.fever"] <- "yes"
visit3.info[visit3.info$infant.fever=="0", "infant.fever"] <- "no"
table(visit3.info$infant.fever)

visit3.info[visit3.info$infant.vomit > 0, "infant.vomit"] <- "yes"
visit3.info[visit3.info$infant.vomit=="0", "infant.vomit"] <- "no"
table(visit3.info$infant.vomit)

visit3.info[visit3.info$infant.antibiotic > 0, "infant.antibiotic"] <- "yes"
visit3.info[visit3.info$infant.antibiotic=="0", "infant.antibiotic"] <- "no"
table(visit3.info$infant.antibiotic)

visit3.info[visit3.info$infant.laxative > 0, "infant.laxative"] <- "yes"
visit3.info[visit3.info$infant.laxative=="0", "infant.laxative"] <- "no"
table(visit3.info$infant.laxative)

visit3.info[visit3.info$infant.probiotic > 0, "infant.probiotic"] <- "yes"
visit3.info[visit3.info$infant.probiotic=="0", "infant.probiotic"] <- "no"
table(visit3.info$infant.probiotic)

visit4.info <- merge(visit4.morbidities, visit4.medication, by="mid", all=F) 
visit4.info$visit <- 4
colnames(visit4.info)
colnames(visit4.info) <- c("mid", "infant.diarrhea", "infant.fever", "infant.vomit", "infant.antibiotic", "infant.laxative", "infant.probiotic", "visit")
colnames(visit4.info)
visit4.info[visit4.info$infant.diarrhea > 0, "infant.diarrhea"] <- "yes"
visit4.info[visit4.info$infant.diarrhea=="0", "infant.diarrhea"] <- "no"
table(visit4.info$infant.diarrhea)

visit4.info[visit4.info$infant.fever > 0, "infant.fever"] <- "yes"
visit4.info[visit4.info$infant.fever=="0", "infant.fever"] <- "no"
table(visit4.info$infant.fever)

visit4.info[visit4.info$infant.vomit > 0, "infant.vomit"] <- "yes"
visit4.info[visit4.info$infant.vomit=="0", "infant.vomit"] <- "no"
table(visit4.info$infant.vomit)

visit4.info[visit4.info$infant.antibiotic > 0, "infant.antibiotic"] <- "yes"
visit4.info[visit4.info$infant.antibiotic=="0", "infant.antibiotic"] <- "no"
table(visit4.info$infant.antibiotic)

visit4.info[visit4.info$infant.laxative > 0, "infant.laxative"] <- "yes"
visit4.info[visit4.info$infant.laxative=="0", "infant.laxative"] <- "no"
table(visit4.info$infant.laxative)

visit4.info[visit4.info$infant.probiotic > 0, "infant.probiotic"] <- "yes"
visit4.info[visit4.info$infant.probiotic=="0", "infant.probiotic"] <- "no"
table(visit4.info$infant.probiotic)

medication.morbidity.table <- rbind(visit2.info, visit3.info, visit4.info)
metadata <- metadata[, -1]
metadata.v2 <- merge(x=metadata, y=medication.morbidity.table, by=c("mid", "visit"), all=T)

#STEP 10: make columns to indicate if the child EVER (in any of the three visits) had:
#diarrhea, fever, and/or vomit
infant.diarrhea.in.any.visit <- unique(metadata.v2[metadata.v2$infant.diarrhea=="yes","mid"])
metadata.v2$infant.diarrhea.ever.in.study <- "no"
metadata.v2[c(which(metadata.v2$mid %in% infant.diarrhea.in.any.visit)),"infant.diarrhea.ever.in.study"] <- "yes"
length(infant.diarrhea.in.any.visit)
length(unique(metadata.v2[metadata.v2$infant.diarrhea.ever.in.study=="yes","mid"]))

infant.fever.in.any.visit <- unique(metadata.v2[metadata.v2$infant.fever=="yes","mid"])
metadata.v2$infant.fever.ever.in.study <- "no"
metadata.v2[c(which(metadata.v2$mid %in% infant.fever.in.any.visit)),"infant.fever.ever.in.study"] <- "yes"
length(infant.fever.in.any.visit)
length(unique(metadata.v2[metadata.v2$infant.fever.ever.in.study=="yes","mid"]))

infant.vomit.in.any.visit <- unique(metadata.v2[metadata.v2$infant.vomit=="yes","mid"])
metadata.v2$infant.vomit.ever.in.study <- "no"
metadata.v2[c(which(metadata.v2$mid %in% infant.vomit.in.any.visit)),"infant.vomit.ever.in.study"] <- "yes"
length(infant.vomit.in.any.visit)
length(unique(metadata.v2[metadata.v2$infant.vomit.ever.in.study=="yes","mid"]))

infant.sick.in.any.visit <- c(infant.diarrhea.in.any.visit, infant.fever.in.any.visit, infant.vomit.in.any.visit)
infant.sick.in.any.visit <- unique(infant.sick.in.any.visit)
metadata.v2$infant.sick.ever.in.study <- "no"
metadata.v2[c(which(metadata.v2$mid %in% infant.sick.in.any.visit)),"infant.sick.ever.in.study"] <- "yes"
length(infant.sick.in.any.visit)
length(unique(metadata.v2[metadata.v2$infant.sick.ever.in.study=="yes","mid"]))

#STEP 11: make columns for child's morbidity status in each visit. For example,
#if a child had diarrhea in visit 3, then all stool samples (from visits 2, 3 and 4)
#for that child will have the value "yes" for diarrhea in visit 3. 

#for diarrhea
metadata.v2$infant.diarrhea.in.visit2 <- "no"
metadata.v2[c(which(metadata.v2$mid %in% (visit2.diarrhea[visit2.diarrhea$visit2.diarrhea>0,][["mid"]]))), "infant.diarrhea.in.visit2"] <- "yes"
length(unique(metadata.v2[metadata.v2$infant.diarrhea.in.visit2=="yes",][["mid"]]))
table(visit2.diarrhea$visit2.diarrhea)
identical(x=unique(metadata.v2[metadata.v2$infant.diarrhea.in.visit2=="yes",][["mid"]]), y=unique(visit2.diarrhea[visit2.diarrhea$visit2.diarrhea>0,"mid"]))

metadata.v2$infant.diarrhea.in.visit3 <- "no"
metadata.v2[c(which(metadata.v2$mid %in% (visit3.diarrhea[visit3.diarrhea$visit3.diarrhea>0,][["mid"]]))), "infant.diarrhea.in.visit3"] <- "yes"
length(unique(metadata.v2[metadata.v2$infant.diarrhea.in.visit3=="yes",][["mid"]]))
table(visit3.diarrhea$visit3.diarrhea)
identical(x=unique(metadata.v2[metadata.v2$infant.diarrhea.in.visit3=="yes",][["mid"]]), y=unique(visit3.diarrhea[visit3.diarrhea$visit3.diarrhea>0,"mid"]))

metadata.v2$infant.diarrhea.in.visit4 <- "no"
metadata.v2[c(which(metadata.v2$mid %in% (visit4.diarrhea[visit4.diarrhea$visit4.diarrhea>0,][["mid"]]))), "infant.diarrhea.in.visit4"] <- "yes"
length(unique(metadata.v2[metadata.v2$infant.diarrhea.in.visit4=="yes",][["mid"]]))
table(visit4.diarrhea$visit4.diarrhea)
identical(x=unique(metadata.v2[metadata.v2$infant.diarrhea.in.visit4=="yes",][["mid"]]), y=unique(visit4.diarrhea[visit4.diarrhea$visit4.diarrhea>0,"mid"]))

#a column for diarrhea in visit 3 and/or visit 4
metadata.v2$infant.diarrhea.in.visit3.or.visit4 <- "no"
metadata.v2[metadata.v2$infant.diarrhea.in.visit3=="yes", "infant.diarrhea.in.visit3.or.visit4"] <- "yes"
metadata.v2[metadata.v2$infant.diarrhea.in.visit4=="yes", "infant.diarrhea.in.visit3.or.visit4"] <- "yes"
length(unique(metadata.v2[metadata.v2$infant.diarrhea.in.visit3.or.visit4=="yes",][["mid"]]))


#for fever
metadata.v2$infant.fever.in.visit2 <- "no"
metadata.v2[c(which(metadata.v2$mid %in% (visit2.fever[visit2.fever$visit2.fever>0,][["mid"]]))), "infant.fever.in.visit2"] <- "yes"
length(unique(metadata.v2[metadata.v2$infant.fever.in.visit2=="yes",][["mid"]]))
table(visit2.fever$visit2.fever)
identical(x=unique(metadata.v2[metadata.v2$infant.fever.in.visit2=="yes",][["mid"]]), y=unique(visit2.fever[visit2.fever$visit2.fever>0,"mid"]))

metadata.v2$infant.fever.in.visit3 <- "no"
metadata.v2[c(which(metadata.v2$mid %in% (visit3.fever[visit3.fever$visit3.fever>0,][["mid"]]))), "infant.fever.in.visit3"] <- "yes"
length(unique(metadata.v2[metadata.v2$infant.fever.in.visit3=="yes",][["mid"]]))
table(visit3.fever$visit3.fever)
identical(x=unique(metadata.v2[metadata.v2$infant.fever.in.visit3=="yes",][["mid"]]), y=unique(visit3.fever[visit3.fever$visit3.fever>0,"mid"]))

metadata.v2$infant.fever.in.visit4 <- "no"
metadata.v2[c(which(metadata.v2$mid %in% (visit4.fever[visit4.fever$visit4.fever>0,][["mid"]]))), "infant.fever.in.visit4"] <- "yes"
length(unique(metadata.v2[metadata.v2$infant.fever.in.visit4=="yes",][["mid"]]))
table(visit4.fever$visit4.fever)
identical(x=unique(metadata.v2[metadata.v2$infant.fever.in.visit4=="yes",][["mid"]]), y=unique(visit4.fever[visit4.fever$visit4.fever>0,"mid"]))

#a column for fever in visit 3 and/or visit 4
metadata.v2$infant.fever.in.visit3.or.visit4 <- "no"
metadata.v2[metadata.v2$infant.fever.in.visit3=="yes", "infant.fever.in.visit3.or.visit4"] <- "yes"
metadata.v2[metadata.v2$infant.fever.in.visit4=="yes", "infant.fever.in.visit3.or.visit4"] <- "yes"
length(unique(metadata.v2[metadata.v2$infant.fever.in.visit3.or.visit4=="yes",][["mid"]]))

#for vomit
metadata.v2$infant.vomit.in.visit2 <- "no"
metadata.v2[c(which(metadata.v2$mid %in% (visit2.vomit[visit2.vomit$visit2.vomit>0,][["mid"]]))), "infant.vomit.in.visit2"] <- "yes"
length(unique(metadata.v2[metadata.v2$infant.vomit.in.visit2=="yes",][["mid"]]))
table(visit2.vomit$visit2.vomit)
identical(x=unique(metadata.v2[metadata.v2$infant.vomit.in.visit2=="yes",][["mid"]]), y=unique(visit2.vomit[visit2.vomit$visit2.vomit>0,"mid"]))

metadata.v2$infant.vomit.in.visit3 <- "no"
metadata.v2[c(which(metadata.v2$mid %in% (visit3.vomit[visit3.vomit$visit3.vomit>0,][["mid"]]))), "infant.vomit.in.visit3"] <- "yes"
length(unique(metadata.v2[metadata.v2$infant.vomit.in.visit3=="yes",][["mid"]]))
table(visit3.vomit$visit3.vomit)
identical(x=unique(metadata.v2[metadata.v2$infant.vomit.in.visit3=="yes",][["mid"]]), y=unique(visit3.vomit[visit3.vomit$visit3.vomit>0,"mid"]))

metadata.v2$infant.vomit.in.visit4 <- "no"
metadata.v2[c(which(metadata.v2$mid %in% (visit4.vomit[visit4.vomit$visit4.vomit>0,][["mid"]]))), "infant.vomit.in.visit4"] <- "yes"
length(unique(metadata.v2[metadata.v2$infant.vomit.in.visit4=="yes",][["mid"]]))
table(visit4.vomit$visit4.vomit)
identical(x=unique(metadata.v2[metadata.v2$infant.vomit.in.visit4=="yes",][["mid"]]), y=unique(visit4.vomit[visit4.vomit$visit4.vomit>0,"mid"]))

#a column for vomit in visit 3 and/or visit 4
metadata.v2$infant.vomit.in.visit3.or.visit4 <- "no"
metadata.v2[metadata.v2$infant.vomit.in.visit3=="yes", "infant.vomit.in.visit3.or.visit4"] <- "yes"
metadata.v2[metadata.v2$infant.vomit.in.visit4=="yes", "infant.vomit.in.visit3.or.visit4"] <- "yes"
length(unique(metadata.v2[metadata.v2$infant.vomit.in.visit3.or.visit4=="yes",][["mid"]]))


#change the order of columns, and omit the column that was a place holder for row names
metadata.v3 <- metadata.v2[, c(4:5,1:2,6:47)]

write.csv(metadata.v3, "infant.morbidities.medications.csv")

#STEP 12: Gather information on the total reported cases of morbidity, duration of
#morbidities, if they were diagnosed, and if they were treated.

#gather the total number of diarrhea cases.
diarrhea.cases <- redcap.wide[,c(1,grep(pattern="f[234]04_diarrhoea_q5_1$|f[234]04_diarrhoea[1234]_q6_1_[1234]$", colnames(redcap.wide)))]
diarrhea.cases <- pivot_longer(data=diarrhea.cases, cols=c(2:16))
diarrhea.cases <- diarrhea.cases[!is.na(diarrhea.cases$value),]
diarrhea.cases$value <- as.numeric(diarrhea.cases$value)
table(diarrhea.cases$value)
#42 reports of diarrhea.
length(unique(diarrhea.cases[diarrhea.cases$value=="1", ][["mid"]]))
#35 infants
#looking at total diarrhea cases per infant
diarrhea.cases.per.infant <-data.frame(table(diarrhea.cases[diarrhea.cases$value=="1",][["mid"]]))
summary(diarrhea.cases.per.infant$Freq)
#min 1, median 1, max 2

#gather the durations of diarrhea for all three visits. 
diarrhea.duration <- redcap.wide[,c(1,grep(pattern="f[234]04_duratn_q5_1_1|f[234]04_duratn_q6_1_", colnames(redcap.wide)))]
diarrhea.duration <- pivot_longer(data=diarrhea.duration, cols=c(2:16))
diarrhea.duration <- diarrhea.duration[!is.na(diarrhea.duration$value),]
diarrhea.duration$value <- as.numeric(diarrhea.duration$value)
dim(diarrhea.duration)
#42 answers for diarrhea duration
summary(diarrhea.duration$value)
#the median duration was 5 days, and the mean was 6.8. 
#the max was 25, the min was 1 day, the first quartile was
#2 days, third quartile was 10.8 days. 
jpeg(filename="histogram.diarrhea.duration.jpeg", width=6, height=6, res=600, units="in")
hist(diarrhea.duration$value, breaks=c(1:25), axes=F, main=NULL, xlab=NULL)
title(main=NULL, xlab="diarrhea duration (days)")
axis(side=1, at=c(0,2,5,12,14,21,25))
axis(side=2)
dev.off()

#gather the numbers for how many diarrhea cases were diagnosed
diarrhea.diagnosis <- redcap.wide[,c(1,grep(pattern="f[234]04_diag_q5_1_2|f[234]04_diag_q6_1_", colnames(redcap.wide)))]
diarrhea.diagnosis <- pivot_longer(data=diarrhea.diagnosis, cols=c(2:16))
diarrhea.diagnosis <- diarrhea.diagnosis[!is.na(diarrhea.diagnosis$value),]
diarrhea.diagnosis$value <- as.numeric(diarrhea.diagnosis$value)
dim(diarrhea.diagnosis)
#42 answers for diarrhea diagnosis
length(unique(diarrhea.diagnosis$mid))
#and these answers came from 35 infants. 
table(diarrhea.diagnosis$value)
#six of the 42 were diagnosed, and 36 were not. 

#gather the numbers for how many diarrhea cases were treated
diarrhea.treatment <- redcap.wide[,c(1,grep(pattern="f[234]04_treatmt_q5_1_3|f[234]04_treatmt_q6_1_", colnames(redcap.wide)))]
diarrhea.treatment <- pivot_longer(data=diarrhea.treatment, cols=c(2:16))
diarrhea.treatment <- diarrhea.treatment[!is.na(diarrhea.treatment$value),]
diarrhea.treatment$value <- as.numeric(diarrhea.treatment$value)
dim(diarrhea.treatment)
#for only 6 diarrhea cases, there were answers regarding treatment of the diarrhea
length(unique(diarrhea.treatment$mid))
#and these answers came from 6 infants. 
table(diarrhea.treatment$value)
#all six answers were "no". 

#gather the total number of fever cases.
fever.cases <- redcap.wide[,c(1,grep(pattern="f[234]04_fever_q5_5$|f[234]04_fever[1234]_q6_5_[1234]$", colnames(redcap.wide)))]
fever.cases <- pivot_longer(data=fever.cases, cols=c(2:16))
fever.cases <- fever.cases[!is.na(fever.cases$value),]
fever.cases$value <- as.numeric(fever.cases$value)
table(fever.cases$value)
#105 reports of fever. 
length(unique(fever.cases[fever.cases$value=="1",][["mid"]]))
#and these cases came from 65 infants. 
#looking at total fever cases per infant
fever.cases.per.infant <-data.frame(table(fever.cases[fever.cases$value=="1",][["mid"]]))
summary(fever.cases.per.infant$Freq)
#min 1, median 1, max 4

#gather the durations of fever for all three visits. 
fever.duration <- redcap.wide[,c(1,grep(pattern="f[234]04_duratn_q5_5_1|f[234]04_duratn_q6_5_", colnames(redcap.wide)))]
fever.duration <- pivot_longer(data=fever.duration, cols=c(2:16))
fever.duration <- fever.duration[!is.na(fever.duration$value),]
fever.duration$value <- as.numeric(fever.duration$value)
dim(fever.duration)
#105 answers for fever duration
summary(fever.duration$value)
#the median duration was 1 day, and the mean was 1.8. 
#the max was 8.5, the min was 0 day, the first quartile was
#1 day, third quartile was 2 days. 
jpeg(filename="histogram.fever.duration.jpeg", width=6, height=6, res=600, units="in")
hist(fever.duration$value, breaks=c(0:9), axes=F, main=NULL, xlab=NULL)
title(main=NULL, xlab="fever duration (days)")
axis(side=1, at=c(0:9))
axis(side=2)
dev.off()

#gather the numbers for how many fever cases were diagnosed
fever.diagnosis <- redcap.wide[,c(1,grep(pattern="f[234]04_diag_q5_5_2|f[234]04_diag_q6_5_", colnames(redcap.wide)))]
fever.diagnosis <- pivot_longer(data=fever.diagnosis, cols=c(2:16))
fever.diagnosis <- fever.diagnosis[!is.na(fever.diagnosis$value),]
fever.diagnosis$value <- as.numeric(fever.diagnosis$value)
dim(fever.diagnosis)
#104 of the fever cases reported if they were diagnosed or not. 
length(unique(fever.diagnosis$mid))
#answers came from 65 infants. 
table(fever.diagnosis$value)
#22 of the 104 cases were diagnosed, and 82 were not. 

#gather the numbers for how many fever cases were treated
fever.treatment <- redcap.wide[,c(1,grep(pattern="f[234]04_treatmt_q5_5_3|f[234]04_treatmt_q6_5_", colnames(redcap.wide)))]
fever.treatment <- pivot_longer(data=fever.treatment, cols=c(2:16))
fever.treatment <- fever.treatment[!is.na(fever.treatment$value),]
fever.treatment$value <- as.numeric(fever.treatment$value)
dim(fever.treatment)
#for only 22 fever cases, there were answers regarding treatment of the fever
length(unique(fever.treatment$mid))
#and these answers came from 19 infants. 
table(fever.treatment$value)
#eight answered "yes", and 14 answered "no".
#look at the treatment types
redcap.wide[c(which(redcap.wide$mid %in% c(fever.treatment[fever.treatment$value=="1",][["mid"]]))),c(grep(pattern="^mid|f[234]04_sptreatmt_q5_5_4|f[234]04_sptreatmt_q6_5", colnames(redcap.wide)))] %>% print()
#antibiotics, acetaminophen
redcap.wide[redcap.wide$mid=="MD070Y", "f304_sptreatmt_q6_4_1d"]
#antibiotic and painkiller

#gather the total number of vomit cases.
vomit.cases <- redcap.wide[,c(1,grep(pattern="f[234]04_vomitn_q5_2$|f[234]04_vomitn[1234]_q6_2_[1234]$", colnames(redcap.wide)))]
vomit.cases <- pivot_longer(data=vomit.cases, cols=c(2:16))
vomit.cases <- vomit.cases[!is.na(vomit.cases$value),]
vomit.cases$value <- as.numeric(vomit.cases$value)
table(vomit.cases$value)
#44 reports of vomit. 
length(unique(vomit.cases[vomit.cases$value=="1",][["mid"]]))
#and these cases came from 28 infants. 
vomit.cases.per.infant <-data.frame(table(vomit.cases[vomit.cases$value=="1",][["mid"]]))
summary(vomit.cases.per.infant$Freq)
#min 1, median 1, max 5. 

#gather the durations of vomit for all three visits. 
vomit.duration <- redcap.wide[,c(1,grep(pattern="f[234]04_duratn_q5_2_1|f[234]04_duratn_q6_2_", colnames(redcap.wide)))]
vomit.duration <- pivot_longer(data=vomit.duration, cols=c(2:16))
vomit.duration <- vomit.duration[!is.na(vomit.duration$value),]
vomit.duration$value <- as.numeric(vomit.duration$value)
dim(vomit.duration)
#44 reports for duration
summary(vomit.duration$value)
#the median duration was 1 day, and the mean was 1.023 day. 
#the max was 2 days, the min was 0.5 day, the first quartile was
#1 day, third quartile was 1 day. 
jpeg(filename="histogram.vomit.duration.jpeg", width=6, height=6, res=600, units="in")
hist(vomit.duration$value, breaks=c(0,0.5,1,1.5,2), axes=F, main=NULL, xlab=NULL)
title(main=NULL, xlab="vomit duration (days)")
axis(side=1, at=c(0,0.5,1,1.5,2))
axis(side=2)
dev.off()

#gather the numbers for how many vomit cases were diagnosed
vomit.diagnosis <- redcap.wide[,c(1,grep(pattern="f[234]04_diag_q5_2_2|f[234]04_diag_q6_2_", colnames(redcap.wide)))]
vomit.diagnosis <- pivot_longer(data=vomit.diagnosis, cols=c(2:16))
vomit.diagnosis <- vomit.diagnosis[!is.na(vomit.diagnosis$value),]
vomit.diagnosis$value <- as.numeric(vomit.diagnosis$value)
dim(vomit.diagnosis)
#44 answers for diagnosis
length(unique(vomit.diagnosis$mid))
#and these answers came from 28 infants. 
table(vomit.diagnosis$value)
#4 were diagnosed, and 40 were not. 

#gather the numbers for how many vomit cases were treated
vomit.treatment <- redcap.wide[,c(1,grep(pattern="f[234]04_treatmt_q5_2_3|f[234]04_treatmt_q6_2_", colnames(redcap.wide)))]
vomit.treatment <- pivot_longer(data=vomit.treatment, cols=c(2:16))
vomit.treatment <- vomit.treatment[!is.na(vomit.treatment$value),]
vomit.treatment$value <- as.numeric(vomit.treatment$value)
dim(vomit.treatment)
#for only 4 vomit cases, there were answers regarding treatment
length(unique(vomit.treatment$mid))
#and these answers came from 4 infants. 
table(vomit.treatment$value)
#two answers were "no", two were "yes".
#look at the treatment types
redcap.wide[c(which(redcap.wide$mid %in% c(vomit.treatment[vomit.treatment$value=="1",][["mid"]]))),c(grep(pattern="^mid|f[234]04_sptreatmt_q5_2_4|f[234]04_sptreatmt_q6_2", colnames(redcap.wide)))] %>% print()
#antiemetic and painkiller. 

#STEP 13: make another color-coded list of the morbidity by infant within each visit,
#but this time include a Venn diagram of the sick infants within that visit that
#to show how morbidities sometimes co-occurred in the same visit for an infant. 
list.visit2 <- list(diarrhea = medication.morbidity.table[medication.morbidity.table$visit=="2" & medication.morbidity.table$infant.diarrhea=="yes","mid"], fever = medication.morbidity.table[medication.morbidity.table$visit=="2" & medication.morbidity.table$infant.fever=="yes","mid"], vomit=medication.morbidity.table[medication.morbidity.table$visit=="2" & medication.morbidity.table$infant.vomit=="yes","mid"])
venndiag.visit2.morbidities <- ggvenn(data=list.visit2, columns= c("diarrhea", "fever", "vomit"), fill_color=c("#F8766D", "#00BA38", "#619CFF"), set_name_size=c(0,0,0), text_size=4, show_percentage=F) + labs(title="visit 2")+theme(plot.margin=unit(c(-0.5,-1,-0.5,-1), "cm"), plot.title=element_text(vjust=-10))

list.visit3 <- list(diarrhea = medication.morbidity.table[medication.morbidity.table$visit=="3" & medication.morbidity.table$infant.diarrhea=="yes","mid"], fever = medication.morbidity.table[medication.morbidity.table$visit=="3" & medication.morbidity.table$infant.fever=="yes","mid"], vomit=medication.morbidity.table[medication.morbidity.table$visit=="3" & medication.morbidity.table$infant.vomit=="yes","mid"])
venndiag.visit3.morbidities <- ggvenn(data=list.visit3, columns= c("diarrhea", "fever", "vomit"), fill_color=c("#F8766D", "#00BA38", "#619CFF"), set_name_size=0, text_size=4, show_percentage=F) + labs(title="visit 3")+theme(plot.margin=unit(c(-0.5,-1, -0.5,-1), "cm"), plot.title=element_text(vjust=-10))

list.visit4 <- list(diarrhea = medication.morbidity.table[medication.morbidity.table$visit=="4" & medication.morbidity.table$infant.diarrhea=="yes","mid"], fever = medication.morbidity.table[medication.morbidity.table$visit=="4" & medication.morbidity.table$infant.fever=="yes","mid"], vomit=medication.morbidity.table[medication.morbidity.table$visit=="4" & medication.morbidity.table$infant.vomit=="yes","mid"])
venndiag.visit4.morbidities <- ggvenn(data=list.visit4, columns= c("diarrhea", "fever", "vomit"), fill_color=c("#F8766D", "#00BA38", "#619CFF"), set_name_size=0, text_size=4, show_percentage=F) + labs(title="visit 4")+theme(plot.margin=unit(c(-0.5,-1,-0.5,-1), "cm"), plot.title=element_text(vjust=-10))

venndiag.visits234.morbidities <- grid.arrange(venndiag.visit2.morbidities, venndiag.visit3.morbidities, venndiag.visit4.morbidities, ncol=1, nrow=3)

colored.morbidity.list <-ggplot(df.illness.medx.v6, aes(x=visit.morbidity.rate, y=fct_reorder(cid,total.illness), fill=morbidity, alpha=value)) + geom_tile()+scale_alpha(range=c(0,1))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(text=element_text(size=8)) +geom_text(alpha=1,data=df.illness.medx.v6,aes(label=drug), color="black", size=1.5)+facet_wrap(~visit, scales="free_x")+ylab("infant gender, mode of delivery")+xlab("")+guides(alpha="none")+scale_y_discrete(breaks=c(df.illness.medx.v6$cid),labels=c(df.illness.medx.v6$infant.gender.mod))

plot.to.print <- grid.arrange(colored.morbidity.list, venndiag.visits234.morbidities, widths=c(1, 0.48))
ggsave("morbidity.medication.per.infant.noCID.withVennDiagrams.jpeg", plot=plot.to.print, dpi=600, height=10, width=10, units="in")
