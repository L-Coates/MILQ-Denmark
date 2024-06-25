#This R code is designed to gather infant vaccination information reported by moms
#and to relate it to fever events. 

library(tidyverse)

#step 1: read in the metadata file with morbidity information for the samples
#included in analysis. 
metadata <- read.csv("../morbidity_medication_prevalences/infant.morbidities.medications.csv")
metadata <- metadata[,-1]

#step 2: read in the redcap data and correct the two maternal IDs
redcap <- read.csv("../../../../REDCap_downloads/Denmark/MILQMAINSTUDYDENMARK_DATA_2022-10-11_2251.csv")
redcap[redcap$mid=="MD346m","mid"] <- "MD346M"
redcap[redcap$mid=="MD378x","mid"] <- "MD378X"

#step 3: pull out the redcap data for the 109 infants included in the analyses. 
redcap.v2 <- redcap[c(which(redcap$mid %in% metadata$mid)),]

#step 4: See if any moms reported vaccinations for their infants
table(redcap.v2$f204_immz_q11)
#10 infants had immunizations sometime before 3.5 months of life. 
table(redcap.v2$f304_immz_q11)
#107 infants had immunizations sometime between their visit 2 and visit 3.
table(redcap.v2$f404_immz_q11)
#82 infants had immunizations sometime between their visit 3 and visit 4.

#step 5: since so many infants had vaccinations (many more than the number of infants
#with a fever in the same visit as the vaccination) there won't be a simple 
#relationship to draw between fever and vaccinations. And no question was
#asked in the questionnaires for the mothers to state the diagnosis of the fever,
#so can't determine if the fever was due to the vaccine.
#So, will create a bar graph of the number of infants with vaccination
#reported at each visit. And include a percentage of those infants that had fever in the same visit. 
vaccines <- redcap.v2[, c(grep(pattern="f[234]04_immz_q11|^mid", colnames(redcap.v2)))]

vaccines.long <- pivot_longer(data=vaccines, cols=2:4, names_to = "visit", values_to = "vaccinated")
vaccines.long$visit <- gsub(pattern="f204_immz_q11", replacement="2", vaccines.long$visit)
vaccines.long$visit <- gsub(pattern="f304_immz_q11", replacement="3", vaccines.long$visit)
vaccines.long$visit <- gsub(pattern="f404_immz_q11", replacement="4", vaccines.long$visit)

vaccines.v2 <- merge(x=vaccines.long,y=metadata, by=c("mid", "visit"))

vaccines.v2 <- vaccines.v2[!is.na(vaccines.v2$vaccinated),]

#gather percentages of infants that had +/- fever and +/- vaccination in the same visit. 
table(vaccines.v2[vaccines.v2$visit=="2"&vaccines.v2$vaccinated=="1", "infant.fever"]) 
#visit 2, infants that were vaccinated, 40 % had fever, 60 % did not. 
table(vaccines.v2[vaccines.v2$visit=="2"&vaccines.v2$vaccinated=="0", "infant.fever"]) 
#visit 2, infants that weren't vaccinated, 4 % had fever, 96 % did not have fever. 

#gather percentages of infants that had +/- fever and +/- vaccination in the same visit. 
table(vaccines.v2[vaccines.v2$visit=="3"&vaccines.v2$vaccinated=="1", "infant.fever"]) 
#visit 3, infants that were vaccinated, 38 % had fever, 62 % did not. 
table(vaccines.v2[vaccines.v2$visit=="3"&vaccines.v2$vaccinated=="0", "infant.fever"]) 
#visit 3, infants that weren't vaccinated, 100 % did not have fever. 

#gather percentages of infants that had +/- fever and +/- vaccination in the same visit. 
table(vaccines.v2[vaccines.v2$visit=="4"&vaccines.v2$vaccinated=="1", "infant.fever"]) 
#visit 4, infants that were vaccinated, 39 % had fever, 61 % did not. 
table(vaccines.v2[vaccines.v2$visit=="4"&vaccines.v2$vaccinated=="0", "infant.fever"]) 
#visit 4, infants that weren't vaccinated, 37 % had fever, 63 % did not have fever. 

#add these percentages of vaccinated to the dataframe
table(vaccines.long[vaccines.long$visit=="2", "vaccinated"]) #10 vaccinated, 99 not, in visit 2. 
table(vaccines.long[vaccines.long$visit=="3", "vaccinated"]) #107 vaccinated, 2 not, in visit 3. 
table(vaccines.long[vaccines.long$visit=="4", "vaccinated"]) #82 vaccinated, 27 not, in visit 4. 

fever.rates.among.vaccinated <- data.frame(visit=c("2", "2", "3", "3", "4", "4"), total.vaccinated=c(10,99, 107, 2, 82, 27), vaccinated=c("yes", "no"), fever.rate=c("40 %", "4 %", "38 %","0 %", "39 %","37 %"), row.names = NULL)

ggplot(aes(x=visit,y=total.vaccinated, fill=as.factor(vaccinated), label=fever.rate), data=fever.rates.among.vaccinated)+geom_bar(stat="identity")+scale_fill_manual(values=c("black", "brown"))+labs(fill="vaccinated")+ylab("number of infants")+scale_y_continuous(breaks=c(0, 30, 60, 90, 109), limits=c(0,109))+geom_text(size = 3, position = position_stack(vjust = 0.5), color="white")+theme_dark()
ggsave(filename="total_vaccinated_infants_and_fever_rates_by_visit.jpeg",plot=last_plot(), dpi=600, width=4, height=3, units="in")
