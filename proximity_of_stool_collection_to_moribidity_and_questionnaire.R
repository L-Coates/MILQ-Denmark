#Part 1 of this script has the code used to determine how many days passed between the morbidity
#questionnaire and infant stool collection and to determine if stool collection
#happened before or after morbidity questionnaires. 
#Part 2 of this script has the code used to determine how many morbidity cases 
#(for diarrhea, fever or vomiting) definitively occurred within a week of stool sample collection
#and morbidity questionnaire collection. 
#Part 3 of this script has the code used to determine how many morbidity cases
#definitively occurred before stool collection. 

#Part 1...

#read in the Denmark REDCap dataset and look across the data to determine
#the time that passed between stool collection and morbidity questionnaire. 
redcap <- read.csv("../../../../REDCap_downloads/Denmark/MILQMAINSTUDYDENMARK_DATA_2022-10-11_2251.csv")
str(redcap[,c(grep(pattern="^mid|f[234]11_dosc_q5|f[234]04_doi_q3", colnames(redcap)))])

#we already know that MD378X and MD346M have lowercase letters in the redcap dataset
#so changing from lowercase to uppercase
redcap[redcap$mid=="MD378x","mid"] <- "MD378X"
redcap[redcap$mid=="MD346m","mid"] <- "MD346M"

#the dates are currently stored as character strings so changing to date. 
redcap$f204_doi_q3 <- as.Date(redcap$f204_doi_q3, format="%Y-%m-%d")
redcap$f304_doi_q3 <- as.Date(redcap$f304_doi_q3, format="%Y-%m-%d")
redcap$f404_doi_q3 <- as.Date(redcap$f404_doi_q3, format="%Y-%m-%d")

redcap$f211_dosc_q5 <- as.Date(redcap$f211_dosc_q5, format="%Y-%m-%d")
redcap$f311_dosc_q5 <- as.Date(redcap$f311_dosc_q5, format="%Y-%m-%d")
redcap$f411_dosc_q5 <- as.Date(redcap$f411_dosc_q5, format="%Y-%m-%d")

redcap$days_stool_collection_to_morbidity_questionnaire.visit2 <- as.numeric(difftime(time1=redcap$f211_dosc_q5, time2=redcap$f204_doi_q3, units="days"))
redcap$days_stool_collection_to_morbidity_questionnaire.visit3 <- as.numeric(difftime(time1=redcap$f311_dosc_q5, time2=redcap$f304_doi_q3, units="days"))
redcap$days_stool_collection_to_morbidity_questionnaire.visit4 <- as.numeric(difftime(time1=redcap$f411_dosc_q5, time2=redcap$f404_doi_q3, units="days"))

summary(redcap$days_stool_collection_to_morbidity_questionnaire.visit2)
unique(redcap$days_stool_collection_to_morbidity_questionnaire.visit2)
table(redcap$days_stool_collection_to_morbidity_questionnaire.visit2)
summary(redcap$days_stool_collection_to_morbidity_questionnaire.visit3)
unique(redcap$days_stool_collection_to_morbidity_questionnaire.visit3)
table(redcap$days_stool_collection_to_morbidity_questionnaire.visit3)
summary(redcap$days_stool_collection_to_morbidity_questionnaire.visit4)
unique(redcap$days_stool_collection_to_morbidity_questionnaire.visit4)
table(redcap$days_stool_collection_to_morbidity_questionnaire.visit4)

#selecting just the samples that were used in the infant stool microbe associations
#with morbidity (i.e. the subset of 327 stool samples)

metadata <- read.csv("../metadata-infants-complete-stool-set-after-6886-read-cutoff-withoutInfantCD144Y.csv")
metadata <- metadata[,-1]
metadata$mid <- gsub(metadata$extraction.id, pattern="[.][234][xb]*", replacement="")
metadata$mid <- gsub(metadata$mid, pattern="^C", replacement="M")
unique(metadata$mid)
redcap.v2 <- redcap[c(which(redcap$mid %in% metadata$mid)),]
length(unique(redcap.v2$mid)) #109

#most of the stool samples were collected before the morbidity questionnaire. 
#Looking at the dates for the stool samples and morbidity questionnaires taken
#more than a day apart. 

summary(redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit2)
unique(redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit2) 
table(redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit2) 

summary(redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit3)
unique(redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit3)
table(redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit3) 

summary(redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit4)
unique(redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit4)
table(redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit4) 
#almost all stool samples were collected before the morbidity questionnaire
#was collected. 
#But there are some samples from visit 2 and 4 that have either the stool collection date or the morbidity
#questionnaire date wrong because too many days (e.g. >75 days) passed between
#the stool collection and the morbidity questionnaire to be considered from the same
#visit or study. 

#looking at the interview dates, sample collection dates, and birth dates to help guide date correction. 
redcap.v2[c(abs(redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit2)>75 & !is.na(redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit2)),c(grep(pattern="^mid|f211_dosc_q5|f211_doi_q3|f204_doi_q3", colnames(redcap.v2)))]
redcap.v2[c(abs(redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit3)>75 & !is.na(redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit3)),c(grep(pattern="^mid|f311_dosc_q5|f311_doi_q3|f304_doi_q3", colnames(redcap.v2)))]
redcap.v2[c(abs(redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit4)>75 & !is.na(redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit4)),c(grep(pattern="^mid|f411_dosc_q5|f411_doi_q3|f404_doi_q3", colnames(redcap.v2)))]
redcap.v2[redcap.v2$mid=="MD124Q", "f101_dob_q5"] #dob: 2018-08-03
redcap.v2[redcap.v2$mid=="MD199W", "f101_dob_q5"] #dob: 2018-11-06
redcap.v2[redcap.v2$mid=="MD182C", "f101_dob_q5"] #dob: 2018-10-16
#for visit 2, for both interview dates, the year is 2018, whereas the year for sample collection
#date is 2918 or 2019 so it looks like the date of sample collection should be corrected to
#2018. 
#For visit 3, there were no suspected incorrect dates for stool collection or interview date.
#For visit 4, the stool collection date (specifically, the year) is incorrect, 
#because it is earlier than the child's date of birth. It should be changed to the
#year 2019. 

redcap.v2[c(redcap.v2$mid=="MD124Q"&redcap.v2$redcap_event_name=="10__349_months_arm_1"),"f211_dosc_q5"] <- "2018-09-22"
redcap.v2[c(redcap.v2$mid=="MD199W"&redcap.v2$redcap_event_name=="10__349_months_arm_1"),"f211_dosc_q5"] <- "2018-12-15"
redcap.v2[c(redcap.v2$mid=="MD182C"&redcap.v2$redcap_event_name=="60__849_months_arm_1"),"f411_dosc_q5"] <- "2019-05-21"

#recalculate the days between stool sample collection and morbidity questionnaire
#for visits 2 and 4. 
redcap.v2$f211_dosc_q5 <- as.Date(redcap.v2$f211_dosc_q5, format="%Y-%m-%d")
redcap.v2$f411_dosc_q5 <- as.Date(redcap.v2$f411_dosc_q5, format="%Y-%m-%d")


redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit2 <- as.numeric(difftime(time1=redcap.v2$f211_dosc_q5, time2=redcap.v2$f204_doi_q3, units="days"))
redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit4 <- as.numeric(difftime(time1=redcap.v2$f411_dosc_q5, time2=redcap.v2$f404_doi_q3, units="days"))

summary(redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit2)
unique(redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit2)
table(redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit2)
summary(redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit4)
unique(redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit4)
table(redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit4)

#make a histogram of the days between stool collection and morbidity questionnaire.
jpeg(filename="days.from.morbidity.interview.to.stool.collection.jpeg", width=12, height=12, units="cm", res=600)
hist(x=(c(redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit2,redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit3, redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit4)), xlim = c(-35,30), breaks=35, xlab="days from morbidity interview to stool collection", main="", ylab="number of stool samples", xaxt='n', yaxt='n')
axis(side=1, at=seq(-35,30, 5), labels=seq(-35,30,5))
axis(side=2, at=seq(0,180,20), labels=seq(0,180,20))
abline(v=c(-6, 0), col="red",lty="dashed")
dev.off()

#determine the median number of days passed between stool collection and morbidity
#questionnaire for all stool samples (regardless of visit) and for stool samples
#within each visit. 
median(c(redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit2,redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit3, redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit4), na.rm=T)
#-1 day for all stool samples (i.e. stool collected the day before the questionnaire was completed)

median(redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit2, na.rm=T)
#-1 day for stool samples from visit 2. 

median(redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit3, na.rm=T)
#-1 day for stool samples from visit 3

median(redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit4, na.rm=T)
#-1 day for stool samples from visit 4. 

#how many stool samples were collected within a week of the questionnaire?
days.between.stool.and.questionnaire <- as.data.frame(table(c(redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit2,redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit3, redcap.v2$days_stool_collection_to_morbidity_questionnaire.visit4)))
days.between.stool.and.questionnaire$Var1 <- as.numeric(as.character(days.between.stool.and.questionnaire$Var1))
days.between.stool.and.questionnaire$Freq <- as.numeric(days.between.stool.and.questionnaire$Freq)
sum(c(days.between.stool.and.questionnaire[days.between.stool.and.questionnaire$Var1>=-6 & days.between.stool.and.questionnaire$Var1<=6,"Freq"]))
#304 stool samples were collected within the week around the interview (i.e, week prior or week following interview). 
sum(c(days.between.stool.and.questionnaire[days.between.stool.and.questionnaire$Var1>=-6 & days.between.stool.and.questionnaire$Var1<=0,"Freq"]))
#298 stool samples were collected within the week prior to the interview. 

#Part 2 ...

#Looking within each visit at how many infants had a morbidity of interest (diarrhea,
#fever, or vomiting) within the week leading up to the morbidity questionnaire,
#including the day of the interview, and
#determining if their stool sample was also collected within that week. 

#visit 2
table(redcap.v2$f204_diarrhoea_q5_1) #two infants had recent diarrhea in visit 2, out of 8 total infants that had diarrhea in visit 2.
redcap.v2[redcap.v2$f204_diarrhoea_q5_1==1 & !is.na(redcap.v2$f204_diarrhoea_q5_1),c("f204_ill_q5_9", "f204_vomitn_q5_2", "f204_cough_q5_3", "f204_rapdbreath_q5_4", "f204_fever_q5_5", "f204_cold_q5_6", "f204_earache_q5_7", "f204_other_q5_8")]
#only one that had diarrhea on the day of interview
redcap.v2[redcap.v2$f204_diarrhoea_q5_1==1 & !is.na(redcap.v2$f204_diarrhoea_q5_1),c("f204_ill_q5_9","days_stool_collection_to_morbidity_questionnaire.visit2")]
#stool was collected 1 day prior to morbidity questionnaire for the infant that was
#sick on day of interview, and stool was collected 2 days prior to interview for 
#the infant that was not sick on day of interview.

#vomiting in the past week
table(redcap.v2$f204_vomitn_q5_2) #two infants had recent vomiting in visit 2, out of 6 total infants that had vomit in visit 2. 
redcap.v2[redcap.v2$f204_vomitn_q5_2==1 & !is.na(redcap.v2$f204_vomitn_q5_2),c("f204_ill_q5_9", "f204_diarrhoea_q5_1", "f204_cough_q5_3", "f204_rapdbreath_q5_4", "f204_fever_q5_5", "f204_cold_q5_6", "f204_earache_q5_7", "f204_other_q5_8")]
#only one that had an illness on day of interview and it might not have been diarrhea since this infant
#has some "other" illness in the week prior to the interview. 
redcap.v2[redcap.v2$f204_vomitn_q5_2==1 & !is.na(redcap.v2$f204_vomitn_q5_2),c("f204_ill_q5_9","days_stool_collection_to_morbidity_questionnaire.visit2")]
#stool was collected 2 days prior to morbidity questionnaire for the infant that was
#sick on day of interview, and stool was collected 1 day prior to interview for 
#the infant that was not sick on day of interview.

#fever in the past week
table(redcap.v2$f204_fever_q5_5) #five infants had recent fever in visit 2, out of 8 total infants that had fever in visit 2. 
redcap.v2[redcap.v2$f204_fever_q5_5==1 & !is.na(redcap.v2$f204_fever_q5_5),c("f204_ill_q5_9", "f204_diarrhoea_q5_1", "f204_cough_q5_3", "f204_rapdbreath_q5_4", "f204_vomitn_q5_2", "f204_cold_q5_6", "f204_earache_q5_7", "f204_other_q5_8")]
#none had an illness on day of interview. 
redcap.v2[redcap.v2$f204_fever_q5_5==1 & !is.na(redcap.v2$f204_fever_q5_5),"days_stool_collection_to_morbidity_questionnaire.visit2"]
#their stool samples were all collected within 2 days, or less, prior to morbidity questionnaire. 

#visit 3
table(redcap.v2$f304_diarrhoea_q5_1) #three infants had recent diarrhea in visit 3, out of 15 total infants that had diarrhea in visit 3. 
redcap.v2[redcap.v2$f304_diarrhoea_q5_1==1 & !is.na(redcap.v2$f304_diarrhoea_q5_1),c("f304_ill_q5_9", "f304_vomitn_q5_2", "f304_cough_q5_3", "f304_rapdbreath_q5_4", "f304_fever_q5_5", "f304_cold_q5_6", "f304_earache_q5_7", "f304_other_q5_8")]
#only one that had an illness on day of interview and it likely was diarrhea. 
redcap.v2[redcap.v2$f304_diarrhoea_q5_1==1 & !is.na(redcap.v2$f304_diarrhoea_q5_1),c("f304_ill_q5_9","days_stool_collection_to_morbidity_questionnaire.visit3")]
#stool was collected 3 days prior to morbidity questionnaire for the infant that was
#sick on day of interview, and stool was collected 7 days and 2 days prior to interview for 
#the infants that were not sick on day of interview.
#For the infant stool collected 7 days prior to interview, was there also
#diarrhea reported outside of the week prior (i.e. "yes" to diarrhea since last visit)?
redcap.v2[redcap.v2$f304_diarrhoea_q5_1==1 & !is.na(redcap.v2$f304_diarrhoea_q5_1),c("f304_ill_q5_9","f304_diarrhoea_q6_1","days_stool_collection_to_morbidity_questionnaire.visit3")]
#for that particular infant stool, no diarrhea was reported since the last visit 
#(i.e. before the week of interview); so we can say that the stool was collected
#BEFORE the diarrhea. 

#vomiting in the past week
table(redcap.v2$f304_vomitn_q5_2) #four infants had recent vomiting in visit 3, out of 7 total infants that had vomit in visit 3. 
redcap.v2[redcap.v2$f304_vomitn_q5_2==1 & !is.na(redcap.v2$f304_vomitn_q5_2),c("f304_ill_q5_9", "f304_diarrhoea_q5_1", "f304_cough_q5_3", "f304_rapdbreath_q5_4", "f304_fever_q5_5", "f304_cold_q5_6", "f304_earache_q5_7", "f304_other_q5_8")]
#only one that had an illness on day of interview and it could have been a cold, cough
#or rapid breathing. 
redcap.v2[redcap.v2$f304_vomitn_q5_2==1 & !is.na(redcap.v2$f304_vomitn_q5_2),c("f304_ill_q5_9","days_stool_collection_to_morbidity_questionnaire.visit3")]
#stool was collected 2 days prior to morbidity questionnaire for the infant that was
#sick on day of interview, and stool was collected 1, 2, and 7 days prior to interview for 
#the infants that were not sick on day of interview.
#For the infant stool collected 7 days prior to interview, was there also
#vomit reported outside of the week prior (i.e. "yes" to vomit since last visit)?
redcap.v2[redcap.v2$f304_vomitn_q5_2==1 & !is.na(redcap.v2$f304_vomitn_q5_2),c("f304_ill_q5_9","f304_vomitn_q6_2","days_stool_collection_to_morbidity_questionnaire.visit3")]
#for that particular infant stool, no vomit was reported since the last visit 
#(i.e. before the week of interview); so we can't determine relative timing
#of stool collection to vomit.  

#fever in the past week
table(redcap.v2$f304_fever_q5_5) #nine infants had recent fever in visit 3, out of 41 total infants that had fever in visit 3. 
redcap.v2[redcap.v2$f304_fever_q5_5==1 & !is.na(redcap.v2$f304_fever_q5_5),c("f304_ill_q5_9", "f304_diarrhoea_q5_1", "f304_cough_q5_3", "f304_rapdbreath_q5_4", "f304_vomitn_q5_2", "f304_cold_q5_6", "f304_earache_q5_7", "f304_other_q5_8")]
#only one that had an illness on day of interview and it could have been a cold or cough. 
redcap.v2[redcap.v2$f304_fever_q5_5==1 & !is.na(redcap.v2$f304_fever_q5_5),c("f304_ill_q5_9","days_stool_collection_to_morbidity_questionnaire.visit3")]
#stool was collected 4 days prior to morbidity questionnaire for the infant that was
#sick on day of interview, and stool was collected 1, 2, 4, 5, 6, and 7 days prior to interview for 
#the infants that were not sick on day of interview.
#For the infant stools collected 7 days prior to interview, was there also
#fever reported outside of the week prior (i.e. "yes" to fever since last visit)?
redcap.v2[redcap.v2$f304_fever_q5_5==1 & !is.na(redcap.v2$f304_fever_q5_5),c("f304_ill_q5_9","f304_fever_q6_5","days_stool_collection_to_morbidity_questionnaire.visit3")]
#for those particular infant stools, no fever was reported since the last visit 
#(i.e. before the week of interview); so we can say that the stools were collected
#BEFORE the fever. 

#visit 4
table(redcap.v2$f404_diarrhoea_q5_1) #four infants had recent diarrhea in visit 4, out of 17 total infants that had diarrhea in visit 4. 
redcap.v2[redcap.v2$f404_diarrhoea_q5_1==1 & !is.na(redcap.v2$f404_diarrhoea_q5_1),c("f404_ill_q5_9", "f404_vomitn_q5_2", "f404_cough_q5_3", "f404_rapdbreath_q5_4", "f404_fever_q5_5", "f404_cold_q5_6", "f404_earache_q5_7", "f404_other_q5_8")]
#none had illness on the day of interview
redcap.v2[redcap.v2$f404_diarrhoea_q5_1==1 & !is.na(redcap.v2$f404_diarrhoea_q5_1),"days_stool_collection_to_morbidity_questionnaire.visit4"]
#one stool was collected 28 days after the morbidity questionnaire, while the other
#three stool samples were collected within 3 days prior to the questionnaire.

#vomiting in the past week
table(redcap.v2$f404_vomitn_q5_2) #three infants had recent vomiting in visit 4, out of 20 total infants that had vomit in visit 4. 
redcap.v2[redcap.v2$f404_vomitn_q5_2==1 & !is.na(redcap.v2$f404_vomitn_q5_2),c("f404_ill_q5_9", "f404_diarrhoea_q5_1", "f404_cough_q5_3", "f404_rapdbreath_q5_4", "f404_fever_q5_5", "f404_cold_q5_6", "f404_earache_q5_7", "f404_other_q5_8")]
#one had illness on the day of interview, and it could have been fever, cold or "other". 
redcap.v2[redcap.v2$f404_vomitn_q5_2==1 & !is.na(redcap.v2$f404_vomitn_q5_2),c("f404_ill_q5_9","days_stool_collection_to_morbidity_questionnaire.visit4")]
#stool was collected on the same day as the morbidity questionnaire for the infant that was
#sick on day of interview, and stool was collected 34 days prior to interview and 28
#day after interview for the infants that were not sick on day of interview.
#For the infant stool collected 34 days prior to interview, was there also
#vomit reported outside of the week prior (i.e. "yes" to vomit since last visit)?
redcap.v2[redcap.v2$f404_vomitn_q5_2==1 & !is.na(redcap.v2$f404_vomitn_q5_2),c("f404_ill_q5_9","f404_vomitn_q6_2","days_stool_collection_to_morbidity_questionnaire.visit4")]
#for that particular infant stool, yes vomit was reported since the last visit 
#(i.e. before the week of interview); so we don't know if the stool was collected
#BEFORE or AFTER the vomit. 

#fever in the past week
table(redcap.v2$f404_fever_q5_5) #fourteen infants had recent fever in visit 4, out of 42 total infants that had fever in visit 4. 
redcap.v2[redcap.v2$f404_fever_q5_5==1 & !is.na(redcap.v2$f404_fever_q5_5),c("f404_ill_q5_9", "f404_diarrhoea_q5_1", "f404_cough_q5_3", "f404_rapdbreath_q5_4", "f404_vomitn_q5_2", "f404_cold_q5_6", "f404_earache_q5_7", "f404_other_q5_8")]
#seven had illness on the day of interview, and for all of them they had some 
#other illness listed for that same week of the interview so it might not have been
#fever that they had on the day of the interview. 
redcap.v2[redcap.v2$f404_fever_q5_5==1 & !is.na(redcap.v2$f404_fever_q5_5),c("f404_ill_q5_9","days_stool_collection_to_morbidity_questionnaire.visit4")]
#stool was collected on 0, 1, 2, and 5 days prior to questionnaire for the infants
#that were sick on the day of the interview. For the others sick during the week of the 
#interview, but not during the interview, the stool was collected, 1, 2, 5, and 34
#days prior to the interview.
#For the infant stool collected 34 days prior to morbidity interview, was there also
#fever reported outside of the week prior (i.e. "yes" to fever since last visit)?
redcap.v2[redcap.v2$f404_fever_q5_5==1 & !is.na(redcap.v2$f404_fever_q5_5),c("f404_ill_q5_9","f404_fever_q6_5","days_stool_collection_to_morbidity_questionnaire.visit4")]
#for that particular infant stool, no fever was reported since the last visit 
#(i.e. before the week of interview); so we can say that the stool was collected
#BEFORE the fever. 

#Part 3...

#Looking within each visit, to determine how many infants with reported morbidity
#had a stool sample collected, definitively, after the morbidity. 

#visit 2

redcap.v2[redcap.v2$f204_diarrhoea_q5_1==0 & !is.na(redcap.v2$f204_diarrhoea_q5_1) & redcap.v2$f204_diarrhoea_q6_1==1 & !is.na(redcap.v2$f204_diarrhoea_q6_1),"days_stool_collection_to_morbidity_questionnaire.visit2"]
#five diarrhea-afflicted infants (out of 8 total) in visit 2 had stool collected 
#during the week prior to the interview, and one had stool collected one day after
#interview, but diarrhea was not during that week for these infants.
#For that stool that was collected a day after interview, was diarrhea reported
#for the subsequent visit (before the week of visit 3 morbidity interview)?
redcap.v2[redcap.v2$f204_diarrhoea_q5_1==0 & !is.na(redcap.v2$f204_diarrhoea_q5_1) & redcap.v2$f204_diarrhoea_q6_1==1 & !is.na(redcap.v2$f204_diarrhoea_q6_1),c("f304_diarrhoea_q6_1","days_stool_collection_to_morbidity_questionnaire.visit2")]
#no.

redcap.v2[redcap.v2$f204_vomitn_q5_2==0 & !is.na(redcap.v2$f204_vomitn_q5_2) & redcap.v2$f204_vomitn_q6_2==1 & !is.na(redcap.v2$f204_vomitn_q6_2),"days_stool_collection_to_morbidity_questionnaire.visit2"]
#two vomit-afflicted infants (out of 6 total) in visit 2 had stool collected 
#during the week prior to the interview but vomit was not during that week. Two
#other infants had stool collected more than a week before the morbidity interview so
#do not know relative timing of stool collection to vomiting. 

redcap.v2[redcap.v2$f204_fever_q5_5==0 & !is.na(redcap.v2$f204_fever_q5_5) & redcap.v2$f204_fever_q6_5==1 & !is.na(redcap.v2$f204_fever_q6_5),"days_stool_collection_to_morbidity_questionnaire.visit2"]
#three fever-afflicted infants (out of 8 total) in visit 2 had stool collected 
#during the week prior to the interview but fever was not during that week. 

#visit 3

redcap.v2[redcap.v2$f304_diarrhoea_q5_1==0 & !is.na(redcap.v2$f304_diarrhoea_q5_1) & redcap.v2$f304_diarrhoea_q6_1==1 & !is.na(redcap.v2$f304_diarrhoea_q6_1),"days_stool_collection_to_morbidity_questionnaire.visit3"]
#12 diarrhea-afflicted infants (out of 15 total) in visit 3 had stool collected 
#during the week prior to the interview but diarrhea was not during that week. 

redcap.v2[redcap.v2$f304_vomitn_q5_2==0 & !is.na(redcap.v2$f304_vomitn_q5_2) & redcap.v2$f304_vomitn_q6_2==1 & !is.na(redcap.v2$f304_vomitn_q6_2),"days_stool_collection_to_morbidity_questionnaire.visit3"]
#three vomit-afflicted infants (out of 7 total) in visit 3 had stool collected 
#during the week prior to the interview but vomit was not during that week. 

redcap.v2[redcap.v2$f304_fever_q5_5==0 & !is.na(redcap.v2$f304_fever_q5_5) & redcap.v2$f304_fever_q6_5==1 & !is.na(redcap.v2$f304_fever_q6_5),"days_stool_collection_to_morbidity_questionnaire.visit3"]
#30 fever-afflicted infants (out of 41 total) in visit 3 had stool collected 
#during the week prior to the interview or soon after the interview,
#but fever was not during that week. 
#For the infant with stool collected two days after visit 3 morbidity interview,
#was there fever reported in visit 4 (sometime between visit 3 interview and before
#week leading up to visit 4 interview)?
redcap.v2[redcap.v2$f304_fever_q5_5==0 & !is.na(redcap.v2$f304_fever_q5_5) & redcap.v2$f304_fever_q6_5==1 & !is.na(redcap.v2$f304_fever_q6_5),c("f404_fever_q6_5","days_stool_collection_to_morbidity_questionnaire.visit3")]
#no

#visit 4

redcap.v2[redcap.v2$f404_diarrhoea_q5_1==0 & !is.na(redcap.v2$f404_diarrhoea_q5_1) & redcap.v2$f404_diarrhoea_q6_1==1 & !is.na(redcap.v2$f404_diarrhoea_q6_1),"days_stool_collection_to_morbidity_questionnaire.visit4"]
#12 diarrhea-afflicted infants (out of 17 total) in visit 4 had stool collected 
#during the week prior to the interview but diarrhea was not during that week. 

redcap.v2[redcap.v2$f404_vomitn_q5_2==0 & !is.na(redcap.v2$f404_vomitn_q5_2) & redcap.v2$f404_vomitn_q6_2==1 & !is.na(redcap.v2$f404_vomitn_q6_2),"days_stool_collection_to_morbidity_questionnaire.visit4"]
#15 vomit-afflicted infants (out of 20 total) in visit 4 had stool collected 
#during the week prior to the interview but vomit was not during that week. One had
#had stool collected one day after interview; so we can still say that stool was 
#collected after vomiting. 

redcap.v2[redcap.v2$f404_fever_q5_5==0 & !is.na(redcap.v2$f404_fever_q5_5) & redcap.v2$f404_fever_q6_5==1 & !is.na(redcap.v2$f404_fever_q6_5),"days_stool_collection_to_morbidity_questionnaire.visit4"]
#24 fever-afflicted infants (out of 42 total) in visit 4 had stool collected 
#during the week prior to the interview but fever was not during that week. One 
#stool was collected 9 days after morbidity interview but we can still say that 
#stool was collected after the fever. 

