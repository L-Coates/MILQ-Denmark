#This script describes the steps taken to assess the metadata on the Danish
#infant stool collection and storage, including age at stool collection, and time
#between subsequent stool collections. 

#STEP 1: read in the MILQ Denmark metadata (downloaded from REDCap) and select just the 
#metadata related to the stool samples that we analyzed. 
redcap <- read.csv('../../../../REDCap_downloads/Denmark/MILQMAINSTUDYDENMARK_DATA_2022-10-11_2251.csv')

#we already know that there are two maternal IDs (MIDs) in the REDCap dataset that have 
#a lowercase check letter instead of the standard uppercase check letter. Correcting
#these MIDs now
#the closest MIDs in REDCap set are "MD36m" and "MD378x" so need to change those to capital letters
redcap[redcap$mid=="MD346m","mid"] <- "MD346M"
redcap[redcap$mid=="MD378x","mid"] <- "MD378X"

#Read in Kable lab data on the 330 Denmark infant stool samples that are
#being used for analysis. 
metadata <- read.csv("../metadata-infants-complete-stool-set-after-6886-read-cutoff.csv")
#first column is empty
metadata <- metadata[,-1]

#make a maternal ID column and a visit column facilitate matching with the redcap data
metadata$mid <- metadata$extraction.id
metadata$mid <- gsub(pattern="^C", replacement="M", metadata$mid)
metadata$mid <- gsub(pattern="[.][2-4][x|b]*", replacement="", metadata$mid)
metadata$visit <- metadata$extraction.id
metadata$visit <- gsub(pattern="CD[0-9]*[A-Z][.]", replacement="", metadata$visit)
metadata$visit <- gsub(pattern="x|b", replacement="", metadata$visit)

redcap.v2 <- redcap[c(which(redcap$mid %in% metadata$mid)),]

#create one row per sample
unique(redcap.v2$redcap_event_name)
redcap.visit1 <- redcap.v2[redcap.v2$redcap_event_name=="2__3_days_postpart_arm_1",c(1,grep(colnames(redcap.v2), pattern="f101"))]
redcap.visit2 <- redcap.v2[redcap.v2$redcap_event_name=="10__349_months_arm_1",c(1,grep(colnames(redcap.v2), pattern="f2|redcap_event_name"))]
redcap.visit3 <- redcap.v2[redcap.v2$redcap_event_name=="35__599_months_arm_1",c(1,grep(colnames(redcap.v2), pattern="f3|redcap_event_name"))]
redcap.visit4 <- redcap.v2[redcap.v2$redcap_event_name=="60__849_months_arm_1",c(1,grep(colnames(redcap.v2), pattern="f4|redcap_event_name"))]

redcap.v3 <- merge(redcap.visit1, redcap.visit2, by="mid", all=F)
redcap.v3 <- merge(redcap.v3, redcap.visit3, by="mid", all=F)
redcap.v3 <- merge(redcap.v3, redcap.visit4, by="mid", all=F)


#STEP 2: calculate the age at time of stool collection
#the column containing infant age is "f101_dob_q5"
#the column containing date of stool collection in visit 2 (1-3.49 months) is "f211_dosc_q5"
#the column containing date of stool collection in visit 3 (3.5-5.9 months) is "f311_dosc_q5"
#the column containing date of stool collection in visit 4 (6-8.49 months) is "f411_dosc_q5"
str(redcap.v3$f101_dob_q5)
str(redcap.v3$f211_dosc_q5)
str(redcap.v3$f311_dosc_q5)
str(redcap.v3$f411_dosc_q5)

#all date columns are currently class character so need to turn into class date
#but put in a new column just in case the format/date gets altered during conversion
redcap.v3$f101_dob_q5_v2 <- as.Date(redcap.v3$f101_dob_q5, format="%Y-%m-%d")
redcap.v3$f211_dosc_q5_v2 <- as.Date(redcap.v3$f211_dosc_q5, format="%Y-%m-%d")
redcap.v3$f311_dosc_q5_v2 <- as.Date(redcap.v3$f311_dosc_q5, format="%Y-%m-%d")
redcap.v3$f411_dosc_q5_v2 <- as.Date(redcap.v3$f411_dosc_q5, format="%Y-%m-%d")

#use the difftime function to calculate infant age in days
redcap.v3$visit2.stool.age <- difftime(time1=redcap.v3$f211_dosc_q5_v2, time2=redcap.v3$f101_dob_q5_v2, units="days")
redcap.v3$visit2.stool.age <- as.numeric(redcap.v3$visit2.stool.age)

redcap.v3$visit3.stool.age <- difftime(time1=redcap.v3$f311_dosc_q5_v2, time2=redcap.v3$f101_dob_q5_v2, units="days")
redcap.v3$visit3.stool.age <- as.numeric(redcap.v3$visit3.stool.age)

redcap.v3$visit4.stool.age <- difftime(time1=redcap.v3$f411_dosc_q5_v2, time2=redcap.v3$f101_dob_q5_v2, units="days")
redcap.v3$visit4.stool.age <- as.numeric(redcap.v3$visit4.stool.age)

#confirming that the infant ages fit in the required infant age for that visit
#in the MILQ study they used a conversion of 30.4375 days in one month
#so dividing the age in days by 30.4375 to get age in months. 
redcap.v3$visit2.stool.age.months <- redcap.v3$visit2.stool.age/30.4375
redcap.v3$visit3.stool.age.months <- redcap.v3$visit3.stool.age/30.4375
redcap.v3$visit4.stool.age.months <- redcap.v3$visit4.stool.age/30.4375

range(redcap.v3$visit2.stool.age.months)
range(redcap.v3$visit3.stool.age.months)
#there is a missing entry so looking to see which infant stool missed a collection date
redcap.v3[is.na(redcap.v3$visit3.stool.age.months), "mid"]
#MD144Y, so we will look to see if the stool was collected at all
redcap.v3[redcap.v3$mid=="MD144Y", "f311_inffeces_q4"]
#apparently the infant feces wasn't collected so can't be certain
#that the stool sample received that's marked CD144Y.3 is actually
#from infant CD144Y.
#so I am removing CD144Y from our subset of samples for analyses. 
redcap.v4 <- redcap.v3[redcap.v3$mid!="MD144Y",]
range(redcap.v4$visit3.stool.age.months)

range(redcap.v4$visit4.stool.age.months)

#look closely at the dates for the stool collection ages that are out of range
out.of.range <- redcap.v4[redcap.v4$visit2.stool.age.months<1|redcap.v4$visit2.stool.age.months>3.5|redcap.v4$visit3.stool.age.months<3.5|redcap.v4$visit3.stool.age.months>6|redcap.v4$visit4.stool.age.months<6|redcap.v4$visit4.stool.age.months>8.5,]
out.of.range <- out.of.range[,c(1,grep(pattern="f101_d|f211_d|f311_d|f411_d|visit[2-4]", colnames(out.of.range)))]

#in looking through the dates, it's clear for infant CD124Q that the year (2918) 
#put in for the stool collection in visit 2 was wrong, and should have been 2018
#which would have matched with the date for birth and questionnaire interview dates. 
#changing the visit 2 stool collection date for infant CD124Q from 2918-09-22 to 2018-09-22
redcap.v4[redcap.v4$mid=="MD124Q", "f211_dosc_q5_v2"]
redcap.v4[redcap.v4$mid=="MD124Q", "f211_dosc_q5_v2"] <- "2018-09-22"
redcap.v4[redcap.v4$mid=="MD124Q", "f211_dosc_q5_v2"]

#for infant belonging to MD182C, the date of stool collection for visit 4 was wrong
#and should have been 2019. 
redcap.v4[redcap.v4$mid=="MD182C", "f411_dosc_q5_v2"]
redcap.v4[redcap.v4$mid=="MD182C", "f411_dosc_q5_v2"] <- "2019-05-21"
redcap.v4[redcap.v4$mid=="MD182C", "f411_dosc_q5_v2"]

#for infant belonging to MD199W, the date of stool collection for visit 2 was wrong
#and should have been 2018. 
redcap.v4[redcap.v4$mid=="MD199W", "f211_dosc_q5_v2"]
redcap.v4[redcap.v4$mid=="MD199W", "f211_dosc_q5_v2"] <- "2018-12-15"
redcap.v4[redcap.v4$mid=="MD199W", "f211_dosc_q5_v2"]

#recalculate ages now
redcap.v4$f211_dosc_q5_v2 <- as.Date(redcap.v4$f211_dosc_q5_v2, format="%Y-%m-%d")
redcap.v4$f411_dosc_q5_v2 <- as.Date(redcap.v4$f411_dosc_q5_v2, format="%Y-%m-%d")

redcap.v4$visit2.stool.age <- difftime(time1=redcap.v4$f211_dosc_q5_v2, time2=redcap.v4$f101_dob_q5_v2, units="days")
redcap.v4$visit2.stool.age <- as.numeric(redcap.v4$visit2.stool.age)

redcap.v4$visit4.stool.age <- difftime(time1=redcap.v4$f411_dosc_q5_v2, time2=redcap.v4$f101_dob_q5_v2, units="days")
redcap.v4$visit4.stool.age <- as.numeric(redcap.v4$visit4.stool.age)

redcap.v4$visit2.stool.age.months <- redcap.v4$visit2.stool.age/30.4375
redcap.v4$visit4.stool.age.months <- redcap.v4$visit4.stool.age/30.4375


#STEP 3: change data frame into long format and generate plot of infant age at 
#time of stool collection with coloring for visit number
age.at.stool <- redcap.v4[,c(1, grep(pattern="visit[2-4].stool.age.months", colnames(redcap.v4)))]
colnames(age.at.stool) <- c("cid", "visit 2", "visit 3", "visit 4")

library(tidyr)
age.at.stool.v2 <- pivot_longer(data=age.at.stool, cols = c("visit 2", "visit 3", "visit 4"), names_to = "visit", values_drop_na = TRUE)
dim(age.at.stool.v2) #327 samples
str(age.at.stool.v2)
library(ggplot2)
ggplot(data=age.at.stool.v2, aes(x=value,fill=visit))+geom_dotplot()+scale_x_continuous(breaks=c(1.0,3.5,6.0,8.49))+scale_y_continuous(limits=c(0,1), labels = c(0,5,10,15,20))+xlab("infant age at stool collection (months)")+ylab("number of infant stool samples")
#save figure as "Denmark.infant.age.at.stool.collection.jpeg"

#STEP 4: calculate amount of time between subsequent stool samples for each infant
#and plot as histograms for time between visits 2->3 and visits 3->4
age.at.stool.v3 <- age.at.stool
age.at.stool.v3$`visits 2 - 3` <- age.at.stool.v3$`visit 3`-age.at.stool.v3$`visit 2`
age.at.stool.v3$`visits 3 - 4` <- age.at.stool.v3$`visit 4`-age.at.stool.v3$`visit 3`

#make into long format
age.at.stool.v4 <- pivot_longer(data=age.at.stool.v3, cols = c("visits 2 - 3", "visits 3 - 4"), names_to = "sample points", values_drop_na = TRUE)

ggplot(data=age.at.stool.v4, aes(x=value,fill=`sample points`))+geom_dotplot()+scale_x_continuous(breaks=c(1.0,1.5,2,2.5,3,3.5,4))+scale_y_continuous(labels = c(0,5,10,15,20))+xlab("time between stool collections (months)")+ylab("number of infant stool samples")

#STEP 5: look at stool storage and time passed until placement in freezer to estimate
#time stool was left at room temperature before freezer storage. 

#looking at visit 2 stool first
redcap.v4.visit2 <- redcap.v4[,c(1, grep(pattern="f211_doi|f211_dosc|f211_tosc|f211_locatn|f211_clinic|f211_hosp|f211_other|f211_hmfrez|f211_timefrez|f211_comments", colnames(redcap.v4)))]
redcap.v4$visit2.stool.days.at.room.temp <- "tbd"

#for the stool samples collected the same day as the interview, set the time
#at room temp to zero
redcap.v4[c(redcap.v4$f211_doi_q3==redcap.v4$f211_dosc_q5),"visit2.stool.days.at.room.temp"] <- "0"
table(redcap.v4$visit2.stool.days.at.room.temp)
#for the stool samples collected at the clinic, hospital, or trial site, set the
#time at room temp to zero as well since the sample was frozen on site. 
redcap.v4[c(redcap.v4$f211_locatn_q7==1|redcap.v4$f211_locatn_q7==2|redcap.v4$f211_locatn_q7==4),"visit2.stool.days.at.room.temp"] <- "0"
table(redcap.v4$visit2.stool.days.at.room.temp)

#for the stool samples collected at home, set the time at room temp to "0" if 
#the time of freezing is the same as the time of stool collection. 
redcap.v4[redcap.v4$f211_locatn_q7==3 & redcap.v4$f211_tosc_q6==redcap.v4$f211_timefrez_q7_5,"visit2.stool.days.at.room.temp"] <- "0" 
table(redcap.v4$visit2.stool.days.at.room.temp)

#among the stool samples that were collected at home but not frozen the same time
#as collection, calculating how long the stool sample was at room temp before freezing.
redcap.v4.tbd <- redcap.v4[redcap.v4$visit2.stool.days.at.room.temp=="tbd",c(1,grep(pattern="f211_tosc|f211_timefr|f211_doi|f211_dosc|f211_comments", colnames(redcap.v4)))]
redcap.v4.tbd$visit2.stool.room.temp.min <-difftime(time1=strptime(redcap.v4.tbd$f211_tosc_q6,format="%H:%M", tz=""), time2=strptime(redcap.v4.tbd$f211_timefrez_q7_5, format="%H:%M", tz=""), units="min")
range(redcap.v4.tbd$visit2.stool.room.temp.min)
#some of the samples don't have a time listed for when frozen
range(redcap.v4.tbd$visit2.stool.room.temp.min, na.rm=T)
#it was 5 to 15 minutes between when some of the stool samples were collected and 
#when they were frozen. That is short enough that I will assign those samples
#a duration of "0" hours between collection and freezing. 
redcap.v4[c(which(redcap.v4$mid %in% redcap.v4.tbd[!is.na(redcap.v4.tbd$visit2.stool.room.temp.min),][["mid"]])),"visit2.stool.days.at.room.temp"] <- "0"
table(redcap.v4$visit2.stool.days.at.room.temp)

#for those stool samples that were collected at home but apparently not frozen, 
#determining if they were actually frozen
redcap.v4[c(which(redcap.v4$mid %in% redcap.v4.tbd[is.na(redcap.v4.tbd$visit2.stool.room.temp.min),][["mid"]])),"f211_hmfrez_q7_4"]
#all answers were "NA" so probably means they were not stored in the freezer
#checking to see if the date of interview was close to the date of stool collection
#and also checking the comments list. 
redcap.v4[c(which(redcap.v4$mid %in% redcap.v4.tbd[is.na(redcap.v4.tbd$visit2.stool.room.temp.min),][["mid"]])),c(grep(pattern="^mid|f211_dosc|f211_doi|f211_hmfrez_q7_4|f211_comments", colnames(redcap.v4)))]

#I can see there are some stool collected more than one day before the date of 
#interview. So need to insert the number of days that the stool sample was potentially
#"left out" at room temperature. 
redcap.v4.tbd <- redcap.v4.tbd[is.na(redcap.v4.tbd$visit2.stool.room.temp.min),]
redcap.v4.tbd$f211_doi_q3 <- as.Date(redcap.v4.tbd$f211_doi_q3, format="%Y-%m-%d")
redcap.v4.tbd$f211_dosc_q5_v2 <- as.Date(redcap.v4.tbd$f211_dosc_q5_v2, format="%Y-%m-%d")

redcap.v4.tbd$visit2.stool.days.at.room.temp <-difftime(time1=redcap.v4.tbd$f211_doi_q3, time2=redcap.v4.tbd$f211_dosc_q5_v2, units="days")
range(redcap.v4.tbd$visit2.stool.days.at.room.temp)
#the range was -1 to 7 days. 
#looking specifically at the stool collected one day after the interview date.
print(redcap.v4.tbd[redcap.v4.tbd$visit2.stool.days.at.room.temp==-1,])
#according to the comments column for this sample, the stool was turned in with 
#a scale. So sounds like the stool sample was frozen the same day it was collected.
#Therefore I am going to put in "0" days as the time at room temp. 
redcap.v4.tbd[redcap.v4.tbd$visit2.stool.days.at.room.temp==-1,"visit2.stool.days.at.room.temp"] <- 0

#now putting in these time values into the original data frame under the column 
#header "visit2.stool.days.at.room.temp"
redcap.v4[c(which(redcap.v4$mid %in% redcap.v4.tbd$mid)),"visit2.stool.days.at.room.temp"] <- redcap.v4.tbd$visit2.stool.days.at.room.temp
redcap.v4$visit2.stool.days.at.room.temp
str(redcap.v4$visit2.stool.days.at.room.temp)

###now, applying the same steps to estimate the time spent at room temperature
#for visit 3 stool samples
#gathering the stool sample info into a small data frame to glance at some of the date trends
redcap.v4.visit3 <- redcap.v4[,c(1, grep(pattern="f311_doi|f311_dosc|f311_tosc|f311_locatn|f311_clinic|f311_hosp|f311_other|f311_hmfrez|f311_timefrez|f311_comments", colnames(redcap.v4)))]
#look at the stool samples that were not collected the same day as the interview 
# because the stool collected the same day as the interview would have been frozen
#at the trial center right away. 
redcap.v4$visit3.stool.days.at.room.temp <- "tbd"

#for the stool samples collected the same day as the interview, set the time
#at room temp to zero
redcap.v4[c(redcap.v4$f311_doi_q3==redcap.v4$f311_dosc_q5),"visit3.stool.days.at.room.temp"] <- "0"
table(redcap.v4$visit3.stool.days.at.room.temp)
#for the stool samples collected at the clinic, hospital, or trial site, set the
#time at room temp to zero as well since the sample was frozen on site. 
redcap.v4[c(redcap.v4$f311_locatn_q7==1|redcap.v4$f311_locatn_q7==2|redcap.v4$f311_locatn_q7==4),"visit3.stool.days.at.room.temp"] <- "0"
table(redcap.v4$visit3.stool.days.at.room.temp)

#for the stool samples collected at home, set the time at room temp to "0" if 
#the time of freezing is the same as the time of stool collection. 
redcap.v4[redcap.v4$f311_locatn_q7==3 & redcap.v4$f311_tosc_q6==redcap.v4$f311_timefrez_q7_5,"visit3.stool.days.at.room.temp"] <- "0" 
table(redcap.v4$visit3.stool.days.at.room.temp)

#among the stool samples that were collected at home but not frozen the same time
#as collection, calculating how long the stool sample was at room temp before freezing.
redcap.v4.tbd <- redcap.v4[redcap.v4$visit3.stool.days.at.room.temp=="tbd",c(1,grep(pattern="f311_tosc|f311_timefr|f311_doi|f311_dosc|f311_comments", colnames(redcap.v4)))]
redcap.v4.tbd$visit3.stool.room.temp.min <-difftime(time1=strptime(redcap.v4.tbd$f311_tosc_q6,format="%H:%M", tz=""), time2=strptime(redcap.v4.tbd$f311_timefrez_q7_5, format="%H:%M", tz=""), units="min")
range(redcap.v4.tbd$visit3.stool.room.temp.min)
#some of the samples don't have a time listed for when frozen
range(redcap.v4.tbd$visit3.stool.room.temp.min, na.rm=T)
#it was 5 to 15 minutes between when some of the stool samples were collected and 
#when they were frozen. That is short enough that I will assign those samples
#a duration of "0" days between collection and freezing. 
redcap.v4[c(which(redcap.v4$mid %in% redcap.v4.tbd[!is.na(redcap.v4.tbd$visit3.stool.room.temp.min),][["mid"]])),"visit3.stool.days.at.room.temp"] <- "0"
table(redcap.v4$visit3.stool.days.at.room.temp)

#for those stool samples that were collected at home but apparently not frozen, 
#determining if they were actually frozen
redcap.v4[c(which(redcap.v4$mid %in% redcap.v4.tbd[is.na(redcap.v4.tbd$visit3.stool.room.temp.min),][["mid"]])),"f311_hmfrez_q7_4"]
#all answers were "NA" so probably means they were not stored in the freezer
#checking to see if the date of interview was close to the date of stool collection
#and also checking the comments list. 
redcap.v4[c(which(redcap.v4$mid %in% redcap.v4.tbd[is.na(redcap.v4.tbd$visit3.stool.room.temp.min),][["mid"]])),c(grep(pattern="^mid|f311_dosc|f311_doi|f311_hmfrez_q7_4|f311_comments", colnames(redcap.v4)))]

#I can see there is a stool collected more than one day before the date of 
#interview. So need to insert the number of days that the stool sample was potentially
#"left out" at room temperature. 
redcap.v4.tbd <- redcap.v4.tbd[is.na(redcap.v4.tbd$visit3.stool.room.temp.min),]
redcap.v4.tbd$f311_doi_q3 <- as.Date(redcap.v4.tbd$f311_doi_q3, format="%Y-%m-%d")
redcap.v4.tbd$f311_dosc_q5_v2 <- as.Date(redcap.v4.tbd$f311_dosc_q5_v2, format="%Y-%m-%d")

redcap.v4.tbd$visit3.stool.days.at.room.temp <-difftime(time1=redcap.v4.tbd$f311_doi_q3, time2=redcap.v4.tbd$f311_dosc_q5_v2, units="days")
range(redcap.v4.tbd$visit3.stool.days.at.room.temp)
#the range was 1 to 6 days. 

#now putting in these time values into the original data frame under the column 
#header "visit3.stool.days.at.room.temp"
redcap.v4[c(which(redcap.v4$mid %in% redcap.v4.tbd$mid)),"visit3.stool.days.at.room.temp"] <- redcap.v4.tbd$visit3.stool.days.at.room.temp
redcap.v4$visit3.stool.days.at.room.temp
str(redcap.v4$visit3.stool.days.at.room.temp)

###now, applying the same steps to estimate the time spent at room temperature
#for visit 4 stool samples
#gathering the stool sample info into a small data frame to glance at some of the date trends
redcap.v4.visit4 <- redcap.v4[,c(1, grep(pattern="f411_doi|f411_dosc|f411_tosc|f411_locatn|f411_clinic|f411_hosp|f411_other|f411_hmfrez|f411_timefrez|f411_comments", colnames(redcap.v4)))]
#look at the stool samples that were not collected the same day as the interview 
# because the stool collected the same day as the interview would have been frozen
#at the trial center right away. 
redcap.v4$visit4.stool.days.at.room.temp <- "tbd"

#for the stool samples collected the same day as the interview, set the time
#at room temp to zero
redcap.v4[c(redcap.v4$f411_doi_q3==redcap.v4$f411_dosc_q5),"visit4.stool.days.at.room.temp"] <- "0"
table(redcap.v4$visit4.stool.days.at.room.temp)
#for the stool samples collected at the clinic, hospital, or trial site, set the
#time at room temp to zero as well since the sample was frozen on site. 
redcap.v4[c(redcap.v4$f411_locatn_q7==1|redcap.v4$f411_locatn_q7==2|redcap.v4$f411_locatn_q7==4),"visit4.stool.days.at.room.temp"] <- "0"
table(redcap.v4$visit4.stool.days.at.room.temp)

#for the stool samples collected at home, set the time at room temp to "0" if 
#the time of freezing is the same as the time of stool collection. 
redcap.v4[redcap.v4$f411_locatn_q7==3 & redcap.v4$f411_tosc_q6==redcap.v4$f411_timefrez_q7_5,"visit4.stool.days.at.room.temp"] <- "0" 
table(redcap.v4$visit4.stool.days.at.room.temp)

#among the stool samples that were collected at home but not frozen the same time
#as collection, calculating how long the stool sample was at room temp before freezing.
redcap.v4.tbd <- redcap.v4[redcap.v4$visit4.stool.days.at.room.temp=="tbd",c(1,grep(pattern="f411_tosc|f411_timefr|f411_doi|f411_dosc|f411_comments", colnames(redcap.v4)))]
redcap.v4.tbd$visit4.stool.room.temp.min <-difftime(time1=strptime(redcap.v4.tbd$f411_tosc_q6,format="%H:%M", tz=""), time2=strptime(redcap.v4.tbd$f411_timefrez_q7_5, format="%H:%M", tz=""), units="min")
range(redcap.v4.tbd$visit4.stool.room.temp.min)
#some of the samples don't have a time listed for when frozen
range(redcap.v4.tbd$visit4.stool.room.temp.min, na.rm=T)
#it was 0 to 15 minutes between when some of the stool samples were collected and 
#when they were frozen. That is short enough that I will assign those samples
#a duration of "0" days between collection and freezing. 
redcap.v4[c(which(redcap.v4$mid %in% redcap.v4.tbd[!is.na(redcap.v4.tbd$visit4.stool.room.temp.min),][["mid"]])),"visit4.stool.days.at.room.temp"] <- "0"
table(redcap.v4$visit4.stool.days.at.room.temp)

#for those stool samples that were collected at home but apparently not frozen, 
#determining if they were actually frozen
redcap.v4[c(which(redcap.v4$mid %in% redcap.v4.tbd[is.na(redcap.v4.tbd$visit4.stool.room.temp.min),][["mid"]])),"f411_hmfrez_q7_4"]
#all answers were "1" so means they were stored in the freezer
#checking to see if the date of interview was close to the date of stool collection
#and also checking the comments list. 
redcap.v4[c(which(redcap.v4$mid %in% redcap.v4.tbd[is.na(redcap.v4.tbd$visit4.stool.room.temp.min),][["mid"]])),c(grep(pattern="^mid|f411_dosc|f411_doi|f411_tosc|f411_timefr|f411_hmfrez_q7_4|f411_comments", colnames(redcap.v4)))]

#it looks like the samples were all stored in the freezer soon or right after collection. 
#so setting number of days at room temp for these stool samples to "0". Now putting
#in these time values into the original data frame under the column header 
#"visit4.stool.days.at.room.temp" 
redcap.v4[c(which(redcap.v4$mid %in% redcap.v4.tbd$mid)),"visit4.stool.days.at.room.temp"] <- "0"
redcap.v4$visit4.stool.days.at.room.temp
str(redcap.v4$visit4.stool.days.at.room.temp)

#STEP 6: see if the stool samples were scored for the stool consistency. 
table(redcap.v4$f211_const_q8_2)
sum(is.na(redcap.v4$f211_const_q8_2)) #one visit 2 stool sample without consistency score 
table(redcap.v4$f311_const_q8_2)
sum(is.na(redcap.v4$f311_const_q8_2)) #six visit 3 stool samples without consistency score
table(redcap.v4$f411_const_q8_2)
sum(is.na(redcap.v4$f411_const_q8_2)) #two visit 3 stool samples without consistency score

#STEP 7: match up the sample.id's with the respective stool sample information
#and write to file for use with other analyses. 
stool.info <- redcap.v4[,c(1,grep(pattern="visit[2-4][.]stool|f[2-4]11_const_q8_2", colnames(redcap.v4)))]
stool.info.vt2 <- stool.info[,c(1, grep("f2|visit2", colnames(stool.info)))]
colnames(stool.info.vt2) <- c("mid","stool_consistency_const_q8_2", "stool_age_days", "stool_age_months", "stool_days_at_room_temperature")
stool.info.vt2$visit <- "2"
stool.info.vt3 <- stool.info[,c(1, grep("f3|visit3", colnames(stool.info)))]
colnames(stool.info.vt3) <- c("mid","stool_consistency_const_q8_2", "stool_age_days", "stool_age_months", "stool_days_at_room_temperature")
stool.info.vt3$visit <- "3"
stool.info.vt4 <- stool.info[,c(1, grep("f4|visit4", colnames(stool.info)))]
colnames(stool.info.vt4) <- c("mid","stool_consistency_const_q8_2", "stool_age_days", "stool_age_months", "stool_days_at_room_temperature")
stool.info.vt4$visit <- "4"

stool.info.v2 <- rbind(stool.info.vt2, stool.info.vt3, stool.info.vt4)
metadata.v2 <- merge(metadata, stool.info.v2, by=c("mid", "visit"), all=F)

#write to file
write.csv(metadata.v2, "metadata-327-samples-with-select-stool-information.csv")
