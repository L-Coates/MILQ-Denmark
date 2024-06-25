#This script describes the steps taken and code involved in summarizing
#demographic and certain metadata features of the Denmark infants among those
#being used for analyses. 

#STEP 1: read in the metadata file of the infant samples being used for analysis, 
#and use to select REDCap data. 

metadata <- read.csv("../metadata-infants-complete-stool-set-after-6886-read-cutoff-withoutInfantCD144Y.csv")
#make a maternal ID column to match up with the redcap dataset
metadata$mid <- metadata$extraction.id
metadata$mid <- gsub(pattern="^C", replacement="M", metadata$mid)
metadata$mid <- gsub(pattern="[.][2-4][x|b]*", replacement="", metadata$mid)

redcap <- read.csv("../../python_scripts/translated_Denmark_REDCap_metadata.csv")
#we know a couple maternal id's have a lowercase check letter which needs to be
#switched to uppercase for purposes of matching. 
redcap[redcap$mid=="MD346m","mid"] <- "MD346M"
redcap[redcap$mid=="MD378x","mid"] <- "MD378X"
redcap <- redcap[,-1]

redcap.v2 <- redcap[c(which(redcap$mid %in% metadata$mid)),]

#STEP 2: look at birth mode, infant gender, parity
unique(redcap.v2$redcap_event_name)
redcap.v2.f201 <- redcap.v2[redcap.v2$redcap_event_name=="10__349_months_arm_1",c(1,grep(pattern="f201|f101", colnames(redcap.v2)))]
length(unique(redcap.v2.f201$mid))
table(redcap.v2.f201$f201_mod_q11)
#102 vaginal deliveries, and 7 c-section deliveries. ("1" = vaginal, "2"=planned c-section, "3" = emergency c-section) 
table(redcap.v2.f201$f201_infgend_q10)
#47 boys, 62 girls ("1"=boy, "2"=girl)

table(redcap.v2$f003_frstpreg_q6) #"1" = yes (first pregnancy), "0" = no (not first pregnancy)
#for 80 of the mothers this was the first pregnancy, for 29 mothers this was not 
#the first pregnancy

table(redcap.v2$f003_chdrnbrn_q6_1)
#for those 29 mothers, 24 had one other child, and 5 had two other children.

#look at number of people living in household
table(redcap.v2$f003_chdrn_q8+redcap.v2$f003_adults_q7)
#for 1 mother, it's just one two people in household (mom and baby)
#for 71 mothers, it's three people in household. 
#for 28 mothers, it's four people in household. 
#for 8 mothers, it's five people in household. 
#for 1 mother, it's six people in household. 

#STEP 3: look at some of mother's conditions like pre-eclampsia and gestational diabetes
#and infant pre-term birth
table(redcap.v2$f101_preec_q24) #none had pre-eclampsia
table(redcap.v2$f101_gestdiab_q25) #none had gestational diabetes
table(redcap.v2$f101_preterm_q5_3) #none were pre-term

#STEP 4: match up to stool sample IDs and write to file. 

redcap.data.of.interest <- redcap.v2.f201[, c(grep(pattern="^mid|f201_mod_q11$|f201_infgend_q10$", colnames(redcap.v2.f201)))]
colnames(redcap.data.of.interest) <- c("mid", "infant.gender", "mode.of.delivery")
redcap.data.of.interest$infant.gender <- gsub(pattern="1", replacement="boy", redcap.data.of.interest$infant.gender)
redcap.data.of.interest$infant.gender <- gsub(pattern="2", replacement="girl", redcap.data.of.interest$infant.gender)

redcap.data.of.interest$mode.of.delivery <- gsub(pattern="1", replacement="vaginal", redcap.data.of.interest$mode.of.delivery)
redcap.data.of.interest$mode.of.delivery <- gsub(pattern="2|3", replacement="c-section", redcap.data.of.interest$mode.of.delivery)

redcap.visit2 <- redcap.v2[redcap.v2$redcap_event_name=="10__349_months_arm_1",]
redcap.visit2[is.na(redcap.visit2$f003_chdrnbrn_q6_1), "f003_chdrnbrn_q6_1"] <- 0
redcap.visit2$parity <- redcap.visit2$f003_frstpreg_q6 + redcap.visit2$f003_chdrnbrn_q6_1

#also create a categorical variable called "first.pregnancy" with values of "yes"
#for parity = 1, and "no" for parity >1
redcap.visit2$first.pregnancy <- ""
redcap.visit2[redcap.visit2$parity=="1", "first.pregnancy"] <- "yes"
redcap.visit2[redcap.visit2$parity!="1", "first.pregnancy"] <- "no"
table(redcap.visit2$first.pregnancy)

#look at number of people living in household and create a new categorical variable
#called "household.size" which has the value "<= 3" for households with 3 or fewer people,
#and the value ">3" for households with more than three people. 
redcap.visit2$number.people.in.household <- redcap.visit2$f003_chdrn_q8+redcap.visit2$f003_adults_q7
redcap.visit2$household.size <- ""
redcap.visit2[redcap.visit2$number.people.in.household<=3,"household.size"] <- "<= 3"
redcap.visit2[redcap.visit2$number.people.in.household>3,"household.size"] <- "> 3"
table(redcap.visit2$household.size)

redcap.data.of.interest <- merge(x=redcap.data.of.interest, y=redcap.visit2[, c(grep(pattern="^mid|parity|number.people.in.household|first|household", colnames(redcap.visit2)))], by="mid", all=F)
metadata.v2 <- merge(x=metadata, y=redcap.data.of.interest, by="mid", all=F)

#change order of columns and remove the placeholder column
metadata.v3 <- metadata.v2[,c(3:27,1)]

write.csv(metadata.v3, "birthmode.gender.parity.csv", row.names = F)

