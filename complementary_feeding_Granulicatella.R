#This script contains the code used to investigate associations between 
#age of complementary feeding and Granulicatella abundance.

library(ggplot2)
library(ggpubr)
library(tidyverse)

#PART 1: determine infant age at which something other than breastmilk was given
#to the infant. 

#step 1: read in all the relevant metadata

#read in the metadata file with the list of infant stool samples included in analyses
metadata <- read.csv("../morbidity_medication_prevalences/infant.morbidities.medications.csv", header=T)
metadata <- metadata[,-1]


#extract the 109 infant IDs 
cids <- metadata$extraction.id
cids <- cids[grep(pattern="^CD", cids)]
cids <- gsub(pattern="x$|b$", replacement="", cids)
cids <- gsub(pattern="[.][234]", replacement="", cids)
cids <- unique(cids)
mids <- gsub(pattern="^C", replacement="M", cids)

#read in the REDCap data
redcap <- read.csv("../../../../REDCap_downloads/Denmark/MILQMAINSTUDYDENMARK_DATA_2022-10-11_2251.csv")

#correct two maternal IDs that are known to have lower case letters
redcap[redcap$mid=="MD346m","mid"] <- "MD346M"
redcap[redcap$mid=="MD378x","mid"] <- "MD378X"

#extract the REDCap data for just the subset of 109 infants
redcap.v2 <- redcap[which(redcap$mid %in% mids),]
length(unique(redcap.v2$mid)) #109 infants

#change data frame so each infant only has one row
redcap.v2.visit1 <- redcap.v2[redcap.v2$redcap_event_name=="2__3_days_postpart_arm_1",]
redcap.v2.visit1 <- redcap.v2.visit1[, grep(pattern="^mid|f101", colnames(redcap.v2.visit1))]

redcap.v2.visit2 <- redcap.v2[redcap.v2$redcap_event_name=="10__349_months_arm_1",]
redcap.v2.visit2 <- redcap.v2.visit2[, grep(pattern="^mid|f201", colnames(redcap.v2.visit2))]

redcap.v2.visit3 <- redcap.v2[redcap.v2$redcap_event_name=="35__599_months_arm_1",]
redcap.v2.visit3 <- redcap.v2.visit3[, grep(pattern="^mid|f306", colnames(redcap.v2.visit3))]

redcap.v2.visit4 <- redcap.v2[redcap.v2$redcap_event_name=="60__849_months_arm_1",]
redcap.v2.visit4 <- redcap.v2.visit4[, grep(pattern="^mid|f406", colnames(redcap.v2.visit4))]

redcap.v3 <- merge(redcap.v2.visit1, redcap.v2.visit2, by="mid", all=T)
redcap.v3 <- merge(redcap.v3, redcap.v2.visit3, by="mid", all=T)
redcap.v3 <- merge(redcap.v3, redcap.v2.visit4, by="mid", all=T)

#step 2: identify infants that got something other than breastmilk (could be water
#or formula) during the first 3.5 months of life

#see if any mothers started giving something other than breastmilk to their babies
#at some point before 3.5 months. 
table(redcap.v3$f201_infconsump_q6) 
#4 babies did receive some drink/food other than breastmilk, while 105 did not.
print(redcap.v3[redcap.v3$f201_infconsump_q6=="1", "mid"])
#mothers "MD078U" "MD124Q" "MD183V" "MD261C"

#looking at whether it was consumed in the first week of life only. 
table(redcap.v3$f201_consumfw_q6_3) 
#no it was not exclusive to the first week of life -- it was consumed at other times
#outside of the first week of life.

#was it infant formula? 
table(redcap.v3$f201_infform_q6_1)
redcap.v3[redcap.v3$f201_infform_q6_1=="1","mid"]
#two babies (MD183V and MD261C) had infant formula, 2 babies did not. 

#looking at what else was given to the infant.
table(redcap.v3$f201_other_q6_2_1)
#"a bit of fennel tea", and "sugar water approximately 5 times during the first 3 weeks of life (max 5 ml at a time)"
redcap.v3[c(redcap.v3$f201_other_q6_2_1!=""), c("mid", "f201_other_q6_2_1")]
#MD078U received the fennel tea, and MD124Q received the sugar water. 

#look at whether this matches up with responses from questionnaire 1.01 ("newborn screening information")
#here they ask if the infant has received anything other than breastmilk since birth
#(not just outside of the first week of life.)
table(redcap.v3$f101_infconsother_q15)
#20 infants received something other than breastmilk early on. 
#89 infants did not. 
redcap.v3[redcap.v3$f101_infconsother_q15=="1","mid"]
#"MD007M" "MD014Y" "MD018W" "MD050Q" "MD063M" "MD080P" "MD093L" "MD130J" "MD140A" "MD146K" "MD183V" "MD273F" "MD293N" "MD296S"
#"MD299X" "MD312P" "MD337O" "MD339A" "MD356D" "MD382E"

#looking at whether it was given in the first seven days after birth. 
table(redcap.v3$f101_aftbirth_q15_3)
identical(unique(redcap.v3[redcap.v3$f101_infconsother_q15=="1","mid"]),unique(redcap.v3[c(redcap.v3$f101_aftbirth_q15_3=="1"& !is.na(redcap.v3$f101_aftbirth_q15_3)),"mid"]))
#all 20 infants had it in the first week of life. 

#did it include breast milk substitute?
table(redcap.v3$f101_infform_q15_1) 
#17 babies had breast milk substitute, and 3 had something else. 
redcap.v3[(redcap.v3$f101_infform_q15_1=="1" & !is.na(redcap.v3$f101_infform_q15_1)), "f101_infform_q15_1_1"]
#most didn't provide more details, some said it was only given right after birth, 
#and some said the types of formula. 
redcap.v3[(redcap.v3$f101_infform_q15_1=="1" & !is.na(redcap.v3$f101_infform_q15_1)), "f101_other_q15_2_1"]
#a baby that had formula also had donor milk. 

#looking at what those three babies received, other than breast milk substitute..
redcap.v3[(redcap.v3$f101_infform_q15_1=="0" & !is.na(redcap.v3$f101_infform_q15_1)), c("f101_other_q15_2_1", "mid")]
redcap.v3[(redcap.v3$f101_infform_q15_1=="0" & !is.na(redcap.v3$f101_infform_q15_1)), "f101_other_q15_2_2"]
#sugar water and donor milk

#making a list of these 19 infants that received something other than breastmilk (or donor milk)
#in the first week of life. This means that infant MD050Q which received donor milk
#would not be among the list. 
comp.food.first.week <- c(redcap.v3[(redcap.v3$f101_infform_q15_1=="1" & !is.na(redcap.v3$f101_infform_q15_1)), "mid"],"MD014Y", "MD130J")

#see if any mothers started giving something other than breast milk to their babies
#at some point in 3.5 - 5.99 months. 
table(redcap.v3$f306_otherfd_q10) 
#55 babies did receive some drink/food other than breastmilk, while 54 did not. 

#does this include the four babies that apparently had something other than breastmilk
#during visit 2. 
redcap.v3[redcap.v3$f306_otherfd_q10==1,"mid"]
#yes, those four babies were among the 55

#looking at whether infants had beverages and food
table(redcap.v3$f306_liq_q11) #49 had liquid, 5 did not. 
table(redcap.v3$f306_semisolid_q12) #46 had semisolids, 9 did not. 
#so most of these babies that had something, had something that was a liquid and 
#also something that was a semisolid. 

#other comments?
table(redcap.v3$f306_comments_q13)

#step 3: look at what time frame other infants started consuming something other 
#than breastmilk. 

#see if all mothers started giving something other than breast milk by visit 4 (6 - 8.49 months)
table(redcap.v3$f406_otherfd_q10) #all mothers answered yes to this question. 

#looking at whether infants had beverages and food
table(redcap.v3$f406_liq_q11) #108 had liquid, 1 did not. 
table(redcap.v3$f406_semisolid_q12) #108 had semisolids, 1 did not. 
#so most of these babies that had something, had something that was a liquid and 
#also something that was a semisolid. 

#other comments?
table(redcap.v3$f406_comments_q13)

#step 4: create a variable which is the earliest age at which an infant received
#anything other than breastmilk. 

#gather the ages listed for when each food type was introduced and select the
#youngest age to mark as the start of complementary feeding. 
age <- redcap.v3[,grep(pattern="^mid|f306_agemnths|f306_agewk|f406_agemnths|f406_agewk", colnames(redcap.v3))]

#write a function to calculate the age in weeks of the infant when each liquid/food
#type was introduced. This involves changing the age in months to weeks, and then
#adding the weeks column (which is a value between 0-3). 

age.months = c(grep(pattern="agemnths", colnames(age)))
calculate.age <- function(age.months){
    age.v2 <- data.frame(age[1])
    for (i in age.months) {
        age
        collapsed <- data.frame((as.numeric(age[[i]])*4.345)+as.numeric(age[[(i+1)]]))
        colnames(collapsed) <- colnames(age[,c(i,i+1)])[2]
        age.v2<- cbind(age.v2, collapsed)
    }
    print(age.v2)
}

age.v2 <- calculate.age(age.months)

age.v2$complementary.feeding.age.weeks <- apply(X=age.v2[,c(grep(pattern="agewks", colnames(age.v2)))], MARGIN=1, FUN=min, na.rm=TRUE)

age.v2$complementary.feeding.age.weeks.notincluding.water <- apply(X=age.v2[,-c(grep(pattern="f306_agewks_q11_1_1b$|f406_agewks_q11_1_1b$|mid$|complementary", colnames(age.v2)))], MARGIN=1, FUN=min, na.rm=TRUE)

#MD238G does not have an age of complementary feeding listed. Looking at whether 
#this infant was ever given complementary foods during the time of study. 
redcap.v3[redcap.v3$mid=="MD238G",] 
#the infant was introduced to foods but didn't consume any, and it was later in infancy (6-8.5 months). 
#so setting the age of complementary feeding to NA
age.v2[age.v2$mid=="MD238G", "complementary.feeding.age.weeks"] <- NA
age.v2[age.v2$mid=="MD238G", "complementary.feeding.age.weeks.notincluding.water"] <- NA

#printout out the names of the infants that received something other than breastmilk
#before the age of 3.5 months (i.e. 15.203 weeks) and confirming, "manually", that the 
#responses show they were given food before 3.5 months. 
print(age.v2[age.v2$complementary.feeding.age.weeks<=15.203,"mid"])
#"MD007M" "MD050Q" "MD119Q" "MD124Q" "MD140A" "MD223Y" "MD259H""MD261C" "MD335C" "MD378X"
#there are two babies (MD078U, and MD183V) that had something other than breast milk during visit 2, but
#later on their responses didn't reflect this. 

#looking at what those two babies' ages were for the specific responses
print(age.v2[age.v2$mid=="MD078U"|age.v2$mid=="MD183V", "complementary.feeding.age.weeks"])
#18.38 weeks and 19.38 weeks were their ages and these ages are greater than 3.5 months, 
#so these ages of food introduction need to be manually changed according to their 
#visit 2 responses about complementary food consumption. 

#for the 19 infants (+ MD124Q) that received something other than
#breastmilk during the first week of life, changing the age in the "age.v2" data frame
#to 1 week. And for the infant MD078U which received fennel tea in visit 2 (but not
#in the first week of life and also no age specified within visit 2), setting the age
#of complementary feeding to the date of the interview for form 2.01. 
age.v2[c(which(age.v2$mid %in% c(comp.food.first.week, "MD124Q"))),"complementary.feeding.age.weeks"] <- 1
age.v2[c(which(age.v2$mid %in% c(comp.food.first.week, "MD124Q"))),"complementary.feeding.age.weeks.notincluding.water"] <- 1
redcap[redcap$mid=="MD078U", c("f101_dob_q5","f201_doi_q3")]
difftime("2018-05-09","2018-08-13")
#infant MD078U was 96 days old at time of interview for f2.01. using 30.4375 days/month
#4.345 weeks/month, the infants age in weeks at the time of interview was 13.7. 
age.v2[age.v2$mid=="MD078U","complementary.feeding.age.weeks"] <- 13.7
age.v2[age.v2$mid=="MD078U","complementary.feeding.age.weeks.notincluding.water"] <- 13.7

#looking at the range of ages when complementary feeding started (or at least when 
#something other than breastmilk was given). 
summary(age.v2$complementary.feeding.age.weeks)
summary(age.v2$complementary.feeding.age.weeks.notincluding.water)
#some babies got something other than breastmilk before 3.5 months of age, and it
#was outside the window of the first week of life when giving water/sugar water/formula is
#permitted in this study. 

#creating a variable that is "ExclusiveBreastfeedingBefore4months" and 
#writing to file. 4.345 * 4. Note: infant MD238G had "NA" for age at complementary feeding
#because it never consumed any during the study, so giving this infant the value "yes"
#for exclusive breastfeeding before 4 months of age.
complementary.feeding.before.after.4months <- age.v2[,grep(pattern="^mid|complemen", colnames(age.v2))]
complementary.feeding.before.after.4months$ExclusiveBreastfeeding.before.4months <-""
complementary.feeding.before.after.4months[(complementary.feeding.before.after.4months$complementary.feeding.age.weeks<17.38)&!is.na(complementary.feeding.before.after.4months$complementary.feeding.age.weeks), "ExclusiveBreastfeeding.before.4months"] <-"no"
complementary.feeding.before.after.4months[complementary.feeding.before.after.4months$complementary.feeding.age.weeks>=17.38 &!is.na(complementary.feeding.before.after.4months$complementary.feeding.age.weeks), "ExclusiveBreastfeeding.before.4months"] <-"yes"
complementary.feeding.before.after.4months[complementary.feeding.before.after.4months$mid=="MD238G", "ExclusiveBreastfeeding.before.4months"] <-"yes"
table(complementary.feeding.before.after.4months$ExclusiveBreastfeeding.before.4months)

write.csv(complementary.feeding.before.after.4months, row.names=F, file="Exclusive.breastfeeding.before.4months.csv")

#calculate the age in days at time of complementary feeding, then
#make a histogram of age (in days) at complementary feeding. 
age.v2$complementary.feeding.age.days <- (age.v2$complementary.feeding.age.weeks)*7
age.v2$complementary.feeding.age.days.notincluding.water <- (age.v2$complementary.feeding.age.weeks.notincluding.water)*7
summary(age.v2$complementary.feeding.age.days)
summary(age.v2$complementary.feeding.age.days.notincluding.water)
jpeg(filename="histogram.age.complementary.feeding.days.jpeg")
hist(age.v2$complementary.feeding.age.days, breaks=20, xlab="age at first complementary fluid/food (days)", main="")
dev.off()

#PART 2: explore the potential relationship between age of complementary feeding
#and Granulicatella abundance. 

#Step 1: make a figure of Granulicatella abundance (rarefied counts) over time 
#for each infant with coloring of the lines or points by age of complementary feeding. 

library(qiime2R)
library(phyloseq)

#read in the metadata file with the 327 stool samples to analyze. 
metadata <- read.csv("../morbidity_medication_prevalences/infant.morbidities.medications.csv")
metadata <- metadata[,-1]

#add the info related to complementary feeding to the metadata
metadata.v2 <- merge(x=metadata, y=age.v2, by="mid", all=F)

#reading in the rarefied counts of Granulicatella
feature.table.rarefied <- read_qza("../feature-table-rarefied-6963.qza")
feature.table.rarefied <- feature.table.rarefied$data

taxonomy <- read_qza("../taxonomy-classification.qza")
taxonomy2 <- as.data.frame(taxonomy$data) %>% column_to_rownames("Feature.ID") %>% separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>% as.matrix()
#remove "Confidence" column
taxonomy2 <- taxonomy2[,-8]

#trim taxonomy to include only ASVs in the ASV table
taxonomy3 <- taxonomy2[c(rownames(taxonomy2)%in%rownames(feature.table.rarefied)),]

phylobj.rarefied <- phyloseq(otu_table(feature.table.rarefied, taxa_are_rows=T),tax_table(taxonomy3), sample_data(metadata.v2 %>% as.data.frame() %>% column_to_rownames("sample.id")))
#get rid of the taxa that are not present at all
phylobj.rarefied<- filter_taxa(phylobj.rarefied, function(x) sum(x)>0, TRUE)
#look at the number of ASVs among the feature table
dim(otu_table(phylobj.rarefied)) #1585 ASVs among 327 samples

#collapse the counts to genus level designations. 
phylobj.rarefied.genus <- tax_glom(phylobj.rarefied, 'Genus', NArm=F)

melted.phylgenus.rarefied <- psmelt(phylobj.rarefied.genus)
unique(melted.phylgenus.rarefied$Genus)
granulicatella.rarefied <- melted.phylgenus.rarefied[c(grep(pattern="Granulica",melted.phylgenus.rarefied$Genus)),]
table(granulicatella.rarefied$Abundance>0) #60 samples had a count above zero. 

granulicatella.rarefied$Relative.Abundance <- (granulicatella.rarefied$Abundance/6963)*100

library(ggpubr)
ggplot(aes(x=stool_age_days, y=Relative.Abundance,group=mid, col=complementary.feeding.age.days), data=granulicatella.rarefied)+
    geom_line()+
    scale_color_gradient2(low="red", mid="green", high="blue", midpoint=(3.49*30.4375))+
    geom_point()+
    theme_bw()+
    xlab("age at stool collection (days)")+
    ylab("Granulicatella relative abundance %")+
    labs(color="age at first consumption of\nnon-breastmilk fluid/food (days)")
ggsave("age.complementary.feeding.vs.granulicatella.jpeg", plot=last_plot(), dpi=600, height=6,width=10, units="in")

ggplot(aes(x=stool_age_days, y=Relative.Abundance,group=mid, col=complementary.feeding.age.days.notincluding.water), data=granulicatella.rarefied)+
    geom_line()+
    scale_color_gradient2(low="red", mid="green", high="blue", midpoint=(3.49*30.4375))+
    geom_point()+
    theme_bw()+
    xlab("age at stool collection (days)")+
    ylab("Granulicatella relative abundance %")+
    labs(color="age of fluid/food introduction,\nnot including water (days)")

#STEP 5: make a box plot of age of complementary feeding for infants that didn't
#have Granulicatella above 0.00% versus infants that did have Granulicatella. 

total.granulicatella.counts.per.infants <- data.frame(tapply(X=granulicatella.rarefied$Abundance, INDEX=granulicatella.rarefied$mid, FUN=sum))
total.granulicatella.counts.per.infants$mid <- row.names(total.granulicatella.counts.per.infants)
colnames(total.granulicatella.counts.per.infants) <- c("total.granulicatella.counts", "mid")

infants.with.granulicatella <- total.granulicatella.counts.per.infants[total.granulicatella.counts.per.infants$total.granulicatella.counts>0, "mid"]
infants.without.granulicatella <- total.granulicatella.counts.per.infants[total.granulicatella.counts.per.infants$total.granulicatella.counts==0, "mid"]

age.v3 <- age.v2
age.v3$Granulicatella.detected <- ""
age.v3[c(which(age.v3$mid %in% infants.with.granulicatella)), "Granulicatella.detected"] <- "yes"
age.v3[c(which(age.v3$mid %in% infants.without.granulicatella)), "Granulicatella.detected"] <- "no"

#is there a statistically significant difference in age of complementary fluid/feed introduction
#for infants that had Granulicatella detected in their stool versus those that did not? 
wilcox.test(x=age.v3[age.v3$Granulicatella.detected=="no","complementary.feeding.age.days"], y=age.v3[age.v3$Granulicatella.detected=="yes","complementary.feeding.age.days"], paired=F)
#not significant difference
wilcox.test(x=age.v3[age.v3$Granulicatella.detected=="no","complementary.feeding.age.days"], y=age.v3[age.v3$Granulicatella.detected=="yes","complementary.feeding.age.days"], paired=F, alternative = "less")
#lower age of complementary fluids/food introduction among the infants that didn't
#have Granulicatella detected in their stool. 

ggplot(aes(x=Granulicatella.detected, y=complementary.feeding.age.days), data=age.v3) +
    geom_boxplot()+
    theme_bw()+
    geom_jitter()+
    xlab("Granulicatella relative abundance > 0.00 %")+ylab("age at first consumption of\nnon-breastmilk fluid/food (days)")
ggsave("age.complementary.feeding.vs.Granulicatella.detection.jpeg", plot=last_plot())

ggplot(aes(x=Granulicatella.detected, y=complementary.feeding.age.notincluding.water.days), data=age.v3) +
    geom_boxplot()+
    theme_bw()+
    geom_jitter() +
    xlab("Granulicatella relative abundance > 0.00 %")+
    ylab("age at first consumption of\nnon-breastmilk fluid/food (days)")

#making a grid of two plots: 1) the earlier scatter plot showing Granulicatella
#abundance over time with coloring by age of complementary feeding; and 2) this
#recent box plot of age at complementary feeding in infants that never had Granulicatella
#in their stool and those that did at least in one stool. 
a <- ggplot(aes(x=stool_age_days, y=Relative.Abundance,group=mid, col=complementary.feeding.age.days), data=granulicatella.rarefied)+
    geom_line()+
    theme_bw()+
    scale_color_gradient2(low="red", mid="green", high="blue", midpoint=(3.49*30.4375))+
    geom_point()+
    xlab("age at stool collection (days)")+
    ylab("Granulicatella relative abundance %")+
    labs(color="age at first consumption of\nnon-breastmilk fluid/food (days)")
b <- ggplot(aes(x=Granulicatella.detected, y=complementary.feeding.age.days), data=age.v3) + 
    geom_boxplot() +
    geom_jitter() +
    theme_bw()+
    xlab("Granulicatella relative abundance > 0.00 %")+
    ylab("age at first consumption of\nnon-breastmilk fluid/food (days)")
ggarrange(a,b, ncol=2, nrow=1, widths=c(2,1))
ggsave("scatterplot.complementary.feeding.vs.Granulicatella.boxplot.Granulicatella.presence.vs.complementary.feeding.tiff", device="tiff", plot=last_plot(), dpi=600, height=6,width=11, units="in")

