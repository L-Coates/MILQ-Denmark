#This script contains the code used to identify core microbial taxa present
#all/most samples and to determine the taxa that were particularly abundant
#among this cohort. 

library(qiime2R)
library(phyloseq)
library(tidyverse)


#STEP 1: read in the qiime2 feature table and taxonomy and combine with metadata
#to build phyloseq object
metadata <- read.csv("../metadata-infants-complete-stool-set-after-6886-read-cutoff-withoutInfantCD144Y.csv") 
metadata <- metadata[,-1]

feature.table <- read_qza("../feature-table-no-chloroplast-eukarya-mitochondria.qza")
feature.table <- feature.table$data

taxonomy <- read_qza("../taxonomy-classification.qza")
taxonomy2 <- as.data.frame(taxonomy$data) %>% column_to_rownames("Feature.ID") %>% separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>% as.matrix()
#remove "Confidence" column
taxonomy2 <- taxonomy2[,-8]

#trim taxonomy to include only ASVs in the ASV table
taxonomy3 <- taxonomy2[c(rownames(taxonomy2)%in%rownames(feature.table)),]

phylobj <- phyloseq(otu_table(feature.table, taxa_are_rows=T),tax_table(taxonomy3), sample_data(metadata %>% as.data.frame() %>% column_to_rownames("sample.id")))
#get rid of the taxa that are not present at all
phylobj.filtered <- filter_taxa(phylobj, function(x) sum(x)>0, TRUE)
#look at the number of ASVs among the feature table
dim(otu_table(phylobj.filtered))
#there are 1782 ASVs
dim(sample_data(phylobj.filtered)) #327 samples. 


#STEP 2: Summarize feature counts at ASV, species, genus, and family levels. 

ASV.counts <- psmelt(phylobj.filtered)

#summarizing counts at species level now
phylobj.species <- tax_glom(phylobj.filtered, 'Species', NArm=F)

#Create a list of unique taxa names that can be assigned to the features within the phyloseq object
taxa.names <- tax_table(phylobj.species)[,7]

#replace empty names with genus level identities
taxa.names1 <- tax_table(phylobj.species)[,6] 
taxa.names[c(which(is.na(taxa.names)), grep("s__$", taxa.names))] <-taxa.names1[c(which(is.na(taxa.names)),grep("s__$", taxa.names))]

#replace empty names with family level identities
taxa.names2 <- tax_table(phylobj.species)[,5] 
taxa.names[c(which(is.na(taxa.names1)), grep("s__$", taxa.names))] <-taxa.names2[c(which(is.na(taxa.names1)),grep("s__$", taxa.names))]

#replace empty names with order level identities
taxa.names3 <-tax_table(phylobj.species)[,4] 
taxa.names[c(which(is.na(taxa.names2)), grep( "s__$", taxa.names))] <- taxa.names3[c(which(is.na(taxa.names2)), grep("s__$", taxa.names))]

#replace empty names with class level identities
taxa.names4 <- tax_table(phylobj.species)[,3]
taxa.names[c(which(is.na(taxa.names3)), grep("s__$", taxa.names))] <-taxa.names4[c(which(is.na(taxa.names3)),grep("s__$", taxa.names))]

#replace empty names with phylum level identities
taxa.names5 <- tax_table(phylobj.species)[,2] 
taxa.names[c(which(is.na(taxa.names4)),grep("s__$", taxa.names))] <-taxa.names5[c(which(is.na(taxa.names4)),grep("s__$", taxa.names))]

#replace empty names with kingdom level identities
taxa.names6 <- tax_table(phylobj.species)[,1] 
taxa.names[c(which(is.na(taxa.names5)),grep("s__$", taxa.names))] <-taxa.names6[c(which(is.na(taxa.names5)),grep("s__$", taxa.names))]

#remove leading white space
taxa.names <- gsub(" ","", taxa.names)

#change the taxa designations to names compatible with R
taxa.names <- make.names(taxa.names, unique = TRUE)

#Make all of the feature names within the phyloseq object correspond to the species level designation
taxa_names(phylobj.species) <- taxa.names

species.counts <- psmelt(phylobj.species)

#gather counts at the genus level now. 
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

genus.counts <- psmelt(phylobj.genus)

#gathering counts at the family level
phylobj.family <- tax_glom(phylobj.filtered, 'Family', NArm=F)

#Create a list of unique taxa names that can be assigned to the features within the phyloseq object
taxa.names <- tax_table(phylobj.family)[,5]

#replace empty names with order level identities
taxa.names2 <- tax_table(phylobj.family)[,4] 
taxa.names[c(which(is.na(taxa.names)), grep("f__$", taxa.names))] <-taxa.names2[c(which(is.na(taxa.names)),grep("f__$", taxa.names))]

#replace empty names with class level identities
taxa.names3 <-tax_table(phylobj.family)[,3] 
taxa.names[c(which(is.na(taxa.names2)), grep( "f__$", taxa.names))] <- taxa.names3[c(which(is.na(taxa.names2)), grep("f__$", taxa.names))]

#replace empty names with phylum level identities
taxa.names4 <- tax_table(phylobj.family)[,2]
taxa.names[c(which(is.na(taxa.names3)), grep("f__$", taxa.names))] <-taxa.names4[c(which(is.na(taxa.names3)),grep("f__$", taxa.names))]

#replace empty names with kingdom level identities
taxa.names5 <- tax_table(phylobj.family)[,1] 
taxa.names[c(which(is.na(taxa.names4)),grep("f__$", taxa.names))] <-taxa.names5[c(which(is.na(taxa.names4)),grep("f__$", taxa.names))]

#remove leading white space
taxa.names <- gsub(" ","", taxa.names)

#change the taxa designations to names compatible with R
taxa.names <- make.names(taxa.names, unique = TRUE)

#Make all of the feature names within the phyloseq object correspond to the family level designation
taxa_names(phylobj.family) <- taxa.names

family.counts <- psmelt(phylobj.family)


#STEP 3: Determine ASVs, species, genera, and families present in all 327 samples. 
ASV.prevalence <- ASV.counts[ASV.counts$Abundance>0,]
ASV.prevalence <- data.frame(table(ASV.prevalence$OTU))

species.prevalence <- species.counts[species.counts$Abundance>0,]
species.prevalence <- data.frame(table(species.prevalence$OTU))

genus.prevalence <- genus.counts[genus.counts$Abundance>0,]
genus.prevalence <- data.frame(table(genus.prevalence$OTU))

family.prevalence <- family.counts[family.counts$Abundance>0,]
family.prevalence <- data.frame(table(family.prevalence$OTU))

#STEP 4: read in the rarefied counts and determine some of the top most abundant taxa by 
#average relative abundance across the entire dataset. 

rarefied.feature.table <- read_qza("../feature-table-rarefied-6963.qza")
rarefied.feature.table <- rarefied.feature.table$data

#trim taxonomy to include only ASVs in the ASV table
taxonomy3 <- taxonomy2[c(rownames(taxonomy2)%in%rownames(rarefied.feature.table)),]

phylobj.rarefied <- phyloseq(otu_table(rarefied.feature.table, taxa_are_rows=T),tax_table(taxonomy3), sample_data(metadata %>% as.data.frame() %>% column_to_rownames("sample.id")))
#get rid of the taxa that are not present at all
phylobj.rarefied.filtered <- filter_taxa(phylobj.rarefied, function(x) sum(x)>0, TRUE)
#look at the number of ASVs among the feature table
dim(otu_table(phylobj.rarefied.filtered))
#there are 1585 ASVs
dim(sample_data(phylobj.rarefied.filtered)) #327 samples. 

#gather counts at the genus level now. 
phylobj.genus.rarefied <- tax_glom(phylobj.rarefied.filtered, 'Genus', NArm=F)

#Create a list of unique taxa names that can be assigned to the features within the phyloseq object
taxa.names <- tax_table(phylobj.genus.rarefied)[,6]

#replace empty names with family level identities
taxa.names2 <- tax_table(phylobj.genus.rarefied)[,5] 
taxa.names[c(which(is.na(taxa.names)), grep("g__$", taxa.names))] <-taxa.names2[c(which(is.na(taxa.names)),grep("g__$", taxa.names))]

#replace empty names with order level identities
taxa.names3 <-tax_table(phylobj.genus.rarefied)[,4] 
taxa.names[c(which(is.na(taxa.names2)), grep( "g__$", taxa.names))] <- taxa.names3[c(which(is.na(taxa.names2)), grep("g__$", taxa.names))]

#replace empty names with class level identities
taxa.names4 <- tax_table(phylobj.genus.rarefied)[,3]
taxa.names[c(which(is.na(taxa.names3)), grep("g__$", taxa.names))] <-taxa.names4[c(which(is.na(taxa.names3)),grep("g__$", taxa.names))]

#replace empty names with phylum level identities
taxa.names5 <- tax_table(phylobj.genus.rarefied)[,2] 
taxa.names[c(which(is.na(taxa.names4)),grep("g__$", taxa.names))] <-taxa.names5[c(which(is.na(taxa.names4)),grep("g__$", taxa.names))]

#replace empty names with kingdom level identities
taxa.names6 <- tax_table(phylobj.genus.rarefied)[,1] 
taxa.names[c(which(is.na(taxa.names5)),grep("g__$", taxa.names))] <-taxa.names6[c(which(is.na(taxa.names5)),grep("g__$", taxa.names))]

#remove leading white space
taxa.names <- gsub(" ","", taxa.names)

#change the taxa designations to names compatible with R
taxa.names <- make.names(taxa.names, unique = TRUE)

#Make all of the feature names within the phyloseq object correspond to the genus level designation
taxa_names(phylobj.genus.rarefied) <- taxa.names

genus.counts.rarefied <- psmelt(phylobj.genus.rarefied)

genus.counts.rarefied.RelAbund <- data.frame(tapply(X=genus.counts.rarefied$Abundance, INDEX=genus.counts.rarefied$Genus, FUN=median))

#gathering counts at the family level
phylobj.family.rarefied <- tax_glom(phylobj.rarefied.filtered, 'Family', NArm=F)

#Create a list of unique taxa names that can be assigned to the features within the phyloseq object
taxa.names <- tax_table(phylobj.family.rarefied)[,5]

#replace empty names with order level identities
taxa.names2 <- tax_table(phylobj.family.rarefied)[,4] 
taxa.names[c(which(is.na(taxa.names)), grep("f__$", taxa.names))] <-taxa.names2[c(which(is.na(taxa.names)),grep("f__$", taxa.names))]

#replace empty names with class level identities
taxa.names3 <-tax_table(phylobj.family.rarefied)[,3] 
taxa.names[c(which(is.na(taxa.names2)), grep( "f__$", taxa.names))] <- taxa.names3[c(which(is.na(taxa.names2)), grep("f__$", taxa.names))]

#replace empty names with phylum level identities
taxa.names4 <- tax_table(phylobj.family.rarefied)[,2]
taxa.names[c(which(is.na(taxa.names3)), grep("f__$", taxa.names))] <-taxa.names4[c(which(is.na(taxa.names3)),grep("f__$", taxa.names))]

#replace empty names with kingdom level identities
taxa.names5 <- tax_table(phylobj.family.rarefied)[,1] 
taxa.names[c(which(is.na(taxa.names4)),grep("f__$", taxa.names))] <-taxa.names5[c(which(is.na(taxa.names4)),grep("f__$", taxa.names))]

#remove leading white space
taxa.names <- gsub(" ","", taxa.names)

#change the taxa designations to names compatible with R
taxa.names <- make.names(taxa.names, unique = TRUE)

#Make all of the feature names within the phyloseq object correspond to the family level designation
taxa_names(phylobj.family.rarefied) <- taxa.names

family.counts.rarefied <- psmelt(phylobj.family.rarefied)

family.counts.rarefied.RelAbund <- data.frame(tapply(X=family.counts.rarefied$Abundance, INDEX=family.counts.rarefied$Family, FUN=median))
