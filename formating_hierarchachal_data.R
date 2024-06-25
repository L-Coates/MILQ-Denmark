#This file contains the code used to format the data and the metadata for taxaHFE use.

#loading libraries
library(qiime2R)
library(tidyverse)
library(phyloseq)
library(tidyverse)

#Step 1: reading in the rarefied feature file, taxonomy, and metadata file and
#creating a phyloseq object (for all 327 samples).
feature.table.rarefied <- read_qza("../feature-table-rarefied-6963.qza")
feature.table.rarefied <- feature.table.rarefied$data

taxonomy <- read_qza("../taxonomy-classification.qza")
taxonomy2 <- as.data.frame(taxonomy$data) %>% column_to_rownames("Feature.ID") %>% separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>% as.matrix()
#remove "Confidence" column
taxonomy2 <- taxonomy2[,-8]

#trim taxonomy to include only ASVs in the ASV table
taxonomy3 <- taxonomy2[c(rownames(taxonomy2)%in%rownames(feature.table.rarefied)),]

metadata <- read.csv("../morbidity_medication_prevalences/infant.morbidities.medications.csv")
metadata <- metadata[,-1]

phylobj <- phyloseq(otu_table(feature.table.rarefied, taxa_are_rows=T),tax_table(taxonomy3), sample_data(metadata %>% as.data.frame() %>% column_to_rownames("sample.id")))

#remove any taxa that aren't present in any samples. 
phylobj <- filter_taxa(phylobj, function(x) sum(x)>0, TRUE)

#Step 2: gather counts at species level and transform to relative abundance
phylobj.species <- tax_glom(phylobj, "Species", NArm=F)
phylobj.species.relabund <- transform_sample_counts(phylobj.species, function(x) x/sum(x) * 100)

feature.table.rarefied.species <- as.data.frame(otu_table(phylobj.species.relabund))
taxonomy.species <- as.data.frame(tax_table(phylobj.species.relabund))
#looking to see if any species names are non-unique
table(taxonomy.species$Species)[table(taxonomy.species$Species)>1]
#there are quite a few species (e.g. s__uncultured_bacterium) without unique names,
#so going to concatenate the higher level taxonomic names to the species name. taxaHFE output only prints the species name, not the full clade name
#so for these species that have the same name, we will be able to determine
#taxonomy if we append the higher level taxonomic names.  
taxonomy.species$Species <- apply(X=taxonomy.species,MARGIN=1, function(x) paste(x, collapse="_"))
#confirm that all species names are unique now that they have been combined with 
#genus names. 
table(taxonomy.species$Species)[table(taxonomy.species$Species)>1]

#for the taxa that had "NA" for species, the entry will need to be changed back
#to "NA". 
taxonomy.species[c(grep(pattern="_NA$", taxonomy.species$Species)), "Species"] <- NA

taxonomy.species$clade_name <- apply(X=taxonomy.species,MARGIN=1, function(x) paste(x, collapse="|"))
#get rid of the space before each taxonomic name 
taxonomy.species$clade_name <- gsub(pattern=" ", replacement="", taxonomy.species$clade_name)
unique(taxonomy.species$clade_name)
#get rid of "NA" taxonomic names
taxonomy.species$clade_name <- gsub(pattern="\\|NA", replacement="", taxonomy.species$clade_name)
unique(taxonomy.species$clade_name)

#Step 3: add the clade_name column to the feature table, then print out the feature table. 
feature.table.rarefied.species$ASV <- rownames(feature.table.rarefied.species) 
taxonomy.species$ASV <- rownames(taxonomy.species)

feature.table.rarefied.species.v2 <- merge(x=feature.table.rarefied.species, y=taxonomy.species[,c(8:9)], by="ASV", all=F)

#get rid of ASV column and place clade_name in first column of data frame. 
feature.table.rarefied.species.v3 <- feature.table.rarefied.species.v2[,c(329,2:328)]

write.csv(feature.table.rarefied.species.v3, row.names = F, "taxaHFE.feature.table.rarefied.counts.species.csv")

#write the metadata to this folder too. 
write.csv(metadata[,c(grep(pattern="^sample.id|infant.diarrhea.in.visit3.or.visit4|infant.fever.in.visit3.or.visit4|infant.vomit.in.visit3.or.visit4", colnames(metadata)))], "taxaHFE.metadata.csv", row.names = F)

#Step 4: now, running the same steps again but just for visit 2 samples (109 total)
phylobj.species.visit2 <- prune_samples(x=phylobj.species.relabund, samples=c(metadata[metadata$visit=="2","sample.id"]))
phylobj.species.visit2 <- filter_taxa(phylobj.species.visit2, function(x) sum(x)>0, TRUE)

feature.table.rarefied.species.visit2<- as.data.frame(otu_table(phylobj.species.visit2))
taxonomy.species.visit2 <- as.data.frame(tax_table(phylobj.species.visit2))
#looking to see if any species names are non-unique
table(taxonomy.species.visit2$Species)[table(taxonomy.species.visit2$Species)>1]
#there are quite a few species (e.g. s__uncultured_bacterium) without unique names,
#so going to concatenate the higher level taxonomic names to the species name. taxaHFE output only prints the species name, not the full clade name
#so for these species that have the same name, we will be able to determine
#taxonomy if we append the higher level taxonomic names.  
taxonomy.species.visit2$Species <- apply(X=taxonomy.species.visit2,MARGIN=1, function(x) paste(x, collapse="_"))
#confirm that all species names are unique now that they have been combined with 
#genus names. 
table(taxonomy.species.visit2$Species)[table(taxonomy.species.visit2$Species)>1]

#for the taxa that had "NA" for species, the entry will need to be changed back
#to "NA". 
taxonomy.species.visit2[c(grep(pattern="_NA$", taxonomy.species.visit2$Species)), "Species"] <- NA

taxonomy.species.visit2$clade_name <- apply(X=taxonomy.species.visit2,MARGIN=1, function(x) paste(x, collapse="|"))
#get rid of the space before each taxonomic name 
taxonomy.species.visit2$clade_name <- gsub(pattern=" ", replacement="", taxonomy.species.visit2$clade_name)
unique(taxonomy.species.visit2$clade_name)
#get rid of "NA" taxonomic names
taxonomy.species.visit2$clade_name <- gsub(pattern="\\|NA", replacement="", taxonomy.species.visit2$clade_name)
unique(taxonomy.species.visit2$clade_name)

#Step 5: add the clad_name column to the feature table, then print out the feature table. 
feature.table.rarefied.species.visit2$ASV <- rownames(feature.table.rarefied.species.visit2) 
taxonomy.species.visit2$ASV <- rownames(taxonomy.species.visit2)

feature.table.rarefied.species.visit2.v2 <- merge(x=feature.table.rarefied.species.visit2, y=taxonomy.species.visit2[,c(8:9)], by="ASV", all=F)

#get rid of ASV column and place clade_name in first column of data frame. 
feature.table.rarefied.species.visit2.v3 <- feature.table.rarefied.species.visit2.v2[,c(111,2:110)]

write.csv(feature.table.rarefied.species.visit2.v3, row.names = F, "taxaHFE.feature.table.rarefied.counts.species.visit2.csv")

#write the metadata to this folder too, but only the columns of interest (i.e. the 
#morbidity outcome, sample.id). 
metadata.visit2 <- data.frame(sample_data(phylobj.species.visit2))
metadata.visit2$sample.id <- rownames(metadata.visit2)
write.csv(metadata.visit2[,c(grep(pattern="^sample.id|infant.diarrhea.in.visit3.or.visit4|infant.fever.in.visit3.or.visit4|infant.vomit.in.visit3.or.visit4", colnames(metadata.visit2)))], "taxaHFE.metadata.visit2.csv", row.names = F)

