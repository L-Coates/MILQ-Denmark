#This script contains the code used to identify microbial correlations in the 
#MILQ Danish stool samples. 

##load needed libraries
library(qiime2R)
library(tidyverse)
library(phyloseq)
library(vegan)
library(NetCoMi)
library(LaplacesDemon)

#STEP 1: read in the metadata, (unrarefied) feature table, and the taxonomy table
#and make into phyloseq object
metadata <- read.csv("../morbidity_medication_prevalences/infant.morbidities.medications.csv", header = T)
metadata <- metadata[,-1]

feature.table <- read_qza("../feature-table-no-chloroplast-eukarya-mitochondria.qza")
feature.table <- feature.table$data

taxonomy <- read_qza("../taxonomy-classification.qza")
taxonomy2 <- as.data.frame(taxonomy$data) %>% column_to_rownames("Feature.ID") %>% separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>% as.matrix()
#remove "Confidence" column
taxonomy2 <- taxonomy2[,-8]

#trim taxonomy to include only ASVs in the ASV table
taxonomy3 <- taxonomy2[c(rownames(taxonomy2)%in%rownames(feature.table)),]

#make phyloseq object
phylobj <- phyloseq(otu_table(feature.table, taxa_are_rows=T),tax_table(taxonomy3), sample_data(metadata %>% as.data.frame() %>% column_to_rownames("sample.id")))
#get rid of the taxa that are not present at all
phylobj.filtered <- filter_taxa(phylobj, function(x) sum(x)>0, TRUE)
#look at the number of ASVs among the feature table
dim(otu_table(phylobj.filtered))
#there are 1782 ASVs
dim(sample_data(phylobj.filtered)) #327 samples. 

#STEP 2: Gather counts at genus level. 
phylobj.filtered.genus <- tax_glom(phylobj.filtered, 'Genus', NArm=F)
dim(otu_table(phylobj.filtered.genus)) #There are 244 genera and 327 samples

#Create a list of unique genus names that can be assigned to the features within the phyloseq object
taxa.names <- tax_table(phylobj.filtered.genus)[,6]

#replace empty names with family level identities
taxa.names2 <- tax_table(phylobj.filtered.genus)[,5] #make a list of family level designations
taxa.names[c(which(is.na(taxa.names)), grep("g__$", taxa.names))] <-taxa.names2[c(which(is.na(taxa.names)),grep("g__$", taxa.names))]

#replace empty names with order level identities
taxa.names3 <-tax_table(phylobj.filtered.genus)[,4] #make a list of order level designations
taxa.names[c(which(is.na(taxa.names2)), grep( "f__$", taxa.names))] <- taxa.names3[c(which(is.na(taxa.names2)), grep("f__$", taxa.names))]

#replace empty names with class level identities
taxa.names4 <- tax_table(phylobj.filtered.genus)[,3] #make a list of class level designations
taxa.names[c(which(is.na(taxa.names3)), grep("o__$", taxa.names))] <-taxa.names4[c(which(is.na(taxa.names3)),grep("o__$", taxa.names))]

#replace empty names with phylum level identities
taxa.names5 <- tax_table(phylobj.filtered.genus)[,2] #make a list of phylum level designations
taxa.names[c(which(is.na(taxa.names4)),grep("c__$", taxa.names))] <-taxa.names5[c(which(is.na(taxa.names4)),grep("c__$", taxa.names))]

#replace empty names with kingdom level identities
taxa.names6 <- tax_table(phylobj.filtered.genus)[,1] #make a list of kingdom level designations
taxa.names[c(which(is.na(taxa.names5)),grep("p__$", taxa.names))] <-taxa.names6[c(which(is.na(taxa.names5)),grep("p__$", taxa.names))]

#remove leading white space
taxa.names <- gsub(" ","", taxa.names)

#change the taxa designations to names compatible with R
taxa.names <- make.names(taxa.names, unique = TRUE)

#Make all of the feature names within the phyloseq object correspond to the genus level designation
taxa_names(phylobj.filtered.genus) <- taxa.names

#important to note that when extracting the sample data (metadata) from the phyloseq
#object, need to use the function data.frame(), not as.data.frame(), to change it
#into a dataframe object that maaslin2 will recognize, otherwise it will still be 
#recognized as a sample_data object if you use as.data.frame(). 
feature.table.genus <- as.data.frame(t(otu_table(phylobj.filtered.genus)))

#remove the genera that have an average of less than 2 reads. 
genus.count.average <- data.frame(apply(X=feature.table.genus,MARGIN=2, FUN=mean))
genus.count.average$genus <- rownames(genus.count.average)
genus.count.average.keep <- genus.count.average[genus.count.average$apply.X...feature.table.genus..MARGIN...2..FUN...mean.>=2,]

#filter the feature table
feature.table.genus.filtered <- feature.table.genus[,c(which(colnames(feature.table.genus) %in% genus.count.average.keep$genus)),]

#STEP 3: Determining microbial correlations with "spieceasi" measure

colnames(feature.table.genus.filtered) <- gsub(pattern="^g__|^f__",replacement="",colnames(feature.table.genus.filtered))

spieceasi.network <- netConstruct(data=feature.table.genus.filtered, dataType="counts", measure="spieceasi", adjust="lfdr",seed=12)
spiec.easi.network.analyzed <- netAnalyze(spieceasi.network)

jpeg("microbial.correlations.network.jpeg", res=600, width=20, height=20, units="cm")    
plot(spiec.easi.network.analyzed, 
     nodeColor = "cluster",
     repulsion = 0.8,
     rmSingles = TRUE,
     labelScale = FALSE,
     cexLabels = 0.8,
     nodeSize="normCounts",
     nodeSizeSpread = 5,
     cexNodes = 2,
     title1 = "Network on genus level correlations with SPIEC-EASI", 
     showTitle = TRUE,
     hubBorderCol = "black",
     cexTitle = 1)
legend(0.7, 0.95, cex = 0.8, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"),
       bty = "n", horiz = TRUE)
dev.off()


#STEP 4: import the rarefied feature table and run multiple spearman correlations
#for Granulicatella abundance in relation to other microbes. 

feature.table.rarefied <- read_qza("../feature-table-rarefied-6963.qza")
feature.table.rarefied <- feature.table.rarefied$data

phylobj.rarefied <- phyloseq(otu_table(feature.table.rarefied, taxa_are_rows=T),tax_table(taxonomy3), sample_data(metadata %>% as.data.frame() %>% column_to_rownames("sample.id")))

#gather the counts at the genus level. 
phylobj.rarefied.filtered <- filter_taxa(phylobj.rarefied, function(x) sum(x)>0, TRUE)
#look at the number of ASVs among the rarefied feature table
dim(otu_table(phylobj.rarefied.filtered))
#there are 1585 ASVs
dim(sample_data(phylobj.rarefied.filtered)) #327 samples. 

#Gather counts at genus level. 
phylobj.rarefied.filtered.genus <- tax_glom(phylobj.rarefied.filtered, 'Genus', NArm=F)
dim(otu_table(phylobj.rarefied.filtered.genus)) #There are 222 genera and 327 samples

#Create a list of unique genus names that can be assigned to the features within the phyloseq object
taxa.names <- tax_table(phylobj.rarefied.filtered.genus)[,6]

#replace empty names with family level identities
taxa.names2 <- tax_table(phylobj.rarefied.filtered.genus)[,5] #make a list of family level designations
taxa.names[c(which(is.na(taxa.names)), grep("g__$", taxa.names))] <-taxa.names2[c(which(is.na(taxa.names)),grep("g__$", taxa.names))]

#replace empty names with order level identities
taxa.names3 <-tax_table(phylobj.rarefied.filtered.genus)[,4] #make a list of order level designations
taxa.names[c(which(is.na(taxa.names2)), grep( "f__$", taxa.names))] <- taxa.names3[c(which(is.na(taxa.names2)), grep("f__$", taxa.names))]

#replace empty names with class level identities
taxa.names4 <- tax_table(phylobj.rarefied.filtered.genus)[,3] #make a list of class level designations
taxa.names[c(which(is.na(taxa.names3)), grep("o__$", taxa.names))] <-taxa.names4[c(which(is.na(taxa.names3)),grep("o__$", taxa.names))]

#replace empty names with phylum level identities
taxa.names5 <- tax_table(phylobj.rarefied.filtered.genus)[,2] #make a list of phylum level designations
taxa.names[c(which(is.na(taxa.names4)),grep("c__$", taxa.names))] <-taxa.names5[c(which(is.na(taxa.names4)),grep("c__$", taxa.names))]

#replace empty names with kingdom level identities
taxa.names6 <- tax_table(phylobj.rarefied.filtered.genus)[,1] #make a list of kingdom level designations
taxa.names[c(which(is.na(taxa.names5)),grep("p__$", taxa.names))] <-taxa.names6[c(which(is.na(taxa.names5)),grep("p__$", taxa.names))]

#remove leading white space
taxa.names <- gsub(" ","", taxa.names)

#change the taxa designations to names compatible with R
taxa.names <- make.names(taxa.names, unique = TRUE)

#Make all of the feature names within the phyloseq object correspond to the genus level designation
taxa_names(phylobj.rarefied.filtered.genus) <- taxa.names

melted.phylgenus.rarefied <- psmelt(phylobj.rarefied.filtered.genus)
length(unique(melted.phylgenus.rarefied$OTU))

#run a function that tests for a correlation between Granulicatella and all the 
#other 221 genera. 

melted.phylgenus.rarefied.wide <- pivot_wider(data=melted.phylgenus.rarefied, id_cols=Sample, names_from=OTU, values_from=Abundance)

correlation.func <- function(genus){
    correlations <- data.frame()
    for(i in genus){
        x.vector <- melted.phylgenus.rarefied.wide[["g__Granulicatella"]]
        y.vector <- melted.phylgenus.rarefied.wide[[i]]
        corr.output <-cor.test(x=x.vector,y=y.vector, method="spearman")
        corr.output.2 <- data.frame(genus_1="g__Granulicatella", genus_2=i, estimate=corr.output[["estimate"]], pvalue=corr.output[["p.value"]])
        correlations <- rbind(correlations, corr.output.2)   
    }
    return(correlations)
}

genus<-colnames(melted.phylgenus.rarefied.wide[,2:223])
spearman.correlations <- correlation.func(genus)
spearman.correlations <- spearman.correlations[spearman.correlations$genus_2!="g__Granulicatella",]

spearman.correlations$pvalue.fdr.adjusted <- p.adjust(p=spearman.correlations$pvalue, method="fdr")
#many, except one, of the significant correlations were positive. 

#visualize some the first three correlations (by adjusted p-value)
ggplot(aes(x=`g__Granulicatella`, y=`g__Gemella`), data=melted.phylgenus.rarefied.wide)+geom_point()
ggplot(aes(x=`g__Granulicatella`, y=`g__Haemophilus`), data=melted.phylgenus.rarefied.wide)+geom_point()
ggplot(aes(x=`g__Granulicatella`, y=`g__Streptococcus`), data=melted.phylgenus.rarefied.wide)+geom_point()
