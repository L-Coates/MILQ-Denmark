# This is the data repository for QIIME2 commands and R code used to analyze microbial data from Danish infant stool in relation to diarrhea in the Mother infant Lactation Quality project.  

## The QIIME2 commands file "qiime2.commands.txt"

## The R code file "differential-abundance.R"
This file contains the R code used to test for differential abundance of stool microbes with diarrhea. Analysis of Compositions of Microbiomes with Bias Correction 2 (ANCOM-BC2) is a stastical method for conducting differential abundance analysis of microbial counts and can be applied to longitudinal data. ANCOM-BC2 was used to determine differential abundance of stool microbes at different taxonomic levels with diarrhea. The various parameters and covariates included in models are shown in the code file. For the microbial taxa found differentially abundant with diarhea or a covariate, wilcoxon rank sum test was performed on the rarefied counts to confirm differential abundance with the particular condition. Graphs were also created to display statistically significant relationships between microbial abundance and diarrhea or covariates.   
 