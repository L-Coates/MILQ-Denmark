# This is the data repository for QIIME2 commands and R code used to analyze microbial data from Danish infant stool in relation to diarrhea, fever, and vomiting in the Mother Infant Lactation Quality project. 

## R code files:
### "alpha-diversity.R"

### "batch.effect.R"

### "beta-diversity.R"

### "complementary_feeding_Granulicatella.R"

### "core_microbes.R"

### "demographics_birth_mode_etc.R"

### "Denmark_morbidity_frequency.R"

### "Denmark_stool_collection_info.R"

### "differential-abundance.R"
This file contains the R code used to test for differential abundance of stool microbes with diarrhea, fever, or vomiting. Analysis of Compositions of Microbiomes with Bias Correction 2 (ANCOM-BC2) is a stastical method for conducting differential abundance analysis of microbial counts and can be applied to longitudinal data. ANCOM-BC2 was used to determine differential abundance of stool microbes at different taxonomic levels with diarrhea, fever, or vomiting. The various parameters and covariates included in models are shown in the code file. For the microbial taxa found differentially abundant with diarhea, fever, vomiting, or a covariate, wilcoxon rank sum test was performed on the rarefied counts to confirm differential abundance with the particular condition. Graphs were also created to display statistically significant relationships between microbial abundance and diarrhea, fever, vomiting, or covariates.

### "exploring_taxa_from_taxaHFE_output.R"

### "formatting_hierarchichal_data.R"

### "network_analysis_NetCoMi.R"

### "predicting_morbidity_outcomes_from_alpha_diversity.R"

### "relating_fever_to_vaccinations.R"

### "taxaHFE.commands.txt"
 
