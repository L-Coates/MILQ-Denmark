# This is the data repository for QIIME2 commands and R code used to analyze microbial data from Danish infant stool in relation to diarrhea, fever, and vomiting in the Mothers, Infants, and Lactation Quality study. 

## R code files:
### "alpha-diversity.R"
This file contains R code used to determine if infant stool alpha-diversity differed with infant morbidity from the same visit. Linear mixed effects modeling was employed to determine if diarrhea, fever or vomiting were significant predictors of alpha-diversity. The various parameters and covariates included in models are shown in the code file. 

### "batch.effect.R"
This file contains R code used (1) to assess the variation in infant stool [DNA] and read depth by extraction batch or bead basher type (old vs. new) via one-way ANOVA; (2) to determine read depth varied by DNA extraction batch, bead basher, plate batch, or sequencing batch via student's t-test or wilcoxon rank sum test; (3) to determine if beta-diversity varied by DNA extraction batch, PCR/plate batch, sequencing batch, read depth, or [DNA] using permutational multivariate analysis of variance.

### "beta-diversity.R"
This file contains R code used to determine if beta-diversity differed by infant morbidity. Permutational multivariate analysis of variances was employed to determine whether weighted or unweighted UniFrac distances associated with infant diarrhea, fever, or vomiting. We investigated associations within each visit (cross-sectional approach) and across all visits (longitudinal approach). Co-variates were included in the full models. Reduced models were also generated that did excluded co-variates which were not found to be significantly associated with beta-diversity in the full model. For the morbidities and co-variates that were significantly associated with beta-diversity, principal coordinate analysis was performed to visualize clustering of samples by the variable of interest.

### "complementary_feeding_Granulicatella.R"
This script contains the code used to determine infant age at first introduction to complementary food/fluids, and code used to investigate associations between infant age at first complementary feeding and infant stool Granulicatella abundance. 

### "core_microbes.R"
This script contains the code used to identify core microbial taxa (ASVs, species, genera, families) present in all/most samples and to determine the taxa (genera and families) that were particularly abundant among this cohort. 

### "Denmark_morbidity_frequency.R"
This file contains the R code used to determine the frequency of diarrhea, fever, vomiting, and medications in the infant Denmark cohort. 

### "Denmark_stool_collection_info.R"
This file contains the code used to assess the metadata on the Danish infant stool collection and storage, including age at stool collection, and time between subsequent stool collections.

### "differential-abundance.R"
This file contains the R code used to test for differential abundance of stool microbes with diarrhea, fever, or vomiting. Analysis of Compositions of Microbiomes with Bias Correction 2 (ANCOM-BC2) is a stastical method for conducting differential abundance analysis of microbial counts and can be applied to longitudinal data. ANCOM-BC2 was used to determine differential abundance of stool microbes at different taxonomic levels with diarrhea, fever, or vomiting. The various parameters and covariates included in models are shown in the code file. For the microbial taxa found differentially abundant with diarhea, fever, vomiting, or a covariate, wilcoxon rank sum test was performed on the rarefied counts to confirm differential abundance with the particular condition. Graphs were also created to display statistically significant relationships between microbial abundance and diarrhea, fever, vomiting, or covariates.

### "exploring_taxa_from_taxaHFE_output.R"
This file contains the code used to explore the microbial taxa from visit 2 that were identified with "taxaHFE 2.0" as potentially associated with diarrhea, fever, or vomitting in visits 3 and/or 4. Taxa identified as important in predicting infant morbidity were included in random forest models. Model performance was determined from a confusion matrix and ROC AUC. Feature importance and shapley values was determined to assess the contribution of variables to the predictive power of the model. Taxa identified as important by taxaHFE were also investigated for differential abundance by infant morbidity (using wilcoxon rank sum test).  

### "network_analysis_NetCoMi.R"
This script contains the code used to identify microbial genera correlations in the infant stool, across all visits, using NetCoMi. 

### "predicting_morbidity_outcomes_from_alpha_diversity.R"
This script contains the code used to determine if infant stool alpha diversity in visit 2 is predictive of infant diarrhea, fever, or vomiting in visits 3 and/or 4. Logistic regression models were generated and the performance of the models were assessed with confusion matrices and ROC AUC.  

### "qiime2-scripts.rtf"
This file contains a collation of qiime2 scripts executed on a high performance computing cluster. 

### "taxaHFE.commands.txt"
These are the commands used with taxaHFE version 2.0 to identify taxa from visit 2 stool samples potentially associated with diarrhea, fever, or vomiting in visit 3 and/or visit 4.
 
