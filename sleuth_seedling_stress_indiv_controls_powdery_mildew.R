# Sleuth for differential expression analysis of RNA-seq quantified by Kallisto
# Philippa Borrill
# 20.9.2015 EDITED 07.07.2016

# Install sleuth
source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5")
install.packages("devtools")
install.packages("stringi")
devtools::install_github("pachterlab/sleuth")
install.packages("scales")

# Run sleuth
library("sleuth")
setwd("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis")

# set base directory 
base_dir <- "Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis"

# set directory of files to work with
samples_to_use <- "SRP041017_powdery_mildew"

# get sample IDs (will need to change "test_sleuth_septoria" when using a new folder)
sample_id <- dir(file.path(base_dir,samples_to_use))
sample_id

kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, samples_to_use, id))
kal_dirs

# swap the / for \\ to have windows paths
kal_dirs<- gsub("/","\\\\", kal_dirs)
kal_dirs

# load sample info (should have 2 columns: sample and condition)
s2c <- read.table ("Y:\\expression_browser\\TGAC_assembly\\analysis\\seedling_stress_analysis\\SRP041017_powdery_mildew.txt",sep="\t", header = TRUE, stringsAsFactors=FALSE)
s2c

# load kallisto processed data into the sleuth object (so)
so <- sleuth_prep(kal_dirs, s2c, ~ condition,target_mapping = NULL, max_bootstrap = 30)

# estimate parameters for sleuth response error measurement model
so <- sleuth_fit(so)

#see which conditions you have
models(so)

# perform differential analysis for each condition vs none
so <- sleuth_test(so, which_beta = 'conditionPowdery mildew pathogen E09 24 hours')

# look at results in Shiny
#sleuth_live(so)

# save results tables

results_table <- sleuth_results(so, 'conditionPowdery mildew pathogen E09 24 hours') 
head(results_table)
write.csv(results_table, "Powdery mildew pathogen E09 24 hours_sorted_abundances.csv")

edited_results_table <- results_table[results_table$qval<0.05,]
edited_results_table <- (edited_results_table[complete.cases(edited_results_table),])
dim(edited_results_table)
write.csv(edited_results_table, "Powdery mildew pathogen E09 24 hours_sorted_abundances_q_value_under_0.05.csv")


# perform differential analysis for each condition vs none
so <- sleuth_test(so, which_beta = 'conditionPowdery mildew pathogen E09 48 hours')

# save results tables

results_table <- sleuth_results(so, 'conditionPowdery mildew pathogen E09 48 hours') 
head(results_table)
write.csv(results_table, "Powdery mildew pathogen E09 48 hours_sorted_abundances.csv")

edited_results_table <- results_table[results_table$qval<0.05,]
edited_results_table <- (edited_results_table[complete.cases(edited_results_table),])
dim(edited_results_table)
write.csv(edited_results_table, "Powdery mildew pathogen E09 48 hours_sorted_abundances_q_value_under_0.05.csv")


# perform differential analysis for each condition vs none
so <- sleuth_test(so, which_beta = 'conditionPowdery mildew pathogen E09 72 hours')

# save results tables

results_table <- sleuth_results(so, 'conditionPowdery mildew pathogen E09 72 hours') 
head(results_table)
write.csv(results_table, "Powdery mildew pathogen E09 72 hours_sorted_abundances.csv")

edited_results_table <- results_table[results_table$qval<0.05,]
edited_results_table <- (edited_results_table[complete.cases(edited_results_table),])
dim(edited_results_table)
write.csv(edited_results_table, "Powdery mildew pathogen E09 72 hours_sorted_abundances_q_value_under_0.05.csv")

# to save so

saveRDS(so, file = 'condition SRP041017_powdery_mildew.rds')

# to reload  the file later:
  
so <- readRDS('shiny/sleuth/so.rds')
sleuth_live(so)



#### meanings of columns from results_table - all except pval and qval are LOG VALUES (LOG, NOT LOG2)
#- target_id transcript name
#- pval - pvalue
#- qval - FDR adjusted pvalue using benjamini-hochberg
#- b - the 'beta' value (analogous to fold change, though technically a bias estimator which has to do with the transformation)
#- se_b - the standard error of beta
#- mean_obs - the mean of the observations. this is used for the smoothing
#- var_obs - the variance of the observations
#- tech_var - the technical variance from the bootstraps
#- sigma_sq - the raw estimator of the variance once the tech_var has been removed
#- smooth_sigma_sq - the smooth regression fit for the shrinkage estimation
#- final_sigma_sq - max(sigma_sq, smooth_sigma_sq). this is the one used for covariance estimation of beta (in addition to tech_var)


