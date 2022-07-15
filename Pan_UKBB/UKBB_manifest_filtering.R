library(openxlsx)
library(tidyverse)

path <- './Pan_UKBB/'
meta_pheno <- readxl::read_xlsx(paste0(path,"Pan_UK_Biobank_phenotype_manifest.xlsx"))
# meta_pheno <- read_csv(paste0(path,"Pan_UK_Biobank_phenotype_manifest.csv"))

table(meta_pheno$trait_type)
meta_pheno_cleaned <- meta_pheno %>% #7221
                filter(saige_heritability_EUR>0.05) %>% #2725
                filter(pheno_sex=="both_sexes") %>% #2681
                filter(grepl("EUR",pops)) %>% #2681
                filter(!grepl("Country of Birth|Country of birth",description)) %>% #2653
                filter(!grepl("ethnic|Ethnic",description)) %>% #2647
                filter(!grepl("Ethnicity",category)) %>% #2646
                filter(!grepl("Male|male|female|Female",category)) %>% #2630
                filter(n_cases_full_cohort_females!=0,
                       n_cases_full_cohort_males!=0) %>% #2400
                filter(!grepl("Type milk consumed",description)) #2397

write.csv(meta_pheno_cleaned,paste0(path,"ukbb_manifest_filtered_phenos.csv"),row.names = F)

