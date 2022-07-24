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

write.csv(meta_pheno,paste0(path,"ukbb_manifest_filtered_phenos.csv"),row.names = F)

dat <- read.csv(paste0(path,"ukbb_manifest_filtered_phenos.csv"))

# Adding the new filtered QC on the GWAS ----

dat_filtered <- read.csv(paste0(path,"ukbb_manifest_filtered_phenos.csv"))
dat_manifest <- read.csv(paste0(path,"Pan-UK Biobank phenotype manifest - phenotype_manifest.csv"))

dat_manifest <- dat_manifest %>% select(lambda_gc_EUR,phenotype_qc_EUR,aws_link_tabix,pops_pass_qc,sldsc_25bin_h2_observed_EUR)

dat_filtered <- merge(dat_filtered,dat_manifest,by="aws_link_tabix")

table(dat_filtered$phenotype_qc_EUR)

write.csv(dat_filtered,paste0(path,"ukbb_manifest_filtered_phenos.csv"),row.names = F)
