import pandas as pd
import os
import numpy as np

ukbb_manifest = pd.read_csv("./Data/ukbb_manifest_filtered_phenos.csv")
##############################################################################
#############################   |            |   #############################                                    
#############################   | ANNOTATION |   #############################
#############################   |            |   #############################
##############################################################################

phenocode_annotate_lst = []
for rows in range(1087):
   print(rows)
   row = ukbb_manifest.iloc[rows,:]
   if row['trait_type'] =='phecode':
      phenocode_annotate = 'PH_'+ row.description.replace(' ','_') + '_h' + str("{:.5f}".format(row.saige_heritability_EUR)) + '_n' + row.n_cases_full_cohort_both_sexes.astype(str)
   if row['trait_type'] =='prescriptions': 
      phenocode_annotate = 'PR_'+ row.phenocode.replace(' ','_').replace('|','_') + '_h' + str("{:.5f}".format(row.saige_heritability_EUR)) + '_n' + row.n_cases_full_cohort_both_sexes.astype(str)
   if row['trait_type'] =='icd10':
      phenocode_annotate = 'IC_'+ row.category.split("|",)[-1][5:].replace(' ','_') + '_h' + str("{:.5f}".format(row.saige_heritability_EUR)) + '_n' + row.n_cases_full_cohort_both_sexes.astype(str)
   if row['trait_type'] =='categorical':
      phenocode_annotate = 'CA_'+ row.category.split(">",)[-1].replace(' ','_') + '_h' + str("{:.5f}".format(row.saige_heritability_EUR)) + '_n' + row.n_cases_full_cohort_both_sexes.astype(str)
   if row['trait_type'] =='continuous':
      phenocode_annotate = 'CO_'+ row.category.split(">",)[-1].replace(' ','_')+ '_' +row.modifier+'_' + '_h' + str("{:.5f}".format(row.saige_heritability_EUR)) + '_n' + row.n_cases_full_cohort_both_sexes.astype(str)
   if row['trait_type'] == 'biomarkers':
      phenocode_annotate = 'BI_'+ row.description.replace(' ','_') + '_h' + str("{:.5f}".format(row.saige_heritability_EUR)) + '_n' + row.n_cases_full_cohort_both_sexes.astype(str)
   
   phenocode_annotate_lst = phenocode_annotate_lst +[phenocode_annotate]

#phenocode_annotate_lst = [i.replace('/','_') for i in phenocode_annotate_lst]
ukbb_manifest['phenocode_annotate_lst'] = phenocode_annotate_lst

# Now that we have manually categorized the phenotypes, we need to change some of the annotations we have:
   ## Medical_Conditions --> coding_description column
   ## CA_Medications --> CA_PR_ + coding_description column
   ## Summary_Operations --> coding_description column minus the first 6 strs

ukbb_manifest_main = ukbb_manifest.copy()
new_pheno_annot = list(ukbb_manifest_main.phenocode_annotate_lst)

new_pheno_annot1 = [i.replace("Medical_conditions",str(ukbb_manifest_main['coding_description'].loc[ukbb_manifest_main['phenocode_annotate_lst']==i].values[0])) for i in new_pheno_annot]
new_pheno_annot2 = [i.replace("Medications",'PR_'+str(ukbb_manifest_main['coding_description'].loc[ukbb_manifest_main['phenocode_annotate_lst']==i].values)) for i in new_pheno_annot1]
new_pheno_annot3 = [i.replace("Operations",str(ukbb_manifest_main['coding_description'].loc[ukbb_manifest_main['phenocode_annotate_lst']==i].values)) for i in new_pheno_annot2]
new_pheno_annot4 = [i.replace("Summary_Operations",str(ukbb_manifest_main['coding_description'].loc[ukbb_manifest_main['phenocode_annotate_lst']==i].values)) for i in new_pheno_annot3]
new_pheno_annot = [i.replace("_['",'_') for i in new_pheno_annot4]
new_pheno_annot = [i.replace("']_",'_') for i in new_pheno_annot]
new_pheno_annot = [i.replace(" ",'_') for i in new_pheno_annot]
new_pheno_annot = [i.replace("__",'_') for i in new_pheno_annot]
new_pheno_annot = [i.replace("___",'_') for i in new_pheno_annot]
new_pheno_annot = [i.replace("/",'_') for i in new_pheno_annot]
new_pheno_annot = [i.replace("_Summary_",'_') for i in new_pheno_annot]
new_pheno_annot = [i.replace(")",'') for i in new_pheno_annot]
new_pheno_annot = [i.replace("(",'') for i in new_pheno_annot]
new_pheno_annot = [i[ 0 : i.index("_n")] for i in new_pheno_annot]

ukbb_manifest['new_pheno_annot'] = new_pheno_annot

# Saving this file
ukbb_manifest.to_csv('./Data/ukbb_manifest_filtered_phenos.csv')
