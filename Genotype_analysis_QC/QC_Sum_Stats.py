import pandas as pd
import numpy as np
import os
from scipy.stats import chisquare
import scipy.stats

path = '/external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/PRS_plink/sum_stats_1' ######################## 
m = [s for s in os.listdir(path) if ".tsv" in s][0]
dat = path +'/'+ m
sum_stat_dat = pd.read_csv(dat,sep='\t')

path1 = '/external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/'

rsid_data = pd.read_table(path1+'full_variant_qc_metrics.txt',sep='\t')

varid = sum_stat_dat.chr.apply(str)+':'+sum_stat_dat.pos.apply(str)+'_'+sum_stat_dat.ref+'_'+sum_stat_dat.alt

sum_stat_dat['varid']=varid

# Adding the rsid of each SNPs
def update_type(t1, t2, dropna=False):
    return t1.map(t2).dropna() if dropna else t1.map(t2).fillna(t1)
rsid = update_type(sum_stat_dat.varid, rsid_data.set_index('varid').rsid, dropna=True)
sum_stat_dat['SNP']=rsid
sum_stat_dat = sum_stat_dat.dropna()

# Removing duplicated SNPs
sum_stat_dat = sum_stat_dat.drop_duplicates(subset='SNP',keep='last')

# Removing SNPs with a minor allele count less than 20
sum_stat_dat = sum_stat_dat[sum_stat_dat['low_confidence_EUR'] == False]

# Change all alleles to upper case
sum_stat_dat['ref']= sum_stat_dat['ref'].str.upper()
sum_stat_dat['alt']= sum_stat_dat['alt'].str.upper()

# renaming the columns of the summary stats file

sum_stat_dat = sum_stat_dat.rename(columns={'ref':'A1','alt':'A2','pval_EUR':'P','beta_EUR':'BETA'})

## Saving the updated summary stat file

path_to_save = '/external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/PRS_plink/sum_stats_1/' + m + '.txt' ####################
sum_stat_dat.to_csv(path_to_save, header=True, index=False, sep='\t', mode='a')

# To save the total number of variants left in a txt file
path_to_save_1 = '/external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/number_of_variants/' + m[:-4] + '.txt'
text_file = open(path_to_save_1, "w")
text_file.write('{}'.format(sum_stat_dat.shape[0]))
text_file.close()

# calculating GClambda

pval = np.array(list(sum_stat_dat.P))
pval = 1 - pval
chisq = scipy.stats.chi2.ppf(pval,1)
lambdagc = np.median(chisq)/scipy.stats.chi2.ppf(0.5, 1)

## To save the gclambda value
path_to_save_2 = '/external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/PRS_plink/GCLambda/' + m[:-4] + '.txt'
text_file = open(path_to_save_2, "w")
text_file.write('{}'.format(lambdagc))
text_file.close()
