# Once we have our scores, we can combine them into 10 matrices alongside with 10 matrices for the SNP counts of PRSs at each P-value threshold.

import os
import pandas as pd
from os import listdir
from os.path import isfile, join


mypath = './Thesis/UK_BioBank_7k_pheno/PRS_PRSice/PRS/'
file_names_all_score = [f for f in listdir(mypath) if isfile(join(mypath, f))]

## P-val <1
dat_Pt_1 = pd.DataFrame()
for item in file_names_all_score:
    if '.all_score' in item:
        tmp_dat = pd.read_csv(mypath+ '/' +item, sep='\s+')
        cols = tmp_dat.columns.tolist()
        if 'Pt_1' in cols:
            dat_Pt_1[item[:-14]] = tmp_dat['Pt_1']
dat_Pt_1['FID'] = tmp_dat['FID']
dat_Pt_1['IID'] = tmp_dat['IID']
print(tmp_dat.head())
cols = dat_Pt_1.columns.tolist()
cols = cols[-2:] + cols[:-2]
dat_Pt_1 = dat_Pt_1[cols]
print(dat_Pt_1.shape)
path = './Thesis/UK_BioBank_7k_pheno/PRS_PRSice/matrices/'
dat_Pt_1.to_csv(path + 'p_val_1.txt',header=True,index=False)

## P-val <0.1
dat_Pt_0_1 = pd.DataFrame()
for item in file_names_all_score:
    if '.all_score' in item:
        tmp_dat = pd.read_csv(mypath+ '/' +item, sep='\s+')
        cols = tmp_dat.columns.tolist()
        if 'Pt_0.1' in cols:
            dat_Pt_0_1[item[:-14]] = tmp_dat['Pt_0.1']
dat_Pt_0_1['FID'] = tmp_dat['FID']
dat_Pt_0_1['IID'] = tmp_dat['IID']
cols = dat_Pt_0_1.columns.tolist()
cols = cols[-2:] + cols[:-2]
dat_Pt_0_1 = dat_Pt_0_1[cols]
print(dat_Pt_0_1.shape)
path = './ak/Thesis/UK_BioBank_7k_pheno/PRS_PRSice/matrices/'
dat_Pt_0_1.to_csv(path + 'p_val_0.1.txt',header=True,index=False)

## P-val < 0.05
dat_Pt_0_05 = pd.DataFrame()
for item in file_names_all_score:
    if '.all_score' in item:
        tmp_dat = pd.read_csv(mypath+ '/' +item, sep='\s+')
        cols = tmp_dat.columns.tolist()
        if 'Pt_0.05' in cols:
            dat_Pt_0_05[item[:-14]] = tmp_dat['Pt_0.05']
dat_Pt_0_05['FID'] = tmp_dat['FID']
dat_Pt_0_05['IID'] = tmp_dat['IID']
cols = dat_Pt_0_05.columns.tolist()
cols = cols[-2:] + cols[:-2]
dat_Pt_0_05 = dat_Pt_0_05[cols]
print(dat_Pt_0_05.shape)
path = './Thesis/UK_BioBank_7k_pheno/PRS_PRSice/matrices/'
dat_Pt_0_05.to_csv(path + 'p_val_0.05.txt',header=True,index=False)

## P-val < 0.01
dat_Pt_0_01 = pd.DataFrame()
for item in file_names_all_score:
    if '.all_score' in item:
        tmp_dat = pd.read_csv(mypath+ '/' +item, sep='\s+')
        cols = tmp_dat.columns.tolist()
        if 'Pt_0.01' in cols:
            dat_Pt_0_01[item[:-14]] = tmp_dat['Pt_0.01']
dat_Pt_0_01['FID'] = tmp_dat['FID']
dat_Pt_0_01['IID'] = tmp_dat['IID']
cols = dat_Pt_0_01.columns.tolist()
cols = cols[-2:] + cols[:-2]
dat_Pt_0_01 = dat_Pt_0_01[cols]
print(dat_Pt_0_01.shape)
path = './Thesis/UK_BioBank_7k_pheno/PRS_PRSice/matrices/'
dat_Pt_0_01.to_csv(path + 'p_val_0.01.txt',header=True,index=False)

## P-val < 0.001
dat_Pt_0_001 = pd.DataFrame()
for item in file_names_all_score:
    if '.all_score' in item:
        tmp_dat = pd.read_csv(mypath+ '/' +item, sep='\s+')
        cols = tmp_dat.columns.tolist()
        if 'Pt_0.001' in cols:
            dat_Pt_0_001[item[:-14]] = tmp_dat['Pt_0.001']
dat_Pt_0_001['FID'] = tmp_dat['FID']
dat_Pt_0_001['IID'] = tmp_dat['IID']
cols = dat_Pt_0_001.columns.tolist()
cols = cols[-2:] + cols[:-2]
dat_Pt_0_001 = dat_Pt_0_001[cols]
print(dat_Pt_0_001.shape)
path = './Thesis/UK_BioBank_7k_pheno/PRS_PRSice/matrices/'
dat_Pt_0_001.to_csv(path + 'p_val_0.001.txt',header=True,index=False)

## P-val < 0.0001
dat_Pt_0_0001 = pd.DataFrame()
for item in file_names_all_score:
    if '.all_score' in item:
        tmp_dat = pd.read_csv(mypath+ '/' +item, sep='\s+')
        cols = tmp_dat.columns.tolist()
        if 'Pt_0.0001' in cols:
            dat_Pt_0_0001[item[:-14]] = tmp_dat['Pt_0.0001']
dat_Pt_0_0001['FID'] = tmp_dat['FID']
dat_Pt_0_0001['IID'] = tmp_dat['IID']
cols = dat_Pt_0_0001.columns.tolist()
cols = cols[-2:] + cols[:-2]
dat_Pt_0_0001 = dat_Pt_0_0001[cols]
print(dat_Pt_0_0001.shape)
path = './Thesis/UK_BioBank_7k_pheno/PRS_PRSice/matrices/'
dat_Pt_0_0001.to_csv(path + 'p_val_0.0001.txt',header=True,index=False)

## P-val < 1e05
dat_Pt_1e_05 = pd.DataFrame()
for item in file_names_all_score:
    if '.all_score' in item:
        tmp_dat = pd.read_csv(mypath+ '/' +item, sep='\s+')
        cols = tmp_dat.columns.tolist()
        if 'Pt_1e-05' in cols:
            dat_Pt_1e_05[item[:-14]] = tmp_dat['Pt_1e-05']
dat_Pt_1e_05['FID'] = tmp_dat['FID']
dat_Pt_1e_05['IID'] = tmp_dat['IID']
cols = dat_Pt_1e_05.columns.tolist()
cols = cols[-2:] + cols[:-2]
dat_Pt_1e_05 = dat_Pt_1e_05[cols]
print(dat_Pt_1e_05.shape)
path = './Thesis/UK_BioBank_7k_pheno/PRS_PRSice/matrices/'
dat_Pt_1e_05.to_csv(path + 'p_val_1e-05.txt',header=True,index=False)

## P-val < 1e06
dat_Pt_1e_06 = pd.DataFrame()
for item in file_names_all_score:
    if '.all_score' in item:
        tmp_dat = pd.read_csv(mypath+ '/' +item,sep='\s+')
        cols = tmp_dat.columns.tolist()
        if 'Pt_1e-06' in cols:
            dat_Pt_1e_06[item[:-14]] = tmp_dat['Pt_1e-06']
dat_Pt_1e_06['FID'] = tmp_dat['FID']
dat_Pt_1e_06['IID'] = tmp_dat['IID']
cols = dat_Pt_1e_06.columns.tolist()
cols = cols[-2:] + cols[:-2]
dat_Pt_1e_06 = dat_Pt_1e_06[cols]
print(dat_Pt_1e_06.shape)
path = './Thesis/UK_BioBank_7k_pheno/PRS_PRSice/matrices/'
dat_Pt_1e_06.to_csv(path + 'p_val_1e-06.txt',header=True,index=False)

## P-val < 1e07
dat_Pt_1e_07 = pd.DataFrame()
for item in file_names_all_score:
    if '.all_score' in item:
        tmp_dat = pd.read_csv(mypath+ '/' +item,sep='\s+')
        cols = tmp_dat.columns.tolist()
        if 'Pt_1e-07' in cols:
            dat_Pt_1e_07[item[:-14]] = tmp_dat['Pt_1e-07']
dat_Pt_1e_07['FID'] = tmp_dat['FID']
dat_Pt_1e_07['IID'] = tmp_dat['IID']
cols = dat_Pt_1e_07.columns.tolist()
cols = cols[-2:] + cols[:-2]
dat_Pt_1e_07 = dat_Pt_1e_07[cols]
print(dat_Pt_1e_07.shape)
path = './Thesis/UK_BioBank_7k_pheno/PRS_PRSice/matrices/'
dat_Pt_1e_07.to_csv(path + 'p_val_1e-07.txt',header=True,index=False)

## P-val < 5e08
dat_Pt_5e_08 = pd.DataFrame()
for item in file_names_all_score:
    if '.all_score' in item:
        tmp_dat = pd.read_csv(mypath+ '/' +item,sep='\s+')
        cols = tmp_dat.columns.tolist()
        if 'Pt_5e-08' in cols:
            dat_Pt_5e_08[item[:-14]] = tmp_dat['Pt_5e-08']
dat_Pt_5e_08['FID'] = tmp_dat['FID']
dat_Pt_5e_08['IID'] = tmp_dat['IID']
cols = dat_Pt_5e_08.columns.tolist()
cols = cols[-2:] + cols[:-2]
dat_Pt_5e_08 = dat_Pt_5e_08[cols]
print(dat_Pt_5e_08.shape)
path = './Thesis/UK_BioBank_7k_pheno/PRS_PRSice/matrices/'
dat_Pt_5e_08.to_csv(path + 'p_val_5e-08.txt',header=True,index=False)

#############################################################################################
# Now, we get the SNP count dataset for PRS of each phenotype at different P-value thresholds:

import os
import pandas as pd
from os import listdir
from os.path import isfile, join


mypath = './Thesis/UK_BioBank_7k_pheno/PRS_PRSice/PRS/'
file_names_snp = [f for f in listdir(mypath) if isfile(join(mypath, f))]

dat_snp_count = pd.DataFrame()
for item in file_names_snp:
    if '.tsv.snp' in item:
        tmp_dat = pd.read_csv(mypath+ '/' +item,sep='\s+')
        n10 = tmp_dat[tmp_dat['P']<=5e-08].shape[0]
        n9 = tmp_dat[(tmp_dat['P']<=1e-07)].shape[0]
        n8 = tmp_dat[(tmp_dat['P']<=1e-06)].shape[0]
        n7 = tmp_dat[(tmp_dat['P']<=1e-05)].shape[0]
        n6 = tmp_dat[(tmp_dat['P']<=0.0001)].shape[0]
        n5 = tmp_dat[(tmp_dat['P']<=0.001)].shape[0]
        n4 = tmp_dat[(tmp_dat['P']<=0.01)].shape[0]
        n3 = tmp_dat[(tmp_dat['P']<=0.05)].shape[0]
        n2 = tmp_dat[(tmp_dat['P']<=0.1)].shape[0]
        n1 = tmp_dat[(tmp_dat['P']<=1)].shape[0]
        dat_snp_count[item[:-8]] = [n1,n2,n3,n4,n5,n6,n7,n8,n9,n10]
path = './Thesis/UK_BioBank_7k_pheno/PRS_PRSice/matrices/'
dat_snp_count.to_csv(path + 'snp_count_cumulative.txt',header=True,index=False,sep='\t')
