# Once we have our scores, we can combine them into 10 matrices alongside with 10 matrices for the SNP counts of PRSs at each P-value threshold.

import os
import pandas as pd
from os import listdir
from os.path import isfile, join


mypath = '/external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/PRS_plink/PRS/'
file_names_all_score = [f for f in listdir(mypath) if isfile(join(mypath, f))]

###############################
#Creating the PRS matrix for P <= 1
dat_P_1 = pd.DataFrame()
dat_snp_count_P_1 = pd.DataFrame()
for item in file_names_all_score:
    if '.P_1.profile' in item:
        tmp_dat = pd.read_csv(mypath +item, sep='\s+')
        dat_P_1[item[:-16]] = tmp_dat['SCORE']
        dat_snp_count_P_1[item[:-16]] = tmp_dat['CNT']
dat_P_1['FID'] = tmp_dat['FID']
dat_P_1['IID'] = tmp_dat['IID']
dat_snp_count_P_1['FID'] = tmp_dat['FID']
dat_snp_count_P_1['IID'] = tmp_dat['IID']

cols = dat_P_1.columns.tolist()
cols = cols[-2:] + cols[:-2]
dat_P_1 = dat_P_1[cols]
print(dat_P_1.shape)
print(dat_snp_count_P_1.shape)
path = '/external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/PRS_plink/matrices/'
dat_P_1.to_csv(path + 'p_val_1.txt',header=True,index=False)
dat_snp_count_P_1.to_csv(path + 'snp_count_p_val_1.txt',header=True,index=False)

###############################
#Creating the PRS matrix for P <= 0.1
dat_P_0_1 = pd.DataFrame()
dat_snp_count_P_0_1 = pd.DataFrame()
for item in file_names_all_score:
    if '.P_0.1.profile' in item:
        tmp_dat = pd.read_csv(mypath +item, sep='\s+')
        dat_P_0_1[item[:-18]] = tmp_dat['SCORE']
        dat_snp_count_P_0_1[item[:-18]] = tmp_dat['CNT']
dat_P_0_1['FID'] = tmp_dat['FID']
dat_P_0_1['IID'] = tmp_dat['IID']
dat_snp_count_P_0_1['FID'] = tmp_dat['FID']
dat_snp_count_P_0_1['IID'] = tmp_dat['IID']

cols = dat_P_0_1.columns.tolist()
cols = cols[-2:] + cols[:-2]
dat_P_1 = dat_P_0_1[cols]
print(dat_P_0_1.shape)
print(dat_snp_count_P_0_1.shape)
path = '/external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/PRS_plink/matrices/'
dat_P_0_1.to_csv(path + 'p_val_0_1.txt',header=True,index=False)
dat_snp_count_P_0_1.to_csv(path + 'snp_count_p_val_0_1.txt',header=True,index=False)

###############################
#Creating the PRS matrix for P <= 0.05
dat_P_0_05 = pd.DataFrame()
dat_snp_count_P_0_05 = pd.DataFrame()
for item in file_names_all_score:
    if '.P_0.05.profile' in item:
        tmp_dat = pd.read_csv(mypath +item, sep='\s+')
        dat_P_0_05[item[:-19]] = tmp_dat['SCORE']
        dat_snp_count_P_0_05[item[:-19]] = tmp_dat['CNT']
dat_P_0_05['FID'] = tmp_dat['FID']
dat_P_0_05['IID'] = tmp_dat['IID']
dat_snp_count_P_0_05['FID'] = tmp_dat['FID']
dat_snp_count_P_0_05['IID'] = tmp_dat['IID']

cols = dat_P_0_05.columns.tolist()
cols = cols[-2:] + cols[:-2]
dat_P_0_05 = dat_P_0_05[cols]
print(dat_P_0_05.shape)
print(dat_snp_count_P_0_05.shape)
path = '/external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/PRS_plink/matrices/'
dat_P_0_05.to_csv(path + 'p_val_0_05.txt',header=True,index=False)
dat_snp_count_P_0_05.to_csv(path + 'snp_count_p_val_0_05.txt',header=True,index=False)

###############################
#Creating the PRS matrix for P <= 0.01
dat_P_0_01 = pd.DataFrame()
dat_snp_count_P_0_01 = pd.DataFrame()
for item in file_names_all_score:
    if '.P_0.01.profile' in item:
        tmp_dat = pd.read_csv(mypath +item, sep='\s+')
        dat_P_0_01[item[:-21]] = tmp_dat['SCORE']
        dat_snp_count_P_0_01[item[:-21]] = tmp_dat['CNT']
dat_P_0_01['FID'] = tmp_dat['FID']
dat_P_0_01['IID'] = tmp_dat['IID']
dat_snp_count_P_0_01['FID'] = tmp_dat['FID']
dat_snp_count_P_0_01['IID'] = tmp_dat['IID']

cols = dat_P_0_01.columns.tolist()
cols = cols[-2:] + cols[:-2]
dat_P_0_01 = dat_P_0_01[cols]
print(dat_P_0_01.shape)
print(dat_snp_count_P_0_01.shape)
path = '/external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/PRS_plink/matrices/'
dat_P_0_01.to_csv(path + 'p_val_0_01.txt',header=True,index=False)
dat_snp_count_P_0_01.to_csv(path + 'snp_count_p_val_0_01.txt',header=True,index=False)

###############################
#Creating the PRS matrix for P <= 0.001
dat_P_0_001 = pd.DataFrame()
dat_snp_count_P_0_001 = pd.DataFrame()
for item in file_names_all_score:
    if '.P_0.001.profile' in item:
        tmp_dat = pd.read_csv(mypath +item, sep='\s+')
        dat_P_0_001[item[:-20]] = tmp_dat['SCORE']
        dat_snp_count_P_0_001[item[:-20]] = tmp_dat['CNT']
dat_P_0_001['FID'] = tmp_dat['FID']
dat_P_0_001['IID'] = tmp_dat['IID']
dat_snp_count_P_0_001['FID'] = tmp_dat['FID']
dat_snp_count_P_0_001['IID'] = tmp_dat['IID']

cols = dat_P_0_001.columns.tolist()
cols = cols[-2:] + cols[:-2]
dat_P_0_001 = dat_P_0_001[cols]
print(dat_P_0_001.shape)
print(dat_snp_count_P_0_001.shape)
path = '/external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/PRS_plink/matrices/'
dat_P_0_001.to_csv(path + 'p_val_0_001.txt',header=True,index=False)
dat_snp_count_P_0_001.to_csv(path + 'snp_count_p_val_0_001.txt',header=True,index=False)

###############################
#Creating the PRS matrix for P <= 0.0001
dat_P_0_0001 = pd.DataFrame()
dat_snp_count_P_0_0001 = pd.DataFrame()
for item in file_names_all_score:
    if '.P_0.0001.profile' in item:
        tmp_dat = pd.read_csv(mypath +item, sep='\s+')
        dat_P_0_0001[item[:-21]] = tmp_dat['SCORE']
        dat_snp_count_P_0_0001[item[:-21]] = tmp_dat['CNT']
dat_P_0_0001['FID'] = tmp_dat['FID']
dat_P_0_0001['IID'] = tmp_dat['IID']
dat_snp_count_P_0_0001['FID'] = tmp_dat['FID']
dat_snp_count_P_0_0001['IID'] = tmp_dat['IID']

cols = dat_P_0_0001.columns.tolist()
cols = cols[-2:] + cols[:-2]
dat_P_0_0001 = dat_P_0_0001[cols]
print(dat_P_0_0001.shape)
print(dat_snp_count_P_0_0001.shape)
path = '/external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/PRS_plink/matrices/'
dat_P_0_0001.to_csv(path + 'p_val_0_0001.txt',header=True,index=False)
dat_snp_count_P_0_0001.to_csv(path + 'snp_count_p_val_0_0001.txt',header=True,index=False)

###############################
#Creating the PRS matrix for P <= 1e-05
dat_P_1e_05 = pd.DataFrame()
dat_snp_count_P_1e_05 = pd.DataFrame()
for item in file_names_all_score:
    if '.P_1e-05.profile' in item:
        tmp_dat = pd.read_csv(mypath +item, sep='\s+')
        dat_P_1e_05[item[:-20]] = tmp_dat['SCORE']
        dat_snp_count_P_1e_05[item[:-20]] = tmp_dat['CNT']
dat_P_1e_05['FID'] = tmp_dat['FID']
dat_P_1e_05['IID'] = tmp_dat['IID']
dat_snp_count_P_1e_05['FID'] = tmp_dat['FID']
dat_snp_count_P_1e_05['IID'] = tmp_dat['IID']

cols = dat_P_1e_05.columns.tolist()
cols = cols[-2:] + cols[:-2]
dat_P_1e_05 = dat_P_1e_05[cols]
print(dat_P_1e_05.shape)
print(dat_snp_count_P_1e_05.shape)
path = '/external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/PRS_plink/matrices/'
dat_P_1e_05.to_csv(path + 'p_val_1e_05.txt',header=True,index=False)
dat_snp_count_P_1e_05.to_csv(path + 'snp_count_p_val_1e_05.txt',header=True,index=False)

###############################
#Creating the PRS matrix for P <= 1e-06
dat_P_1e_06 = pd.DataFrame()
dat_snp_count_P_1e_06 = pd.DataFrame()
for item in file_names_all_score:
    if '.P_1e-06.profile' in item:
        tmp_dat = pd.read_csv(mypath +item, sep='\s+')
        dat_P_1e_06[item[:-20]] = tmp_dat['SCORE']
        dat_snp_count_P_1e_06[item[:-20]] = tmp_dat['CNT']
dat_P_1e_06['FID'] = tmp_dat['FID']
dat_P_1e_06['IID'] = tmp_dat['IID']
dat_snp_count_P_1e_06['FID'] = tmp_dat['FID']
dat_snp_count_P_1e_06['IID'] = tmp_dat['IID']

cols = dat_P_1e_06.columns.tolist()
cols = cols[-2:] + cols[:-2]
dat_P_1e_06 = dat_P_1e_06[cols]
print(dat_P_1e_06.shape)
print(dat_snp_count_P_1e_06.shape)
path = '/external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/PRS_plink/matrices/'
dat_P_1e_06.to_csv(path + 'p_val_1e_06.txt',header=True,index=False)
dat_snp_count_P_1e_06.to_csv(path + 'snp_count_p_val_1e_06.txt',header=True,index=False)

###############################
#Creating the PRS matrix for P <= 1e-07
dat_P_1e_07 = pd.DataFrame()
dat_snp_count_P_1e_07 = pd.DataFrame()
for item in file_names_all_score:
    if '.P_1e-07.profile' in item:
        tmp_dat = pd.read_csv(mypath +item, sep='\s+')
        dat_P_1e_07[item[:-20]] = tmp_dat['SCORE']
        dat_snp_count_P_1e_07[item[:-20]] = tmp_dat['CNT']
dat_P_1e_07['FID'] = tmp_dat['FID']
dat_P_1e_07['IID'] = tmp_dat['IID']
dat_snp_count_P_1e_07['FID'] = tmp_dat['FID']
dat_snp_count_P_1e_07['IID'] = tmp_dat['IID']

cols = dat_P_1e_07.columns.tolist()
cols = cols[-2:] + cols[:-2]
dat_P_1e_07 = dat_P_1e_07[cols]
print(dat_P_1e_07.shape)
print(dat_snp_count_P_1e_07.shape)
path = '/external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/PRS_plink/matrices/'
dat_P_1e_07.to_csv(path + 'p_val_1e_07.txt',header=True,index=False)
dat_snp_count_P_1e_07.to_csv(path + 'snp_count_p_val_1e_07.txt',header=True,index=False)

###############################
#Creating the PRS matrix for P <= 5e-08
dat_P_5e_08 = pd.DataFrame()
dat_snp_count_P_5e_08 = pd.DataFrame()
for item in file_names_all_score:
    if '.P_5e-08.profile' in item:
        tmp_dat = pd.read_csv(mypath +item, sep='\s+')
        dat_P_5e_08[item[:-20]] = tmp_dat['SCORE']
        dat_snp_count_P_5e_08[item[:-20]] = tmp_dat['CNT']
dat_P_5e_08['FID'] = tmp_dat['FID']
dat_P_5e_08['IID'] = tmp_dat['IID']
dat_snp_count_P_5e_08['FID'] = tmp_dat['FID']
dat_snp_count_P_5e_08['IID'] = tmp_dat['IID']

cols = dat_P_5e_08.columns.tolist()
cols = cols[-2:] + cols[:-2]
dat_P_5e_08 = dat_P_5e_08[cols]
print(dat_P_5e_08.shape)
print(dat_snp_count_P_5e_08.shape)
path = '/external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/PRS_plink/matrices/'
dat_P_5e_08.to_csv(path + 'p_val_5e_08.txt',header=True,index=False)
dat_snp_count_P_5e_08.to_csv(path + 'snp_count_p_val_5e_08.txt',header=True,index=False)
