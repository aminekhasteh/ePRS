Generating PRS (PLINK)
================

``` bash
#!/bin/bash --login
#SBATCH --time=50:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --job-name=prsukbb1
#SBATCH --output prs_ukbb_overlap1.out.txt

cd /external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/PRS_plink/sum_stats_1 ########
module load PYTHON/3.6
module load PLINK2/1.90b3.46

START=$(date +%s.%N)
while read url <&3 && read pheno_names <&4; do
    wget $url
    gunzip -c *.bgz > $pheno_names
      rm *.bgz
      python /external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/PRS_plink/py_files/clean_data_overlap_1.py ##########
    rm *.tsv
    plink \
      --bfile /external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/QC_Geno/ROSMAP/ROSMAP_QC_rsid/ROSMAP.QC \
      --clump-p1 1 \
      --clump-r2 0.2 \
      --clump-kb 500 \
      --clump $pheno_names.txt \
      --clump-snp-field SNP \
      --clump-field P \
      --out ROSMAP
    awk 'NR!=1{print $3}' ROSMAP.clumped >  ROSMAP.valid.snp
    awk '{print $12,$9}' $pheno_names.txt > /external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/PRS_plink/snp_pval_all_pheno/$pheno_names.SNP.pvalue
    awk '{print $12,$7}' $pheno_names.txt > /external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/PRS_plink/snp_beta_all_pheno/$pheno_names.SNP.beta
    awk '{print $12,$8}' $pheno_names.txt > /external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/PRS_plink/snp_se_all_pheno/$pheno_names.SNP.se
    awk '{print $12,$5}' $pheno_names.txt > /external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/PRS_plink/snp_afcase_all_pheno/$pheno_names.SNP.afcase
    awk '{print $12,$6}' $pheno_names.txt > /external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/PRS_plink/snp_afcontrol_all_pheno/$pheno_names.SNP.afcontrol
    plink \
      --bfile /external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/QC_Geno/ROSMAP/ROSMAP_QC_rsid/ROSMAP.QC \
      --score $pheno_names.txt 12 3 7 header \
      --q-score-range /external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/PRS_plink/range_list /external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/PRS_plink/snp_pval_all_pheno/$pheno_names.SNP.pvalue \
      --extract ROSMAP.valid.snp \
      --out /external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/PRS_plink/PRS/$pheno_names
    rm *.txt
    rm *.clumped
    rm *.log
    rm *.nopred
    rm *.nosex
    rm *.valid.snp
    rm *.pvalue
done 3< /external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/url/url_aa.txt 4</external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/pheno_names/pheno_names_aa.txt ########
END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo $DIFF
```
