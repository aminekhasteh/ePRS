# A) Libraries -------------------------------
library(ggplot2)
library(plotly)
library(hrbrthemes)
library(RColorBrewer)
library(tidyverse)
library(latex2exp)
library(WGCNA)
library(reshape2)
library(ggrepel)
library(stringr)
library(ggbeeswarm)
library(hrbrthemes)
library(ggpubr)
library(cowplot)
library(gplots)
library(ggdendro)
library(ggraph)
library(stringr)
library(data.table)
library(rms)
library(Hmisc)
library(summarytools)
library(factoextra)
library(epicalc)
library(lmtest)
library(broom)
library(ggpol)
library(RColorBrewer) 
Study = "ROSMAP"

# B) Functions ------
# Using this function in ggplots to avoid repeating labels
stat_box_data <- function(y) {
                return( 
                                data.frame(
                                                y = 0.5+1.1*max(y),  #may need to modify this depending on your data
                                                label = paste('M =', length(y), '\n')
                                )
                )
}

## These functions help set up our data for heatmaps:
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
                cormat[upper.tri(cormat)] <- NA
                return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
                cormat[lower.tri(cormat)]<- NA
                return(cormat)
}
# reorder the correlations
reorder_cormat <- function(cormat){
                # Use correlation between variables as distance
                dd <- as.dist((1-cormat)/2)
                hc <- hclust(dd)
                cormat <-cormat[hc$order, hc$order]
}

# This function will generate association analysis on the WGCNA hubs and AD and AD-related pathology
# Using the rms package and running .236 bootstrap method
# 95% CI (cite): https://discourse.datamethods.org/t/confidence-intervals-for-bootstrap-validated-bias-corrected-performance-estimates/1990/9
OLS_heatmap_dat_GEN_WGCNA <- function(dat,
                                      Study=c("ROSMAP","ADNI"),
                                      GenoType=c("No_APOE","No_MHC","No_MHC_APOE","With_MHC_APOE")){
                
                path_wgcna <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/WGCNA/",GenoType)
                path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/Resid_PRS/",GenoType)
                wgcna_dat_results <- readRDS(file=paste0(path_wgcna,"/wgcnaresults_all",".rds"))
                prs <- readRDS(paste0(path,"/Residual_results_all_p-vals.rds"))
                ROSmaster <- readRDS("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/ROSMAP_Phenotype/ROSmaster.rds")
                # Changing Cogdx variable:
                ROSmaster$cogdx[which((ROSmaster$cogdx==2)|(ROSmaster$cogdx==3)|(ROSmaster$cogdx==5)|(ROSmaster$cogdx==6))] <- NA
                ROSmaster$cogdx[which((ROSmaster$cogdx==4))] <- 2
                
                # form. covatiates only, base_form. IGAP + covariates, form_full. IGAP + eigenPS + covariates
                indepvec.pathology <- c("cogdx","amyloid_sqrt","tangles_sqrt")
                indepvec.cognition <- c("cogn_global_random_slope","cogn_globaln_lv")
                
                index <- 1
                alpha_val = NULL #-----> Y
                r2validated_base <- NULL
                r2validated <- NULL
                r2validated_full <- NULL
                r2validated_base_CI_L <- NULL
                r2validated_CI_L <- NULL
                r2validated_full_CI_L <- NULL
                r2validated_base_CI_U <- NULL
                r2validated_CI_U <- NULL
                r2validated_full_CI_U <- NULL
                nvalues <- NULL
                phenovalues <- NULL
                pcvalues <- NULL
                P_likelihood <- NULL
                coeff_IGAP <- NULL
                coeff_eigenPS <- NULL
                modulevalues <- NULL
                
                for(alpha in names(prs)){
                                # Defining the datasets
                                results <- wgcna_dat_results[[alpha]]$assoc_res
                                net <- wgcna_dat_results[[alpha]]$net
                                df <- prs[[alpha]]$residuals
                                rownames(df) <- df$IID
                                df$IID <- NULL
                                
                                #### Modelling against LOAD 
                                # merge data
                                MEs <- net$MEs
                                MEs$IID <- rownames(MEs)
                                df2 <- df
                                df2$IID <- rownames(df)
                                md2 <- merge(df2,ROSmaster,by="IID")
                                md3 <- merge(md2,MEs,by="IID")
                                
                                for (pheno in c(indepvec.pathology,indepvec.cognition)){
                                                
                                                ph=pheno
                                                ab <- dat %>% filter(phenotype==ph,
                                                                     PVAL==alpha,
                                                                     geno==GenoType)
                                                print(dim(ab))
                                                
                                                if(dim(ab)[1]==0){next}
                                                else{
                                                                for(j in 1:dim(ab)[1]){
                                                                                egene <- paste0("ME",as.character(ab$Module[j]))
                                                                                print(paste(egene,pheno,alpha))
                                                                                if(egene=="grey"){next}
                                                                                if (pheno=="cogdx"){
                                                                                base_form <- formula(paste(pheno,"~","LOAD+msex+age_death+pmi"))
                                                                                form <- formula(paste(pheno,"~","+msex+age_death+pmi"))
                                                                                form_full <- formula(paste(pheno,"~",egene,"+LOAD+msex+age_death+pmi"))
                                                                                mod_base <- lrm(data=md3,base_form, y=T, x=T)
                                                                                mod <- lrm(data=md3,form, y=T, x=T)
                                                                                mod_full <- lrm(data=md3,form_full, y=T, x=T)
                                                                                P_likelihood_test <- rms::lrtest(mod_base, mod_full)$stats[3]
                                                                                valid_base <- rms::validate(mod_base, method=".632",B=200,bw=TRUE) # valid_bsae[2,3] --> R2
                                                                                valid <- rms::validate(mod, method=".632",B=200,bw=TRUE)
                                                                                valid_full <- rms::validate(mod_full, method=".632",B=200,bw=TRUE)
                                                                                
                                                                                R2_base <- valid_base['Dxy', 'index.corrected'] * 0.5 + 0.5
                                                                                R2 <- valid['Dxy', 'index.corrected'] * 0.5 + 0.5
                                                                                R2_full <- valid_full['Dxy', 'index.corrected'] * 0.5 + 0.5
                                                                                
                                                                                # Because the AUC/R^2 values are so close, we calculate the 95% CI here:
                                                                                dxy <- NULL
                                                                                dxy1 <- NULL
                                                                                dxy2 <- NULL
                                                                                for(i in 1 : 50) {
                                                                                                f <- mod_base
                                                                                                f1 <- mod
                                                                                                f2 <- mod_full
                                                                                                n <- nrow(md3)                
                                                                                                g <- update(f, subset=sample(1 : n, n, replace=TRUE))
                                                                                                v <- rms::validate(g, method=".632",B=200,bw=TRUE)
                                                                                                g1 <- update(f1, subset=sample(1 : n, n, replace=TRUE))
                                                                                                v1 <- rms::validate(g1, method=".632",B=200,bw=TRUE)
                                                                                                g2 <- update(f2, subset=sample(1 : n, n, replace=TRUE))
                                                                                                v2 <- rms::validate(g2, method=".632",B=200,bw=TRUE)
                                                                                                dxy[i] <- v['Dxy', 'index.corrected'] * 0.5 + 0.5
                                                                                                dxy1[i] <- v1['Dxy', 'index.corrected'] * 0.5 + 0.5
                                                                                                dxy2[i] <- v2['Dxy', 'index.corrected'] * 0.5 + 0.5
                                                                                }
                                                                                R2_base_CI_L <-  as.numeric(quantile(dxy, c(.025, .975),na.rm=T)[1])
                                                                                R2_base_CI_U <-  as.numeric(quantile(dxy, c(.025, .975),na.rm=T)[2])
                                                                                R2_CI_L <-  as.numeric(quantile(dxy1, c(.025, .975),na.rm=T)[1])
                                                                                R2_CI_U <-  as.numeric(quantile(dxy1, c(.025, .975),na.rm=T)[2])
                                                                                R2_full_CI_L <-  as.numeric(quantile(dxy2, c(.025, .975),na.rm=T)[1])
                                                                                R2_full_CI_U <-  as.numeric(quantile(dxy2, c(.025, .975),na.rm=T)[2])
                                                                } 
                                                                                if (pheno=="cogn_global_random_slope"|pheno=="cogn_globaln_lv"){
                                                                                base_form <- formula(paste(pheno,"~","LOAD+msex+educ"))
                                                                                form <- formula(paste(pheno,"~","+msex+educ"))
                                                                                form_full <- formula(paste(pheno,"~",egene,"+LOAD+msex+educ"))
                                                                                mod_base <- ols(data=md3,base_form, y=T, x=T)
                                                                                mod <- ols(data=md3,form, y=T, x=T)
                                                                                mod_full <- ols(data=md3,form_full, y=T, x=T)
                                                                                P_likelihood_test <- rms::lrtest(mod_base, mod_full)$stats[3]
                                                                                valid_base <- rms::validate(mod_base, method=".632",B=200,bw=TRUE) # valid_base[1,3] --> R2
                                                                                valid <- rms::validate(mod, method=".632",B=200,bw=TRUE)
                                                                                valid_full <- rms::validate(mod_full, method=".632",B=200,bw=TRUE)
                                                                                
                                                                                R2_base <- valid_base['R-square', 'index.corrected']
                                                                                R2 <- valid['R-square', 'index.corrected']
                                                                                R2_full <- valid_full['R-square', 'index.corrected']
                                                                                
                                                                                # Because the AUC/R^2 values are so close, we calculate the 95% CI here:
                                                                                dxy <- NULL
                                                                                dxy1 <- NULL
                                                                                dxy2 <- NULL
                                                                                for(i in 1 : 50) {
                                                                                                f <- mod_base
                                                                                                f1 <- mod
                                                                                                f2 <- mod_full
                                                                                                n <- nrow(md3)                
                                                                                                g <- update(f, subset=sample(1 : n, n, replace=TRUE))
                                                                                                v <- rms::validate(g, method=".632",B=200,bw=TRUE)
                                                                                                g1 <- update(f1, subset=sample(1 : n, n, replace=TRUE))
                                                                                                v1 <- rms::validate(g1, method=".632",B=200,bw=TRUE)
                                                                                                g2 <- update(f2, subset=sample(1 : n, n, replace=TRUE))
                                                                                                v2 <- rms::validate(g2, method=".632",B=200,bw=TRUE)
                                                                                                dxy[i] <-   v['R-square', 'index.corrected']
                                                                                                dxy1[i] <- v1['R-square', 'index.corrected']
                                                                                                dxy2[i] <- v2['R-square', 'index.corrected']
                                                                                }
                                                                                R2_base_CI_L <-  as.numeric(quantile(dxy, c(.025, .975))[1])
                                                                                R2_base_CI_U <-  as.numeric(quantile(dxy, c(.025, .975))[2])
                                                                                R2_CI_L <-  as.numeric(quantile(dxy1, c(.025, .975))[1])
                                                                                R2_CI_U <-  as.numeric(quantile(dxy1, c(.025, .975))[2])
                                                                                R2_full_CI_L <-  as.numeric(quantile(dxy2, c(.025, .975))[1])
                                                                                R2_full_CI_U <-  as.numeric(quantile(dxy2, c(.025, .975))[2])
                                                                } 
                                                                                if (pheno=="tangles_sqrt"|pheno=="amyloid_sqrt"){
                                                                                base_form <- formula(paste(pheno,"~","LOAD+msex+age_death+pmi"))
                                                                                form <- formula(paste(pheno,"~","+msex+age_death+pmi"))
                                                                                form_full <- formula(paste(pheno,"~",egene,"+LOAD+msex+age_death+pmi"))
                                                                                mod_base <- ols(data=md3,base_form, y=T, x=T)
                                                                                mod <- ols(data=md3,form, y=T, x=T)
                                                                                mod_full <- ols(data=md3,form_full, y=T, x=T)
                                                                                P_likelihood_test <- rms::lrtest(mod_base, mod_full)$stats[3]
                                                                                valid_base <- rms::validate(mod_base, method=".632",B=200,bw=TRUE) # valid_base[1,3] --> R2
                                                                                valid <- rms::validate(mod, method=".632",B=200,bw=TRUE)
                                                                                valid_full <- rms::validate(mod_full, method=".632",B=200,bw=TRUE)
                                                                                
                                                                                R2_base <- valid_base['R-square', 'index.corrected']
                                                                                R2 <- valid['R-square', 'index.corrected']
                                                                                R2_full <- valid_full['R-square', 'index.corrected']
                                                                                
                                                                                # Because the AUC/R^2 values are so close, we calculate the 95% CI here:
                                                                                dxy <- NULL
                                                                                dxy1 <- NULL
                                                                                dxy2 <- NULL
                                                                                for(i in 1 : 50) {
                                                                                                f <- mod_base
                                                                                                f1 <- mod
                                                                                                f2 <- mod_full
                                                                                                n <- nrow(md3)                
                                                                                                g <- update(f, subset=sample(1 : n, n, replace=TRUE))
                                                                                                v <- rms::validate(g, method=".632",B=200,bw=TRUE)
                                                                                                g1 <- update(f1, subset=sample(1 : n, n, replace=TRUE))
                                                                                                v1 <- rms::validate(g1, method=".632",B=200,bw=TRUE)
                                                                                                g2 <- update(f2, subset=sample(1 : n, n, replace=TRUE))
                                                                                                v2 <- rms::validate(g2, method=".632",B=200,bw=TRUE)
                                                                                                dxy[i] <-   v['R-square', 'index.corrected']
                                                                                                dxy1[i] <- v1['R-square', 'index.corrected']
                                                                                                dxy2[i] <- v2['R-square', 'index.corrected']
                                                                                }
                                                                                R2_base_CI_L <-  as.numeric(quantile(dxy, c(.025, .975))[1])
                                                                                R2_base_CI_U <-  as.numeric(quantile(dxy, c(.025, .975))[2])
                                                                                R2_CI_L <-  as.numeric(quantile(dxy1, c(.025, .975))[1])
                                                                                R2_CI_U <-  as.numeric(quantile(dxy1, c(.025, .975))[2])
                                                                                R2_full_CI_L <-  as.numeric(quantile(dxy2, c(.025, .975))[1])
                                                                                R2_full_CI_U <-  as.numeric(quantile(dxy2, c(.025, .975))[2])
                                                                }
                                                                                r2validated_base[index] <- R2_base
                                                                                r2validated[index] <- R2
                                                                                r2validated_full[index] <- R2_full
                                                                                r2validated_base_CI_L[index] <- R2_base_CI_L
                                                                                r2validated_CI_L[index] <- R2_CI_L
                                                                                r2validated_full_CI_L[index] <- R2_full_CI_L     
                                                                                r2validated_base_CI_U[index] <- R2_base_CI_U      
                                                                                r2validated_CI_U[index] <- R2_CI_U
                                                                                r2validated_full_CI_U[index] <- R2_full_CI_U
                                                                                r2validated_base[index] <- R2_base
                                                                                r2validated[index] <- R2
                                                                                r2validated_full[index] <- R2_full
                                                                                nvalues[index] <- length(mod$y)
                                                                                modulevalues[index] <- egene
                                                                                phenovalues[index] <- pheno
                                                                                alpha_val[index] <- alpha
                                                                                coeff_IGAP[index] <- mod_full$coefficients[3]
                                                                                coeff_eigenPS[index] <- mod_full$coefficients[2]
                                                                                P_likelihood[index] <- P_likelihood_test
                                                                                index <- index + 1

                                                                }
                                                }
                                }
                                
                                
                                assocres <- data.frame(alpha_level=alpha_val,
                                                       pheno=phenovalues,
                                                       modulevalues=modulevalues,
                                                       r2val=r2validated,
                                                       r2val_base=r2validated_base,
                                                       r2val_full=r2validated_full,
                                                       coeff_IGAP=coeff_IGAP,
                                                       r2validated_base_CI_L=r2validated_base_CI_L,
                                                       r2validated_CI_L=r2validated_CI_L,
                                                       r2validated_full_CI_L=r2validated_full_CI_L,
                                                       r2validated_base_CI_U=r2validated_base_CI_U,
                                                       r2validated_CI_U=r2validated_CI_U,
                                                       r2validated_full_CI_U=r2validated_full_CI_U,
                                                       coeff_eigenPS=coeff_eigenPS,
                                                       P_likelihood=P_likelihood,
                                                       n=nvalues)
                                write.csv(assocres,paste0(path_wgcna,"/assocres_compare_hub_heatmap_dat6.csv"),row.names = F)
                }
}



## 1) ROSMAP demographic Table ------

GenoType=c("No_APOE","No_MHC","No_MHC_APOE","With_MHC_APOE")

# Reading ROS/MAP phenotype dataset
ROSmaster <- readRDS("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/ROSMAP_Phenotype/ROSmaster.rds")
# We want to create some tables in understating the ROS/MAP participants
prs <- fread("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/PRS/No_APOE/p_val_0.0001.txt")
geno_pcs <- read.table("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/PCA_Genotype/geno_qc.eigenvec_new.txt",header=F)
names(geno_pcs) <- c("FID","IID",paste0("genoPC",seq(1:10)))
ROSmaster <- ROSmaster %>% merge(.,prs,by="IID") %>%
                merge(.,geno_pcs,by="IID") %>%
                filter(race7==1)
# IID %in% IIDs$x,
### Batch effect

batch_effect <- fread("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/batch_index.txt",col.names = c("IID","batch"))
ROSmaster <- merge(ROSmaster,batch_effect,by.x="IID",by.y="IID")
table(ROSmaster$batch)

path_IIDs <- "/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/"
IIDs <- read.csv(paste0(path_IIDs,"IIDs.csv"))

path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/Resid_PRS/",GenoType[1])
prs <- readRDS(paste0(path,"/Residual_results_all_p-vals.rds"))
prs_noapoe <- list()
for(alpha in names(prs)){
                prs_noapoe[[alpha]] <- prs[[alpha]]$residuals$IGAP_LOAD        
}
prs_noapoe <- do.call(cbind.data.frame, prs_noapoe) %>%
                rename_all(funs(paste0(.,"_noapoe"))) %>%
                mutate(mean_noapoe = rowMeans(.),
                       IID = prs[[alpha]]$residuals$IID)

path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/Resid_PRS/",GenoType[2])
prs <- readRDS(paste0(path,"/Residual_results_all_p-vals.rds"))
prs_nomhc <- list()
for(alpha in names(prs)){
                prs_nomhc[[alpha]] <- prs[[alpha]]$residuals$IGAP_LOAD        
}
prs_nomhc <- do.call(cbind.data.frame, prs_nomhc) %>%
                rename_all(funs(paste0(.,"_nomhc"))) %>%
                mutate(mean_nomhc = rowMeans(.),
                       IID = prs[[alpha]]$residuals$IID)

path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/Resid_PRS/",GenoType[3])
prs <- readRDS(paste0(path,"/Residual_results_all_p-vals.rds"))
prs_nomhcapoe <- list()
for(alpha in names(prs)){
                prs_nomhcapoe[[alpha]] <- prs[[alpha]]$residuals$IGAP_LOAD        
}
prs_nomhcapoe <- do.call(cbind.data.frame, prs_nomhcapoe)  %>%
                rename_all(funs(paste0(.,"_nomhcapoe"))) %>%
                mutate(mean_nomhcapoe = rowMeans(.),
                       IID = prs[[alpha]]$residuals$IID)

path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/Resid_PRS/",GenoType[4])
prs <- readRDS(paste0(path,"/Residual_results_all_p-vals.rds"))
prs_all <- list()
for(alpha in names(prs)){
                prs_all[[alpha]] <- prs[[alpha]]$residuals$IGAP_LOAD        
}
prs_all <- do.call(cbind.data.frame, prs_all) %>%
                rename_all(funs(paste0(.,"_all"))) %>%
                mutate(mean_all = rowMeans(.),
                       IID = prs[[alpha]]$residuals$IID)

prs <- path <- alpha <- NULL
# Adding the PRS variables

### msex
# 0 --> female
# 1 --> male

### dcfdx - Clinical cognitive diagnosis summary
# cite:https://www.radc.rush.edu/docs/var/detail.htm?category=Clinical+Diagnosis&subcategory=Dementia&variable=dcfdx

# Coding Detail
# 1 NCI, No cognitive impairment (No impaired domains)
# 2 MCI, Mild cognitive impairment (One impaired domain) and NO other cause of CI
# 3 MCI, Mild cognitive impairment (One impaired domain) AND another cause of CI
# 4 AD, Alzheimer's disease and NO other cause of CI (NINCDS PROB AD)
# = cogdx
# *5 AD, Alzheimer's disease AND another
# *6 Other dementia. Other primary cause
# variable = demlewy
# Q: Dementia due to Lewy Body disease/Variant?
#                 Coding Detail *1 Yes 2 No
# blank N/A
# cause of CI (NINCDS POSS AD) of dementia

mean(ROSmaster$age_death,na.rm = T) # 88.98304

# TABLE 1
table1 <- ROSmaster %>% mutate(LOAD = ifelse(dcfdx_lv==4|dcfdx_lv==5,"AD",
                                             ifelse(dcfdx_lv==1,"CN",
                                                    ifelse(dcfdx_lv==2|dcfdx_lv==3,"MCI","Other"))),
                               deceased = ifelse(is.na(age_death),0,1)) %>%
                merge(.,prs_noapoe,by="IID") %>%
                merge(.,prs_nomhc,by="IID") %>%
                merge(.,prs_nomhcapoe,by="IID") %>%
                merge(.,prs_all,by="IID") %>%
                # filter(!is.na(cogdx_new)) %>%
                group_by(deceased,msex,LOAD) %>% 
                summarise(n_sex=length(LOAD),
                          mean_age_death=mean(age_death,na.rm=T),
                          sd_age_death=sd(age_death, na.rm = T),
                          n_age_death = n(),
                          mean_age_bl = mean(age_bl),
                          sd_age_bl = sd(age_bl),
                          mean_education = mean(educ,na.rm=T),
                          sd_education = sd(educ,na.rm=T),
                          PRS_noapoe_mean_all = mean(mean_noapoe),
                          PRS_noapoe_sd_all = sd(mean_noapoe),
                          PRS_noapoe_mean_5e08 = mean(`5e-08_noapoe`),
                          PRS_noapoe_sd_5e08 = sd(`5e-08_noapoe`),
                          PRS_nomhc_mean_all = mean(mean_nomhc),
                          PRS_nomhc_sd_all = sd(mean_nomhc),
                          PRS_nomhc_mean_5e08 = mean(`5e-08_nomhc`),
                          PRS_nomhc_sd_5e08 = sd(`5e-08_nomhc`),
                          PRS_nomhcapoe_mean_all = mean(mean_nomhcapoe),
                          PRS_nomhcapoe_sd_all = sd(mean_nomhcapoe),
                          PRS_nomhcapoe_mean_5e08 = mean(`5e-08_nomhcapoe`),
                          PRS_nomhcapoe_sd_5e08 = sd(`5e-08_nomhcapoe`),
                          PRS_all_mean_all = mean(mean_all),
                          PRS_all_sd_all = sd(mean_all),
                          PRS_all_mean_5e08 = mean(`5e-08_all`),
                          PRS_all_sd_5e08 = sd(`5e-08_all`)) 
# mutate(se = sd / sqrt(n),
#        lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
#        upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se)

write.csv(table1,"/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Thesis/Manuscript/Figures/Table1.csv",row.names = F)


## 1.1) patho/cog variables distribution ----

# Reading ROS/MAP phenotype dataset
ROSmaster <- readRDS("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/ROSMAP_Phenotype/ROSmaster.rds")
indepvec.pathology <- c("amyloid_sqrt","tangles_sqrt")
indepvec.cognition <- c("cogn_global_random_slope","cogn_globaln_lv")
dat <- ROSmaster %>% dplyr::select(indepvec.pathology,indepvec.cognition)

g1 <- ggplot(dat, aes(x=amyloid_sqrt)) + 
                geom_histogram(color="#62906D", 
                               fill="#62906D") +
                theme_classic() +
                theme(plot.title = element_text(hjust = 0.5)) +
                xlab(TeX("$A\\beta$ levels")) +
                ylab("Frequency")

g2 <- ggplot(dat, aes(x=tangles_sqrt)) + 
                geom_histogram(color="#8BBD8C", 
                               fill="#8BBD8C") +
                theme_classic() +
                theme(plot.title = element_text(hjust = 0.5)) +
                xlab(TeX("PHFtau neuropathologies")) +
                ylab("Frequency")

g3 <- ggplot(dat, aes(x=cogn_global_random_slope)) + 
                geom_histogram(color="#7a6c25", 
                               fill="#7a6c25") +
                theme_classic() +
                theme(plot.title = element_text(hjust = 0.5)) +
                xlab("Random slope of global function") +
                ylab("Frequency") 

g4 <- ggplot(dat, aes(x=cogn_globaln_lv)) + 
                geom_histogram(color="#B09E76", 
                               fill="#B09E76") +
                theme_classic() +
                theme(plot.title = element_text(hjust = 0.5)) +
                xlab("Random slope of global function (last visit)") +
                ylab("Frequency") 

p <- plot_grid(g1, g2, g3, g4, labels=c('A', 'B','C','D'),nrow = 2)
ggsave("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/distribution_aging_patho_vars.jpg",p,
       w=11,h=15, dpi=700)

# 2) GWAS QC ----
path_to_save <- "/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/PRS/"

# Reading Filtered PNUKBB manifest
meta_pheno <- read.csv("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Thesis_code/Pan_UKBB/ukbb_manifest_filtered_phenos.csv")

# GClambda_after_QC <- read.csv(paste0(path_to_save,"GClambda_ALL.txt"))
# GClambda_after_QC$X <- NULL
# names(GClambda_after_QC) <- c("GClambda","phenotype")
# GClambda_after_QC$phenotype <- gsub(".tsv.txt","",gsub("/external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/UK_BioBank_7k_pheno/PRS_PRSice/GCLambda/","",GClambda_after_QC$phenotype))
# i1 <- match(GClambda_after_QC$phenotype, meta_pheno$phenocode_annotate_lst)
# i2 <- !is.na(i1) # to take care of non matches which are NA
GClambda_after_QC <- meta_pheno %>%
                mutate(phenotype_new =  new_pheno_annot ,
                       GC_lambda_before =  lambda_gc_EUR.y ,
                       hertiability =  saige_heritability_EUR ,
                       cases_n =  n_cases_EUR ,
                       control_n =  n_controls_EUR ,
                       remove =  to_remove ,
                       phenotype_qc_EUR = phenotype_qc_EUR,
                       phenotype = gsub("pr ","",
                                         gsub("cabg","",
                                              gsub("medadj","",
                                                   tolower(str_sub(gsub("_"," ",
                                                                        sub("\\h0..*", "", 
                                                                            phenotype_new)),4))))),
                       case_control = ifelse(is.na(control_n),"F","T"),
                       n = ifelse(is.na(control_n),cases_n,cases_n+control_n)) %>% 
                filter(remove!=1,
                       !is.na(phenotype_qc_EUR))

a <- GClambda_after_QC %>% filter(case_control=='T',
                                  phenotype_qc_EUR!="")
print(dim(a)[1])
g1 <- ggplot(a,
       aes(x=GC_lambda_before,
           y=hertiability,
           color=n,
           shape=phenotype_qc_EUR)) +
                geom_point(size=1.5) +
                theme_bw() +
                ylab(latex2exp::TeX("Hertiability Ratio ($\\H^2$)")) +
                xlab(latex2exp::TeX("Genomic Inflation Factor ($\\lambda_{GC}$)")) +
                ggtitle(latex2exp::TeX("Relationship between $\\H^2$ and $\\lambda_{GC}$ for categorical traits (n=2,031)")) +
                theme(plot.title = element_text(hjust = 0.5)) +
                scale_shape_discrete(name="Phenotype QC") +
                scale_y_continuous(breaks = c(0.05,
                                              0.1,0.15,
                                              0.20,0.25,
                                              0.30,0.35,
                                              0.40,0.45,
                                              0.50,0.55)) +
                scale_color_gradient(name="Cases and control",
                                     high="cornflowerblue",low="brown1",
                                     breaks = c(min(a$n),
                                                # median(GClambda_after_QC$n)/2,
                                                max(a$n)),
                                     labels = c(signif(min(a$n),2),
                                                # median(GClambda_after_QC$n),
                                                signif(max(a$n),2))) +
                geom_label_repel(aes(GC_lambda_before,
                                     label = phenotype),
                                 # data = a %>% filter(phenotype_qc_EUR=="PASS"),
                                 data = a %>% filter(phenotype_qc_EUR=="PASS") %>% 
                                                     filter(GC_lambda_before==max(.$GC_lambda_before)|
                                                                     hertiability  == max(.$hertiability)),
                                 fill=NA,
                                 color="black", 
                                 max.overlaps = Inf)

b <- GClambda_after_QC %>% filter(case_control=='F',
                                  phenotype_qc_EUR!="")
print(dim(b)[1])
g2 <- ggplot(b,
             aes(x=GC_lambda_before,
                 y=hertiability,
                 color=n,
                 shape=phenotype_qc_EUR)) +
                geom_point(size=1.5) +
                theme_bw() +
                ylab(latex2exp::TeX("Hertiability Ratio ($\\H^2$)")) +
                xlab(latex2exp::TeX("Genomic Inflation Factor ($\\lambda_{GC}$)")) +
                ggtitle(latex2exp::TeX("Relationship between $\\H^2$ and $\\lambda_{GC}$ for continuous traits (n=205)")) +
                theme(plot.title = element_text(hjust = 0.5)) +
                scale_shape_discrete(name="Phenotype QC") +
                scale_y_continuous(breaks = c(0.05,
                                              0.1,0.15,
                                              0.20,0.25,
                                              0.30,0.35,
                                              0.40,0.45,
                                              0.50,0.55)) +
                scale_color_gradient(name="Cases",
                                     high="cornflowerblue",low="brown1",
                                     breaks = c(min(b$n),
                                                # median(GClambda_after_QC$n)/2,
                                                max(b$n)),
                                     labels = c(signif(min(b$n),2),
                                                # median(GClambda_after_QC$n),
                                                signif(max(b$n),2))) +
                geom_label_repel(aes(GC_lambda_before,
                                     label = phenotype),
                                 data = b %>% filter(phenotype_qc_EUR=="PASS") %>% 
                                                 filter(GC_lambda_before==max(.$GC_lambda_before)|
                                                                        hertiability  == max(.$hertiability)),
                                 fill=NA,
                                 colour="black")

p <- plot_grid(g1, g2, labels=c('A', 'B'),nrow = 2)

ggsave("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/GWAS_GC_H2_compare.jpg",p,
       w=8,h=10, dpi=500)
# dat <- meta_pheno %>% dplyr::select(new_pheno_annot,
#                                     lambda_gc_EUR,
#                                     saige_heritability_EUR,
#                                     to_remove) %>%
#                 filter(to_remove!=1,
#                        !is.na(saige_heritability_EUR)) %>%
#                 dplyr::select(-to_remove)
# 
# 
# g1 <- ggplot(dat, aes(x=saige_heritability_EUR)) + 
#                 geom_histogram(color="#706f6f", 
#                                fill="#706f6f",
#                                bins=50) +
#                 theme_classic() +
#                 theme(plot.title = element_text(hjust = 0.5)) +
#                 xlab(latex2exp::TeX("Hertiability Ratio ($\\H^2$)")) +
#                 ylab("Frequency") +
#                 ggtitle("Distribution of hertiablity ratio of each phenotype") +
#                 scale_x_continuous(breaks = c(0.05,
#                                               0.1,0.15,
#                                               0.20,0.25,
#                                               0.30,0.35,
#                                               0.40,0.45,
#                                               0.50,0.55))
# 
# 
# g2 <- ggplot(dat, aes(x=lambda_gc_EUR)) + 
#                 geom_histogram(color="#706f6f", 
#                                fill="#706f6f",
#                                bins=50) +
#                 theme_classic() +
#                 theme(plot.title = element_text(hjust = 0.5)) +
#                 xlab(latex2exp::TeX("$\\lambda_{GC}$")) +
#                 ylab("Frequency") +
#                 geom_vline(xintercept = 1) +
#                 ggtitle("Distribution of genomic control of each phenotype")
# 
# p <- plot_grid(g1, g2, labels=c('A', 'B'),nrow = 2)
# pdf(paste0(path_to_save,"/distribution_GClambda_h2.pdf"),w=10,h=12)
# print(plot_grid(p, ncol=1, rel_heights=c(0.1, 1)))
# dev.off()

# 3) Principal Components plots (PCA) ----

Genotype="With_MHC_APOE"
path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/PCA/",Genotype)
pc_info_dat <- read.csv(paste0(path,"/pc_heatmap_dat.csv"))
pc_info_dat <- pc_info_dat %>% mutate(color=ifelse(pc_vector>5,1,0))
p1 <- ggplot(pc_info_dat, 
             aes(x=factor(alpha_val,levels=c("1","0.1","0.05","0.01","0.005","0.001",
                                             "5e-04","1e-04","5e-05","1e-05","5e-06",
                                             "1e-06","5e-07","1e-07","5e-08")), 
                 y=pc_vector,
                 color=as.numeric(gsub("PC","",pc_name)))) + 
                geom_quasirandom(method = "pseudorandom") +
                scale_color_gradient(name="PC number",
                                     high="cornflowerblue",low="brown1",
                                     breaks = c(min(as.numeric(gsub("PC","",pc_info_dat$pc_name))),
                                                mean(as.numeric(gsub("PC","",pc_info_dat$pc_name))),
                                                max(as.numeric(gsub("PC","",pc_info_dat$pc_name)))),
                                     labels = c(signif(min(as.numeric(gsub("PC","",pc_info_dat$pc_name))),1),
                                                signif(mean(as.numeric(gsub("PC","",pc_info_dat$pc_name))),1),
                                                signif(max(as.numeric(gsub("PC","",pc_info_dat$pc_name))),1))) +
                theme_bw() +
                theme(plot.title = element_text(hjust = 0.5)) +
                geom_hline(yintercept = 1,linetype="dashed", color = "gray60") +
                geom_hline(yintercept = 5,linetype="dashed", color = "black") +
                # geom_text(aes(x = "1", y = 1.1, 
                #               label = "1% of total variation explained"),
                #           color="red") +
                ggtitle(latex2exp::TeX("PCA results for $\\Pi_{ALL}$")) +
                labs(x = latex2exp::TeX("$\\alpha$-value threshold"), 
                     y = latex2exp::TeX("Total variation explained by each PC")) +
                geom_label_repel(aes(factor(alpha_val,levels=c("1","0.1","0.05","0.01","0.005","0.001",
                                                               "5e-04","1e-04","5e-05","1e-05","5e-06",
                                                               "1e-06","5e-07","1e-07","5e-08")),
                                     label = paste(top_pc_explain_80,"PCs")),
                                 data = pc_info_dat %>% filter(pc_name=="PC20") %>% group_by(alpha_val) %>% sample_n(1),
                                 y=3,
                                 na.rm = TRUE,
                                 angle=45,
                                 col="deepskyblue3") +
                geom_text(aes(factor(alpha_val,levels=c("1","0.1","0.05","0.01","0.005","0.001",
                                                        "5e-04","1e-04","5e-05","1e-05","5e-06",
                                                        "1e-06","5e-07","1e-07","5e-08")),
                              label = paste(top_pc_explain_80,"PCs")),
                          data = pc_info_dat %>% filter(pc_name=="PC1",alpha_val=="1e-04") ,
                          label="number of PCs explaining 80% of the total variations",
                          y=0.75+3,
                          size=4,
                          col="deepskyblue3") +
                scale_y_continuous(labels = 1:as.integer(max(pc_info_dat$pc_vector)), 
                                   breaks = 1:as.integer(max(pc_info_dat$pc_vector))) +
                geom_label_repel(data=subset(pc_info_dat, 
                                             pc_name=="PC1"),
                                 aes(label=paste0(str_sub(gsub("_"," ",sub("\\h0..*", "", PC_contributors)),4)," (",
                                                  str_extract(PC_contributors,"\\d+(\\.\\d+){0,1}%"),")")),
                                 col="black",size=3.5,
                                 fontface="bold",
                                 fill=NA)


Genotype="No_APOE"
path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/PCA/",Genotype)
pc_info_dat <- read.csv(paste0(path,"/pc_heatmap_dat.csv"))
pc_info_dat <- pc_info_dat %>% mutate(color=ifelse(pc_vector>5,1,0))
p2 <- ggplot(pc_info_dat, 
             aes(x=factor(alpha_val,levels=c("1","0.1","0.05","0.01","0.005","0.001",
                                             "5e-04","1e-04","5e-05","1e-05","5e-06",
                                             "1e-06","5e-07","1e-07","5e-08")), 
                 y=pc_vector,
                 color=as.numeric(gsub("PC","",pc_name)))) + 
                geom_quasirandom(method = "pseudorandom") +
                scale_color_gradient(name="PC number",
                                     high="cornflowerblue",low="brown1",
                                     breaks = c(min(as.numeric(gsub("PC","",pc_info_dat$pc_name))),
                                                mean(as.numeric(gsub("PC","",pc_info_dat$pc_name))),
                                                max(as.numeric(gsub("PC","",pc_info_dat$pc_name)))),
                                     labels = c(signif(min(as.numeric(gsub("PC","",pc_info_dat$pc_name))),1),
                                                signif(mean(as.numeric(gsub("PC","",pc_info_dat$pc_name))),1),
                                                signif(max(as.numeric(gsub("PC","",pc_info_dat$pc_name))),1))) +
                theme_bw() +
                theme(plot.title = element_text(hjust = 0.5)) +
                geom_hline(yintercept = 1,linetype="dashed", color = "gray60") +
                geom_hline(yintercept = 5,linetype="dashed", color = "black") +
                # geom_text(aes(x = "1", y = 1.1, 
                #               label = "1% of total variation explained"),
                #           color="red") +
                ggtitle(latex2exp::TeX("PCA results for $\\Pi_{APOE-}$")) +
                labs(x = latex2exp::TeX("$\\alpha$-value threshold"), 
                     y = latex2exp::TeX("Total variation explained by each PC")) +
                geom_label_repel(aes(factor(alpha_val,levels=c("1","0.1","0.05","0.01","0.005","0.001",
                                                               "5e-04","1e-04","5e-05","1e-05","5e-06",
                                                               "1e-06","5e-07","1e-07","5e-08")),
                                     label = paste(top_pc_explain_80,"PCs")),
                                 data = pc_info_dat %>% filter(pc_name=="PC20") %>% group_by(alpha_val) %>% sample_n(1),
                                 y=3,
                                 na.rm = TRUE,
                                 angle=45,
                                 col="deepskyblue3") +
                geom_text(aes(factor(alpha_val,levels=c("1","0.1","0.05","0.01","0.005","0.001",
                                                        "5e-04","1e-04","5e-05","1e-05","5e-06",
                                                        "1e-06","5e-07","1e-07","5e-08")),
                              label = paste(top_pc_explain_80,"PCs")),
                          data = pc_info_dat %>% filter(pc_name=="PC1",alpha_val=="1e-04") ,
                          label="number of PCs explaining 80% of the total variations",
                          y=0.75+3,
                          size=4,
                          col="deepskyblue3") +
                scale_y_continuous(labels = 1:as.integer(max(pc_info_dat$pc_vector)), 
                                   breaks = 1:as.integer(max(pc_info_dat$pc_vector))) + 
                geom_label_repel(data=subset(pc_info_dat, 
                                             pc_name=="PC1"),
                                 aes(label=paste0(str_sub(gsub("_"," ",sub("\\h0..*", "", PC_contributors)),4)," (",
                                                  str_extract(PC_contributors,"\\d+(\\.\\d+){0,1}%"),")")),
                                 col="black",size=3.5,
                                 fontface="bold",
                                 fill=NA)

Genotype="No_MHC"
path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/PCA/",Genotype)
pc_info_dat <- read.csv(paste0(path,"/pc_heatmap_dat.csv"))
pc_info_dat <- pc_info_dat %>% mutate(color=ifelse(pc_vector>5,1,0))
p3 <- ggplot(pc_info_dat, 
             aes(x=factor(alpha_val,levels=c("1","0.1","0.05","0.01","0.005","0.001",
                                             "5e-04","1e-04","5e-05","1e-05","5e-06",
                                             "1e-06","5e-07","1e-07","5e-08")), 
                 y=pc_vector,
                 color=as.numeric(gsub("PC","",pc_name)))) + 
                geom_quasirandom(method = "pseudorandom") +
                scale_color_gradient(name="PC number",
                                     high="cornflowerblue",low="brown1",
                                     breaks = c(min(as.numeric(gsub("PC","",pc_info_dat$pc_name))),
                                                mean(as.numeric(gsub("PC","",pc_info_dat$pc_name))),
                                                max(as.numeric(gsub("PC","",pc_info_dat$pc_name)))),
                                     labels = c(signif(min(as.numeric(gsub("PC","",pc_info_dat$pc_name))),1),
                                                signif(mean(as.numeric(gsub("PC","",pc_info_dat$pc_name))),1),
                                                signif(max(as.numeric(gsub("PC","",pc_info_dat$pc_name))),1))) +
                theme_bw() +
                theme(plot.title = element_text(hjust = 0.5)) +
                geom_hline(yintercept = 1,linetype="dashed", color = "gray60") +
                geom_hline(yintercept = 5,linetype="dashed", color = "black") +
                # geom_text(aes(x = "1", y = 1.1, 
                #               label = "1% of total variation explained"),
                #           color="red") +
                ggtitle(latex2exp::TeX("PCA results for $\\Pi_{MHC-}$")) +
                labs(x = latex2exp::TeX("$\\alpha$-value threshold"), 
                     y = latex2exp::TeX("Total variation explained by each PC")) +
                geom_label_repel(aes(factor(alpha_val,levels=c("1","0.1","0.05","0.01","0.005","0.001",
                                                               "5e-04","1e-04","5e-05","1e-05","5e-06",
                                                               "1e-06","5e-07","1e-07","5e-08")),
                                     label = paste(top_pc_explain_80,"PCs")),
                                 data = pc_info_dat %>% filter(pc_name=="PC20") %>% group_by(alpha_val) %>% sample_n(1),
                                 y=2,
                                 na.rm = TRUE,
                                 angle=45,
                                 col="deepskyblue3") +
                geom_text(aes(factor(alpha_val,levels=c("1","0.1","0.05","0.01","0.005","0.001",
                                                        "5e-04","1e-04","5e-05","1e-05","5e-06",
                                                        "1e-06","5e-07","1e-07","5e-08")),
                              label = paste(top_pc_explain_80,"PCs")),
                          data = pc_info_dat %>% filter(pc_name=="PC1",alpha_val=="1e-04") ,
                          label="number of PCs explaining 80% of the total variations",
                          y=0.75+2,
                          size=4,
                          col="deepskyblue3") +
                scale_y_continuous(labels = 1:as.integer(max(pc_info_dat$pc_vector)), 
                                   breaks = 1:as.integer(max(pc_info_dat$pc_vector))) +
                geom_label_repel(data=subset(pc_info_dat, 
                                             pc_name=="PC1"),
                                 aes(label=paste0(gsub("PR ","",str_sub(gsub("_"," ",sub("\\h0..*", "", PC_contributors)),4))," (",
                                                  str_extract(PC_contributors,"\\d+(\\.\\d+){0,1}%"),")")),
                                 col="black",size=3.5,
                                 fontface="bold",
                                 fill=NA)

Genotype="No_MHC_APOE"
path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/PCA/",Genotype)
pc_info_dat <- read.csv(paste0(path,"/pc_heatmap_dat.csv"))
pc_info_dat <- pc_info_dat %>% mutate(color=ifelse(pc_vector>5,1,0))
p4 <- ggplot(pc_info_dat, 
             aes(x=factor(alpha_val,levels=c("1","0.1","0.05","0.01","0.005","0.001",
                                             "5e-04","1e-04","5e-05","1e-05","5e-06",
                                             "1e-06","5e-07","1e-07","5e-08")), 
                 y=pc_vector,
                 color=as.numeric(gsub("PC","",pc_name)))) + 
                geom_quasirandom(method = "pseudorandom") +
                scale_color_gradient(name="PC number",
                                     high="cornflowerblue",low="brown1",
                                     breaks = c(min(as.numeric(gsub("PC","",pc_info_dat$pc_name))),
                                                mean(as.numeric(gsub("PC","",pc_info_dat$pc_name))),
                                                max(as.numeric(gsub("PC","",pc_info_dat$pc_name)))),
                                     labels = c(signif(min(as.numeric(gsub("PC","",pc_info_dat$pc_name))),1),
                                                signif(mean(as.numeric(gsub("PC","",pc_info_dat$pc_name))),1),
                                                signif(max(as.numeric(gsub("PC","",pc_info_dat$pc_name))),1))) +
                theme_bw() +
                theme(plot.title = element_text(hjust = 0.5)) +
                geom_hline(yintercept = 1,linetype="dashed", color = "gray60") +
                geom_hline(yintercept = 5,linetype="dashed", color = "black") +
                # geom_text(aes(x = "1", y = 1.1, 
                #               label = "1% of total variation explained"),
                #           color="red") +
                ggtitle(latex2exp::TeX("PCA results for $\\Pi_{MHC-APOE-}$")) +
                labs(x = latex2exp::TeX("$\\alpha$-value threshold"), 
                     y = latex2exp::TeX("Total variation explained by each PC")) +
                geom_label_repel(aes(factor(alpha_val,levels=c("1","0.1","0.05","0.01","0.005","0.001",
                                                               "5e-04","1e-04","5e-05","1e-05","5e-06",
                                                               "1e-06","5e-07","1e-07","5e-08")),
                                     label = paste(top_pc_explain_80,"PCs")),
                                 data = pc_info_dat %>% filter(pc_name=="PC20") %>% group_by(alpha_val) %>% sample_n(1),
                                 y=2,
                                 na.rm = TRUE,
                                 angle=45,
                                 col="deepskyblue3") +
                geom_text(aes(factor(alpha_val,levels=c("1","0.1","0.05","0.01","0.005","0.001",
                                                        "5e-04","1e-04","5e-05","1e-05","5e-06",
                                                        "1e-06","5e-07","1e-07","5e-08")),
                              label = paste(top_pc_explain_80,"PCs")),
                          data = pc_info_dat %>% filter(pc_name=="PC1",alpha_val=="1e-04") ,
                          label="number of PCs explaining 80% of the total variations",
                          y=0.75+2,
                          size=4,
                          col="deepskyblue3") +
                scale_y_continuous(labels = 1:as.integer(max(pc_info_dat$pc_vector)), 
                                   breaks = 1:as.integer(max(pc_info_dat$pc_vector))) +
                geom_label_repel(data=subset(pc_info_dat, 
                                             pc_name=="PC1"),
                                 aes(label=paste0(gsub("PR ","",str_sub(gsub("_"," ",sub("\\h0..*", "", PC_contributors)),4))," (",
                                                  str_extract(PC_contributors,"\\d+(\\.\\d+){0,1}%"),")")),
                                 col="black",size=3.5,
                                 fontface="bold",
                                 fill=NA)

ggsave("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/PCA/PCA_results_ALL.jpg",p1,
       w=12,h=10, dpi=700)
ggsave("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/PCA/PCA_results_APOE-.jpg",p2,
       w=12,h=10, dpi=700)
ggsave("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/PCA/PCA_results_MHC-.jpg",p3,
       w=12,h=10, dpi=700)
ggsave("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/PCA/PCA_results_MHC-APOE.jpg",p4,
       w=12,h=10, dpi=700)
## 3.1) Proportion of PCs out of total PRSs ----

# Reading all PC results and combining them

Genotype="With_MHC_APOE"
path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/PCA/",Genotype)
pc_info_dat <- read.csv(paste0(path,"/pc_heatmap_dat.csv"))
pc_info_dat_ALL <- pc_info_dat %>% 
                mutate(color=ifelse(pc_vector>5,1,0),
                       genome="ALL")

Genotype="No_MHC"
path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/PCA/",Genotype)
pc_info_dat <- read.csv(paste0(path,"/pc_heatmap_dat.csv"))
pc_info_dat_MHC_ <- pc_info_dat %>% 
                mutate(color=ifelse(pc_vector>5,1,0),
                       genome="MHC-")

Genotype="No_APOE"
path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/PCA/",Genotype)
pc_info_dat <- read.csv(paste0(path,"/pc_heatmap_dat.csv"))
pc_info_dat_APOE_ <- pc_info_dat %>% 
                mutate(color=ifelse(pc_vector>5,1,0),
                       genome="APOE-")

Genotype="No_MHC_APOE"
path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/PCA/",Genotype)
pc_info_dat <- read.csv(paste0(path,"/pc_heatmap_dat.csv"))
pc_info_dat_MHC_APOE_ <- pc_info_dat %>% 
                mutate(color=ifelse(pc_vector>5,1,0),
                       genome="MHC-APOE-")
pc_info_dat <- rbind(pc_info_dat_ALL,pc_info_dat_MHC_,
                     pc_info_dat_APOE_,pc_info_dat_MHC_APOE_)
                
# Reading the dimension of each PRS matrix

Genotype="With_MHC_APOE"
path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/Resid_PRS/",Genotype)
prs_info_dat <- read.csv(paste0(path,"/details.csv"))
prs_info_dat_ALL <- prs_info_dat %>% 
                mutate(genome="ALL")

Genotype="No_MHC"
path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/Resid_PRS/",Genotype)
prs_info_dat <- read.csv(paste0(path,"/details.csv"))
prs_info_dat_MHC_ <- prs_info_dat %>% 
                mutate(genome="MHC-")

Genotype="No_APOE"
path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/Resid_PRS/",Genotype)
prs_info_dat <- read.csv(paste0(path,"/details.csv"))
prs_info_dat_APOE_ <- prs_info_dat %>% 
                mutate(genome="APOE-")

Genotype="No_MHC_APOE"
path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/Resid_PRS/",Genotype)
prs_info_dat <- read.csv(paste0(path,"/details.csv"))
prs_info_dat_MHC_APOE_ <- prs_info_dat %>% 
                mutate(genome="MHC-APOE-")

prs_info_dat <- rbind(prs_info_dat_ALL,prs_info_dat_MHC_,
                      prs_info_dat_APOE_,prs_info_dat_MHC_APOE_) %>%
                mutate(alpha_val=P.value.threshhold,
                       n_prs=Number.of.PRSs.after.QC) %>%
                dplyr::select(alpha_val,
                              n_prs,
                              genome) 

pc_info_dat <- merge(pc_info_dat,prs_info_dat,by=c("genome","alpha_val"))
 
pc_info_dat <- pc_info_dat %>% 
                group_by(genome,alpha_val) %>%
                filter(pc_name %in% paste0("PC",1:top_pc_explain_80)) %>%
                mutate(pc_prop=pc_vector/n_prs) 

p1 <- ggplot(pc_info_dat,aes(y=pc_prop,
                       x=factor(alpha_val,levels=c("1","0.1","0.05","0.01","0.005","0.001",
                                                   "5e-04","1e-04","5e-05","1e-05","5e-06",
                                                   "1e-06","5e-07","1e-07","5e-08"))))+
                geom_quasirandom(alpha=.3,varwidth=TRUE,
                                 color="darkgrey",
                                 size=1.5) +
                # geom_quasirandom(method = "pseudorandom",
                #                  color="lightskyblue3",
                #                  size=1.5)+
                facet_wrap(~genome,ncol=1, scales = "free_x") +
                theme_bw() +
                xlab(latex2exp::TeX("$\\alpha$-value threshold")) +
                ylab(latex2exp::TeX("Proportion of PCs out of number of included PRSs")) +
                ggtitle("Number of PCs required to explain 80% of variance, divided by the number of PRSs in that matrix") +
                theme(plot.title = element_text(hjust = 0.5))
ggsave("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/PCA/pca_proportions_top_80.jpg",p1,
       w=10,h=9, dpi=400)                

## 3.2) PC examples plots and their contributions ----

# ALL:
## PC1, alpha<5e-08
pc_cont_dat <- pc_info_dat %>% 
                filter(genome=="ALL",
                       alpha_val==5e-08,
                       pc_name=="PC1")

contributors <- unlist(strsplit(pc_cont_dat$PC_contributors, "\n"))
pheno <- gsub("pr ","",gsub("cabg","",
                         gsub("medadj","",
                              tolower(str_sub(gsub("_"," ",sub("\\h0..*", "", contributors)),4)))))
conts <- as.numeric(gsub("%","",as.vector(sub(".*Contribution: ", "", contributors))))
pc_cont_dat <- data.frame(pheno,conts) %>%
                rbind(.,c("other",100-sum(.$conts))) %>%
                mutate(conts=signif(as.numeric(conts),2),
                       cont=ceiling(as.numeric(conts)),
                       pheno=paste0(pheno," (",conts,"%)"))

p1 <- ggplot(pc_cont_dat) + 
                geom_parliament(aes(seats = 10*cont, # -->% of contriutions (integer)
                                    fill = pheno), # --> phenotype
                                color = "white") +
                scale_fill_manual(name="PRSs",
                                  values = get_palette(palette = "Set2", 21),
                                  labels = pc_cont_dat$pheno) +
                coord_fixed() + 
                theme_void() +
                theme(title = element_text(size = 18),
                       plot.title = element_text(hjust = 0.5),
                       plot.subtitle = element_text(hjust = 0.5),
                       plot.caption = element_text(vjust = -3,hjust = 0.9),    
                       legend.position = 'bottom',
                       legend.direction = "horizontal",
                       legend.spacing.y = unit(0.05,"cm"),
                       legend.spacing.x = unit(0.05,"cm"),
                       legend.key.size = unit(0.95, 'lines'),
                       legend.text = element_text(margin = margin(r = 1, unit = 'cm'),
                                                  size=10, face="bold"),
                       legend.text.align = 0)+
                annotate("text", x = 0, y = 0.4, label = latex2exp::TeX("PC1 from $\\Pi_{ALL,5e-08}$"),colour = "black",size=6) +
                guides(fill=guide_legend(nrow=14,byrow=TRUE,reverse = FALSE,title=NULL))

# PC28, alpha<5e-08, ALL
pc_cont_dat <- pc_info_dat %>% 
                filter(genome=="ALL",
                       alpha_val==5e-08,
                       grepl("LOAD",PC_contributors)) %>%
                filter(pc_name==paste0("PC",min(as.numeric(substring(pc_name,3)))))

contributors <- unlist(strsplit(pc_cont_dat$PC_contributors, "\n"))
pheno <- gsub("pr ","",gsub("cabg","",
                            gsub("medadj","",
                                 tolower(str_sub(gsub("_"," ",sub("\\h0..*", "", contributors)),4)))))
pheno[which(grepl("LOAD",contributors))] <- "LOAD"
conts <- as.numeric(gsub("%","",as.vector(sub(".*Contribution: ", "", contributors))))
pc_cont_dat <- data.frame(pheno,conts) %>%
                rbind(.,c("other",100-sum(.$conts))) %>%
                mutate(conts=signif(as.numeric(conts),2),
                       cont=ceiling(as.numeric(conts)),
                       pheno=paste0(pheno," (",conts,"%)"))

p2 <- ggplot(pc_cont_dat) + 
                geom_parliament(aes(seats = 10*cont, # -->% of contriutions (integer)
                                    fill = pheno), # --> phenotype
                                color = "white") +
                scale_fill_manual(name="PRSs",
                                  values = get_palette(palette = "Set2", 21),
                                  labels = pc_cont_dat$pheno) +
                coord_fixed() + 
                theme_void() +
                theme(title = element_text(size = 18),
                      plot.title = element_text(hjust = 0.5),
                      plot.subtitle = element_text(hjust = 0.5),
                      plot.caption = element_text(vjust = -3,hjust = 0.9),    
                      legend.position = 'bottom',
                      legend.direction = "horizontal",
                      legend.spacing.y = unit(0.05,"cm"),
                      legend.spacing.x = unit(0.05,"cm"),
                      legend.key.size = unit(0.95, 'lines'),
                      legend.text = element_text(margin = margin(r = 1, unit = 'cm'),
                                                 size=10, face="bold"),
                      legend.text.align = 0)+
                annotate("text", x = 0, y = 0.4, label = latex2exp::TeX("PC28 from $\\Pi_{ALL,5e-08}$"),colour = "black",size=6) +
                guides(fill=guide_legend(nrow=14,byrow=TRUE,reverse = FALSE,title=NULL))


## PC1 0.001, ALL
pc_cont_dat <- pc_info_dat %>% 
                filter(genome=="ALL",
                       alpha_val==0.001,
                       pc_name=="PC1")

contributors <- unlist(strsplit(pc_cont_dat$PC_contributors, "\n"))
pheno <- gsub("pr ","",gsub("cabg","",
                            gsub("medadj","",
                                 tolower(str_sub(gsub("_"," ",sub("\\h0..*", "", contributors)),4)))))
conts <- as.numeric(gsub("%","",as.vector(sub(".*Contribution: ", "", contributors))))
pc_cont_dat <- data.frame(pheno,conts) %>%
                rbind(.,c("other",100-sum(.$conts))) %>%
                mutate(conts=signif(as.numeric(conts),2),
                       cont=ceiling(as.numeric(conts)),
                       pheno=paste0(pheno," (",conts,"%)"))

p3 <- ggplot(pc_cont_dat) + 
                geom_parliament(aes(seats = 10*cont, # -->% of contriutions (integer)
                                    fill = pheno), # --> phenotype
                                color = "white") +
                scale_fill_manual(name="PRSs",
                                  values = get_palette(palette = "Set2", 21),
                                  labels = pc_cont_dat$pheno) +
                coord_fixed() + 
                theme_void() +
                theme(title = element_text(size = 18),
                      plot.title = element_text(hjust = 0.5),
                      plot.subtitle = element_text(hjust = 0.5),
                      plot.caption = element_text(vjust = -3,hjust = 0.9),    
                      legend.position = 'bottom',
                      legend.direction = "horizontal",
                      legend.spacing.y = unit(0.05,"cm"),
                      legend.spacing.x = unit(0.05,"cm"),
                      legend.key.size = unit(0.95, 'lines'),
                      legend.text = element_text(margin = margin(r = 1, unit = 'cm'),
                                                 size=10, face="bold"),
                      legend.text.align = 0)+
                annotate("text", x = 0, y = 0.4, label = latex2exp::TeX("PC1 from $\\Pi_{ALL,0.001}$"),colour = "black",size=6) +
                guides(fill=guide_legend(nrow=14,byrow=TRUE,reverse = FALSE,title=NULL))

## PC1 0.001, MHC-
pc_cont_dat <- pc_info_dat %>% 
                filter(genome=="MHC-",
                       alpha_val==0.001,
                       pc_name=="PC1")

contributors <- unlist(strsplit(pc_cont_dat$PC_contributors, "\n"))
pheno <- gsub("pr ","",gsub("cabg","",
                            gsub("medadj","",
                                 tolower(str_sub(gsub("_"," ",sub("\\h0..*", "", contributors)),4)))))
conts <- as.numeric(gsub("%","",as.vector(sub(".*Contribution: ", "", contributors))))
pc_cont_dat <- data.frame(pheno,conts) %>%
                rbind(.,c("other",100-sum(.$conts))) %>%
                mutate(conts=signif(as.numeric(conts),2),
                       cont=ceiling(as.numeric(conts)),
                       pheno=paste0(pheno," (",conts,"%)"))

p4 <- ggplot(pc_cont_dat) + 
                geom_parliament(aes(seats = 10*cont, # -->% of contriutions (integer)
                                    fill = pheno), # --> phenotype
                                color = "white") +
                scale_fill_manual(name="PRSs",
                                  values = get_palette(palette = "Set2", 21),
                                  labels = pc_cont_dat$pheno) +
                coord_fixed() + 
                theme_void() +
                theme(title = element_text(size = 18),
                      plot.title = element_text(hjust = 0.5),
                      plot.subtitle = element_text(hjust = 0.5),
                      plot.caption = element_text(vjust = -3,hjust = 0.9),    
                      legend.position = 'bottom',
                      legend.direction = "horizontal",
                      legend.spacing.y = unit(0.05,"cm"),
                      legend.spacing.x = unit(0.05,"cm"),
                      legend.key.size = unit(0.95, 'lines'),
                      legend.text = element_text(margin = margin(r = 1, unit = 'cm'),
                                                 size=10, face="bold"),
                      legend.text.align = 0)+
                annotate("text", x = 0, y = 0.4, label = latex2exp::TeX("PC1 from $\\Pi_{MHC-,0.001}$"),colour = "black",size=6) +
                guides(fill=guide_legend(nrow=14,byrow=TRUE,reverse = FALSE,title=NULL))

# PC2, alpha<0.001, ALL
pc_cont_dat <- pc_info_dat %>% 
                filter(genome=="ALL",
                       alpha_val==0.001,
                       pc_name=="PC2")

contributors <- unlist(strsplit(pc_cont_dat$PC_contributors, "\n"))
pheno <- gsub("pr ","",gsub("cabg","",
                            gsub("medadj","",
                                 tolower(str_sub(gsub("_"," ",sub("\\h0..*", "", contributors)),4)))))
conts <- as.numeric(gsub("%","",as.vector(sub(".*Contribution: ", "", contributors))))
pc_cont_dat <- data.frame(pheno,conts) %>%
                rbind(.,c("other",100-sum(.$conts))) %>%
                mutate(conts=signif(as.numeric(conts),2),
                       cont=ceiling(as.numeric(conts)),
                       pheno=paste0(pheno," (",conts,"%)"))

p5 <- ggplot(pc_cont_dat) + 
                geom_parliament(aes(seats = 10*cont, # -->% of contriutions (integer)
                                    fill = pheno), # --> phenotype
                                color = "white") +
                scale_fill_manual(name="PRSs",
                                  values = get_palette(palette = "Set2", 21),
                                  labels = pc_cont_dat$pheno) +
                coord_fixed() + 
                theme_void() +
                theme(title = element_text(size = 18),
                      plot.title = element_text(hjust = 0.5),
                      plot.subtitle = element_text(hjust = 0.5),
                      plot.caption = element_text(vjust = -3,hjust = 0.9),    
                      legend.position = 'bottom',
                      legend.direction = "horizontal",
                      legend.spacing.y = unit(0.05,"cm"),
                      legend.spacing.x = unit(0.05,"cm"),
                      legend.key.size = unit(0.95, 'lines'),
                      legend.text = element_text(margin = margin(r = 1, unit = 'cm'),
                                                 size=10, face="bold"),
                      legend.text.align = 0)+
                annotate("text", x = 0, y = 0.4, label = latex2exp::TeX("PC2 from $\\Pi_{ALL,0.001}$"),colour = "black",size=6) +
                guides(fill=guide_legend(nrow=14,byrow=TRUE,reverse = FALSE,title=NULL))

# PC2, alpha<0.001, MHC-
pc_cont_dat <- pc_info_dat %>% 
                filter(genome=="MHC-",
                       alpha_val==0.001,
                       pc_name=="PC2")

contributors <- unlist(strsplit(pc_cont_dat$PC_contributors, "\n"))
pheno <- gsub("pr ","",gsub("cabg","",
                            gsub("medadj","",
                                 tolower(str_sub(gsub("_"," ",sub("\\h0..*", "", contributors)),4)))))
conts <- as.numeric(gsub("%","",as.vector(sub(".*Contribution: ", "", contributors))))
pc_cont_dat <- data.frame(pheno,conts) %>%
                rbind(.,c("other",100-sum(.$conts))) %>%
                mutate(conts=signif(as.numeric(conts),2),
                       cont=ceiling(as.numeric(conts)),
                       pheno=paste0(pheno," (",conts,"%)"))

p6 <- ggplot(pc_cont_dat) + 
                geom_parliament(aes(seats = 10*cont, # -->% of contriutions (integer)
                                    fill = pheno), # --> phenotype
                                color = "white") +
                scale_fill_manual(name="PRSs",
                                  values = get_palette(palette = "Set2", 21),
                                  labels = pc_cont_dat$pheno) +
                coord_fixed() + 
                theme_void() +
                theme(title = element_text(size = 18),
                      plot.title = element_text(hjust = 0.5),
                      plot.subtitle = element_text(hjust = 0.5),
                      plot.caption = element_text(vjust = -3,hjust = 0.9),    
                      legend.position = 'bottom',
                      legend.direction = "horizontal",
                      legend.spacing.y = unit(0.05,"cm"),
                      legend.spacing.x = unit(0.05,"cm"),
                      legend.key.size = unit(0.95, 'lines'),
                      legend.text = element_text(margin = margin(r = 1, unit = 'cm'),
                                                 size=10, face="bold"),
                      legend.text.align = 0)+
                annotate("text", x = 0, y = 0.4, label = latex2exp::TeX("PC2 from $\\Pi_{MHC-,0.001}$"),colour = "black",size=6) +
                guides(fill=guide_legend(nrow=14,byrow=TRUE,reverse = FALSE,title=NULL))

# PC9, alpha<1e-07, ALL
pc_cont_dat <- pc_info_dat %>% 
                filter(genome=="ALL",
                       alpha_val==1e-07,
                       grepl("LOAD",PC_contributors)) %>%
                filter(pc_name==paste0("PC",min(as.numeric(substring(pc_name,3)))))

contributors <- unlist(strsplit(pc_cont_dat$PC_contributors, "\n"))
pheno <- gsub("pr ","",gsub("cabg","",
                            gsub("medadj","",
                                 tolower(str_sub(gsub("_"," ",sub("\\h0..*", "", contributors)),4)))))
pheno[which(grepl("LOAD",contributors))] <- "LOAD"
conts <- as.numeric(gsub("%","",as.vector(sub(".*Contribution: ", "", contributors))))
pc_cont_dat <- data.frame(pheno,conts) %>%
                rbind(.,c("other",100-sum(.$conts))) %>%
                mutate(conts=signif(as.numeric(conts),2),
                       cont=ceiling(as.numeric(conts)),
                       pheno=paste0(pheno," (",conts,"%)"))

p7 <- ggplot(pc_cont_dat) + 
                geom_parliament(aes(seats = 10*cont, # -->% of contriutions (integer)
                                    fill = pheno), # --> phenotype
                                color = "white") +
                scale_fill_manual(name="PRSs",
                                  values = get_palette(palette = "Set2", 21),
                                  labels = pc_cont_dat$pheno) +
                coord_fixed() + 
                theme_void() +
                theme(title = element_text(size = 18),
                      plot.title = element_text(hjust = 0.5),
                      plot.subtitle = element_text(hjust = 0.5),
                      plot.caption = element_text(vjust = -3,hjust = 0.9),    
                      legend.position = 'bottom',
                      legend.direction = "horizontal",
                      legend.spacing.y = unit(0.05,"cm"),
                      legend.spacing.x = unit(0.05,"cm"),
                      legend.key.size = unit(0.95, 'lines'),
                      legend.text = element_text(margin = margin(r = 1, unit = 'cm'),
                                                 size=10, face="bold"),
                      legend.text.align = 0)+
                annotate("text", x = 0, y = 0.4, label = latex2exp::TeX("PC9 from $\\Pi_{ALL,1e-07}$"),colour = "black",size=6) +
                guides(fill=guide_legend(nrow=14,byrow=TRUE,reverse = FALSE,title=NULL))

# PC12, alpha<1e-05, ALL
pc_cont_dat <- pc_info_dat %>% 
                filter(genome=="ALL",
                       alpha_val==1e-05,
                       grepl("LOAD",PC_contributors)) %>%
                filter(pc_name==paste0("PC",min(as.numeric(substring(pc_name,3)))))

contributors <- unlist(strsplit(pc_cont_dat$PC_contributors, "\n"))
pheno <- gsub("pr ","",gsub("cabg","",
                            gsub("medadj","",
                                 tolower(str_sub(gsub("_"," ",sub("\\h0..*", "", contributors)),4)))))
pheno[which(grepl("LOAD",contributors))] <- "LOAD"
conts <- as.numeric(gsub("%","",as.vector(sub(".*Contribution: ", "", contributors))))
pc_cont_dat <- data.frame(pheno,conts) %>%
                rbind(.,c("other",100-sum(.$conts))) %>%
                mutate(conts=signif(as.numeric(conts),2),
                       cont=ceiling(as.numeric(conts)),
                       pheno=paste0(pheno," (",conts,"%)"))

p8 <- ggplot(pc_cont_dat) + 
                geom_parliament(aes(seats = 10*cont, # -->% of contriutions (integer)
                                    fill = pheno), # --> phenotype
                                color = "white") +
                scale_fill_manual(name="PRSs",
                                  values = get_palette(palette = "Set2", 21),
                                  labels = pc_cont_dat$pheno) +
                coord_fixed() + 
                theme_void() +
                theme(title = element_text(size = 18),
                      plot.title = element_text(hjust = 0.5),
                      plot.subtitle = element_text(hjust = 0.5),
                      plot.caption = element_text(vjust = -3,hjust = 0.9),    
                      legend.position = 'bottom',
                      legend.direction = "horizontal",
                      legend.spacing.y = unit(0.05,"cm"),
                      legend.spacing.x = unit(0.05,"cm"),
                      legend.key.size = unit(0.95, 'lines'),
                      legend.text = element_text(margin = margin(r = 1, unit = 'cm'),
                                                 size=10, face="bold"),
                      legend.text.align = 0)+
                annotate("text", x = 0, y = 0.4, label = latex2exp::TeX("PC12 from $\\Pi_{ALL,1e-05}$"),colour = "black",size=6) +
                guides(fill=guide_legend(nrow=14,byrow=TRUE,reverse = FALSE,title=NULL))



ggsave("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/PCA/PC1_5e-08_ALL.jpg",
       p1,
       w=11,h=7, dpi=400) 
ggsave("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/PCA/PC9_5e-08_ALL.jpg",
       p2,
       w=11,h=7, dpi=400)
ggsave("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/PCA/PC1_0.001_ALL.jpg",
       p3,
       w=11,h=7, dpi=400) 
ggsave("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/PCA/PC1_0.001_MHC-.jpg",
       p4,
       w=11,h=7, dpi=400) 
ggsave("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/PCA/PC2_0.001_ALL.jpg",
       p5,
       w=11,h=7, dpi=400) 
ggsave("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/PCA/PC2_0.001_MHC-.jpg",
       p6,
       w=11,h=7, dpi=400)
ggsave("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/PCA/PC10_1e-05_MHC-.jpg",
       p7,
       w=11,h=7, dpi=400) 
ggsave("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/PCA/PC6_1e-07_MHC-.jpg",
       p8,
       w=11,h=7, dpi=400) 

# 4) WGCNA power analysis plots ----

gg_dat <- list()
for(GenoType in c("With_MHC_APOE","No_APOE","No_MHC","No_MHC_APOE")){
                path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/",GenoType)
                sft_pwr_dat <- readRDS(paste0(path,"/powerstudy_results.rds"))
                
                dftmp1 <- do.call("rbind", sft_pwr_dat)
                dftmp2 <- dftmp1[,2]
                dftmp3 <- bind_rows(dftmp2, .id = "alpha_value")
                
                power_best <- unlist(dftmp1[,1]) %>% 
                                data.frame(alpha_value=names(.),power=.)
                
                dftmp3 <- merge(dftmp3,power_best,by="alpha_value") %>%
                                mutate(pred=ifelse(power==Power,1,0))
                
                gg_dat[[GenoType]] <- dftmp3
}

gg_dat <- bind_rows(gg_dat, .id = "Genome_area")

p <- ggplot(gg_dat,aes(x=factor(alpha_value,levels=c("1","0.1","0.05","0.01","0.005","0.001",
                                                     "0.0005","0.0001","5e-05","1e-05","5e-06",
                                                     "1e-06","5e-07","1e-07","5e-08")), 
                       y=SFT.R.sq,
                       color=mean.k.)) + 
                scale_color_gradient(name="Mean Connectivity",
                                     high="brown1",low="cornflowerblue") +
                geom_quasirandom(method = "pseudorandom") +
                geom_label_repel(data=subset(gg_dat,pred==1),
                                 aes(label=as.character(power)),
                                 direction="both",
                                 max.overlaps=Inf,
                                 size=5,
                                 min.segment.length = 0, # removes the arrow
                                 fontface="bold") +
                facet_wrap(~Genome_area,ncol=1) +
                theme_bw() +
                geom_hline(yintercept=c(0.85,0.9),
                           linetype="dashed",
                           col="red") +
                xlab(latex2exp::TeX("$\\alpha$-value threshold")) +
                ylab(latex2exp::TeX("$\\R^2$ (Signed)"))  +
                geom_vline(xintercept=seq(1.5, length(unique(dftmp3$alpha_value))-0.5, 1), 
                           linetype="dotted", colour="grey28")

ggsave("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/power_analysis_ALL.jpg",p,
       w=13,h=12, dpi=400)

##### all analysis:

# library(WGCNA)
# library(ggplot2)
# library(rms)
# library(stringr)
# library(forcats)
# library(latex2exp)
# Study="ROSMAP"
# 
# for(GenoType in c("With_MHC_APOE","No_MHC","No_APOE","No_MHC_APOE")){
#                 path_to_save <- "/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Thesis/Manuscript/Figures"
#                 path_wgcna <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/WGCNA/",GenoType)
#                 path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/Resid_PRS/",GenoType)
#                 wgcna_dat_results <- readRDS(file=paste0(path_wgcna,"/wgcnaresults_all",".rds"))
#                 prs <- readRDS(paste0(path,"/Residual_results_all_p-vals.rds"))
#                 
# 5) WGCNA module memberships/preservation ----

GenoType="With_MHC_APOE"
path_wgcna <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/WGCNA/",GenoType)
path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/Resid_PRS/",GenoType)
wgcna_dat_results <- readRDS(file=paste0(path_wgcna,"/wgcnaresults_all",".rds"))
prs <- readRDS(paste0(path,"/Residual_results_all_p-vals.rds"))

mostly_connected_phenos <- list()
for(prs in names(wgcna_dat_results)){
                print(prs)
                # print(wgcna_dat_results[[prs]]$hubPRSbarplot)
                hub <- wgcna_dat_results[[prs]]$hub
                
                # find AD scores
                net <- wgcna_dat_results[[prs]]$net
                # print(length(unique(net$colors)))
                ad_col <- net$colors[grep("LOAD",names(net$colors))][1] #  vascular_implants|
                hub$LOAD <- as.numeric(hub$module == ad_col)
                mostly_connected_phenos[[prs]]<-hub
                # print(table(net$colors))
                # print(names(net$colors[which(net$colors==net$colors[grep("Alz|dement|IGAP|Dement|_AD",names(net$colors))][1])]))
}

mostly_connected_phenos <- bind_rows(mostly_connected_phenos, .id = "alpha_value")

mostly_connected_phenos %<>%
                mutate(hub=gsub("cabg","",
                                gsub("medadj","",
                                     tolower(str_sub(gsub("_"," ",sub("\\h0..*", "", hub)),4)))),
                       hub=ifelse(grepl("covid",hub),"covid 19 ",
                                  ifelse(grepl("smoking",hub),"smoking ",
                                         ifelse(grepl("pr ",hub),str_sub(hub,4),
                                                hub))),
                       LOAD=ifelse(LOAD==1,"present","absent")) %>%
                group_by(hub) %>%
                summarise(alpha_value=alpha_value,
                          hub=hub,
                          hex=module,
                          order=length(alpha_value),
                          LOAD=LOAD) %>%
                mutate(order=ifelse(grepl("LOAD|_AD|Alz|dement|Dement",hub),20,order))
                
p1 <- ggplot(mostly_connected_phenos %>% filter(order!=1), 
       aes(y=reorder(hub,-order), 
           x=factor(alpha_value,levels=c("1","0.1","0.05","0.01","0.005","0.001",
                                         "0.0005","0.0001","5e-05","1e-05","5e-06",
                                         "1e-06","5e-07","1e-07","5e-08")), 
           color= hex)) + 
                geom_point(aes(shape=LOAD),
                           size=3) +
                geom_line(aes(group=hub),
                          color="grey28") +
                scale_shape_discrete(name=latex2exp::TeX("$\\PRS_{LOAD}$")) +
                scale_color_manual(values=c("cyan"="cyan","darkgreen"="darkgreen","orange"="orange","saddlebrown"="saddlebrown","sienna3"="sienna3",
                                            "paleturquoise"="paleturquoise","darkolivegreen"="darkolivegreen","palevioletred3","royalblue"="royalblue","brown"="brown","darkred"="darkred",
                                            "magenta"="magenta","navajowhite2"="navajowhite2","green"="green","red"="red","pink"="pink","greenyellow"="greenyellow",
                                            "midnightblue"="midnightblue","salmon"="salmon","yellow"="yellow","salmon4"="salmon4","skyblue3"="skyblue3","darkturquoise"="darkturquoise",
                                            "steelblue"="steelblue","darkgrey"="darkgrey","skyblue"="skyblue","maroon"="maroon","lavenderblush3"="lavenderblush3",
                                            "tan"="tan","lightgreen"="lightgreen","black"="black","lightyellow"="lightyellow","floralwhite"="floralwhite","thistle1"="thistle1",
                                            "plum2"="plum2","turquoise"="turquoise","blue"="blue","ivory"="ivory","coral1"="coral1","darkseagreen4"="darkseagreen4",
                                            "plum1"="plum1","purple"="purple","violet"="violet","bisque4"="bisque4","yellowgreen"="yellowgreen","lightsteelblue1"="lightsteelblue1",
                                            "brown4"="brown4","darkorange"="darkorange","grey60"="grey60","lightcyan1"="lightcyan1","honeydew1"="honeydew1","darkorange2"="darkorange2",
                                            "mediumpurple3"="mediumpurple3","darkslateblue"="darkslateblue","thistle2"="thistle2","darkmagenta"="darkmagenta","lightpink4"="lightpink4","orangered4"="orangered4",
                                            "white"="lightgrey","lightcyan"="lightcyan3")) +
                theme_bw() +
                theme(axis.text.x = element_text(angle=90),
                      axis.text.y = element_text(size=c(7,
                                                        rep(7,150)),
                                                 angle=0,
                                                 face="bold"),
                      plot.title = element_text(hjust = 0.5))+
                guides(color = "none") +
                ggtitle(latex2exp::TeX("Module preservation in $\\Pi_{ALL}")) +
                xlab(latex2exp::TeX("$\\alpha$-value threshold")) +
                ylab("")


GenoType="No_MHC"
path_wgcna <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/WGCNA/",GenoType)
path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/Resid_PRS/",GenoType)
wgcna_dat_results <- readRDS(file=paste0(path_wgcna,"/wgcnaresults_all",".rds"))
prs <- readRDS(paste0(path,"/Residual_results_all_p-vals.rds"))

mostly_connected_phenos <- list()
for(prs in names(wgcna_dat_results)){
                print(prs)
                # print(wgcna_dat_results[[prs]]$hubPRSbarplot)
                hub <- wgcna_dat_results[[prs]]$hub
                
                # find AD scores
                net <- wgcna_dat_results[[prs]]$net
                # print(length(unique(net$colors)))
                ad_col <- net$colors[grep("LOAD",names(net$colors))][1] #  vascular_implants|
                hub$LOAD <- as.numeric(hub$module == ad_col)
                mostly_connected_phenos[[prs]]<-hub
                # print(table(net$colors))
                # print(names(net$colors[which(net$colors==net$colors[grep("Alz|dement|IGAP|Dement|_AD",names(net$colors))][1])]))
}

mostly_connected_phenos <- bind_rows(mostly_connected_phenos, .id = "alpha_value")

mostly_connected_phenos %<>%
                mutate(hub=gsub("cabg","",
                                gsub("medadj","",
                                     tolower(str_sub(gsub("_"," ",sub("\\h0..*", "", hub)),4)))),
                       hub=ifelse(grepl("covid",hub),"covid 19 ",
                                  ifelse(grepl("smoking",hub),"smoking ",
                                         ifelse(grepl("pr ",hub),str_sub(hub,4),
                                                hub))),
                       LOAD=ifelse(LOAD==1,"present","absent")) %>%
                group_by(hub) %>%
                summarise(alpha_value=alpha_value,
                          hub=hub,
                          hex=module,
                          order=length(alpha_value),
                          LOAD=LOAD) %>%
                mutate(order=ifelse(grepl("LOAD|_AD|Alz|dement|Dement",hub),20,order))

p2 <- ggplot(mostly_connected_phenos %>% filter(order!=1), 
            aes(y=reorder(hub,-order), 
                x=factor(alpha_value,levels=c("1","0.1","0.05","0.01","0.005","0.001",
                                              "0.0005","0.0001","5e-05","1e-05","5e-06",
                                              "1e-06","5e-07","1e-07","5e-08")), 
                color= hex)) + 
                geom_point(aes(shape=LOAD),
                           size=3) +
                scale_color_manual(values=c("cyan"="cyan","darkgreen"="darkgreen","orange"="orange","saddlebrown"="saddlebrown","sienna3"="sienna3",
                                            "paleturquoise"="paleturquoise","darkolivegreen"="darkolivegreen","palevioletred3","royalblue"="royalblue","brown"="brown","darkred"="darkred",
                                            "magenta"="magenta","navajowhite2"="navajowhite2","green"="green","red"="red","pink"="pink","greenyellow"="greenyellow",
                                            "midnightblue"="midnightblue","salmon"="salmon","yellow"="yellow","salmon4"="salmon4","skyblue3"="skyblue3","darkturquoise"="darkturquoise",
                                            "steelblue"="steelblue","darkgrey"="darkgrey","skyblue"="skyblue","maroon"="maroon","lavenderblush3"="lavenderblush3",
                                            "tan"="tan","lightgreen"="lightgreen","black"="black","lightyellow"="lightyellow","floralwhite"="floralwhite","thistle1"="thistle1",
                                            "plum2"="plum2","turquoise"="turquoise","blue"="blue","ivory"="ivory","coral1"="coral1","darkseagreen4"="darkseagreen4",
                                            "plum1"="plum1","purple"="purple","violet"="violet","bisque4"="bisque4","yellowgreen"="yellowgreen","lightsteelblue1"="lightsteelblue1",
                                            "brown4"="brown4","darkorange"="darkorange","grey60"="grey60","lightcyan1"="lightcyan1","honeydew1"="honeydew1","darkorange2"="darkorange2",
                                            "mediumpurple3"="mediumpurple3","darkslateblue"="darkslateblue","thistle2"="thistle2","darkmagenta"="darkmagenta","lightpink4"="lightpink4","orangered4"="orangered4",
                                            "white"="lightgrey","lightcyan"="lightcyan3")) +
                geom_line(aes(group=hub),
                          color="grey28") +
                scale_shape_discrete(name=latex2exp::TeX("$\\PRS_{LOAD}$")) +
                theme_bw() +
                theme(axis.text.x = element_text(angle=90),
                      axis.text.y = element_text(size=c(7,
                                                        rep(7,150)),
                                                 angle=0,
                                                 face="bold"),
                      plot.title = element_text(hjust = 0.5))+
                guides(color = "none") +
                ggtitle(latex2exp::TeX("Module preservation in $\\Pi_{MHC-}")) +
                xlab(latex2exp::TeX("$\\alpha$-value threshold")) +
                ylab("")                

GenoType="No_APOE"
path_wgcna <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/WGCNA/",GenoType)
path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/Resid_PRS/",GenoType)
wgcna_dat_results <- readRDS(file=paste0(path_wgcna,"/wgcnaresults_all",".rds"))
prs <- readRDS(paste0(path,"/Residual_results_all_p-vals.rds"))

mostly_connected_phenos <- list()
for(prs in names(wgcna_dat_results)){
                print(prs)
                # print(wgcna_dat_results[[prs]]$hubPRSbarplot)
                hub <- wgcna_dat_results[[prs]]$hub
                
                # find AD scores
                net <- wgcna_dat_results[[prs]]$net
                # print(length(unique(net$colors)))
                ad_col <- net$colors[grep("LOAD",names(net$colors))][1] #  vascular_implants|
                hub$LOAD <- as.numeric(hub$module == ad_col)
                mostly_connected_phenos[[prs]]<-hub
                # print(table(net$colors))
                # print(names(net$colors[which(net$colors==net$colors[grep("Alz|dement|IGAP|Dement|_AD",names(net$colors))][1])]))
}

mostly_connected_phenos <- bind_rows(mostly_connected_phenos, .id = "alpha_value")

mostly_connected_phenos %<>%
                mutate(hub=gsub("cabg","",
                                gsub("medadj","",
                                     tolower(str_sub(gsub("_"," ",sub("\\h0..*", "", hub)),4)))),
                       hub=ifelse(grepl("covid",hub),"covid 19 ",
                                  ifelse(grepl("smoking",hub),"smoking ",
                                         ifelse(grepl("pr ",hub),str_sub(hub,4),
                                                hub))),
                       LOAD=ifelse(LOAD==1,"present","absent")) %>%
                group_by(hub) %>%
                summarise(alpha_value=alpha_value,
                          hub=hub,
                          hex=module,
                          order=length(alpha_value),
                          LOAD=LOAD) %>%
                mutate(order=ifelse(grepl("LOAD|_AD|Alz|dement|Dement",hub),20,order))

p3 <- ggplot(mostly_connected_phenos %>% filter(order!=1), 
             aes(y=reorder(hub,-order), 
                 x=factor(alpha_value,levels=c("1","0.1","0.05","0.01","0.005","0.001",
                                               "0.0005","0.0001","5e-05","1e-05","5e-06",
                                               "1e-06","5e-07","1e-07","5e-08")), 
                 color= hex)) + 
                geom_point(aes(shape=LOAD),
                           size=3) +
                scale_color_manual(values=c("cyan"="cyan","darkgreen"="darkgreen","orange"="orange","saddlebrown"="saddlebrown","sienna3"="sienna3",
                                            "paleturquoise"="paleturquoise","darkolivegreen"="darkolivegreen","palevioletred3","royalblue"="royalblue","brown"="brown","darkred"="darkred",
                                            "magenta"="magenta","navajowhite2"="navajowhite2","green"="green","red"="red","pink"="pink","greenyellow"="greenyellow",
                                            "midnightblue"="midnightblue","salmon"="salmon","yellow"="yellow","salmon4"="salmon4","skyblue3"="skyblue3","darkturquoise"="darkturquoise",
                                            "steelblue"="steelblue","darkgrey"="darkgrey","skyblue"="skyblue","maroon"="maroon","lavenderblush3"="lavenderblush3",
                                            "tan"="tan","lightgreen"="lightgreen","black"="black","lightyellow"="lightyellow","floralwhite"="floralwhite","thistle1"="thistle1",
                                            "plum2"="plum2","turquoise"="turquoise","blue"="blue","ivory"="ivory","coral1"="coral1","darkseagreen4"="darkseagreen4",
                                            "plum1"="plum1","purple"="purple","violet"="violet","bisque4"="bisque4","yellowgreen"="yellowgreen","lightsteelblue1"="lightsteelblue1",
                                            "brown4"="brown4","darkorange"="darkorange","grey60"="grey60","lightcyan1"="lightcyan1","honeydew1"="honeydew1","darkorange2"="darkorange2",
                                            "mediumpurple3"="mediumpurple3","darkslateblue"="darkslateblue","thistle2"="thistle2","darkmagenta"="darkmagenta","lightpink4"="lightpink4","orangered4"="orangered4",
                                            "white"="lightgrey","lightcyan"="lightcyan3")) +
                geom_line(aes(group=hub),
                          color="grey28") +
                scale_shape_discrete(name=latex2exp::TeX("$\\PRS_{LOAD}$")) +
                theme_bw() +
                theme(axis.text.x = element_text(angle=90),
                      axis.text.y = element_text(size=c(7,
                                                        rep(7,150)),
                                                 angle=0,
                                                 face="bold"),
                      plot.title = element_text(hjust = 0.5))+
                guides(color = "none") +
                ggtitle(latex2exp::TeX("Module preservation in $\\Pi_{APOE-}")) +
                xlab(latex2exp::TeX("$\\alpha$-value threshold")) +
                ylab("") 


GenoType="No_MHC_APOE"
path_wgcna <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/WGCNA/",GenoType)
path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/Resid_PRS/",GenoType)
wgcna_dat_results <- readRDS(file=paste0(path_wgcna,"/wgcnaresults_all",".rds"))
prs <- readRDS(paste0(path,"/Residual_results_all_p-vals.rds"))

mostly_connected_phenos <- list()
for(prs in names(wgcna_dat_results)){
                print(prs)
                # print(wgcna_dat_results[[prs]]$hubPRSbarplot)
                hub <- wgcna_dat_results[[prs]]$hub
                
                # find AD scores
                net <- wgcna_dat_results[[prs]]$net
                # print(length(unique(net$colors)))
                ad_col <- net$colors[grep("LOAD",names(net$colors))][1] #  vascular_implants|
                hub$LOAD <- as.numeric(hub$module == ad_col)
                mostly_connected_phenos[[prs]]<-hub
                # print(table(net$colors))
                # print(names(net$colors[which(net$colors==net$colors[grep("Alz|dement|IGAP|Dement|_AD",names(net$colors))][1])]))
}

mostly_connected_phenos <- bind_rows(mostly_connected_phenos, .id = "alpha_value")

mostly_connected_phenos %<>%
                mutate(hub=gsub("cabg","",
                                gsub("medadj","",
                                     tolower(str_sub(gsub("_"," ",sub("\\h0..*", "", hub)),4)))),
                       hub=ifelse(grepl("covid",hub),"covid 19 ",
                                  ifelse(grepl("smoking",hub),"smoking ",
                                         ifelse(grepl("pr ",hub),str_sub(hub,4),
                                                hub))),
                       LOAD=ifelse(LOAD==1,"present","absent")) %>%
                group_by(hub) %>%
                summarise(alpha_value=alpha_value,
                          hub=hub,
                          hex=module,
                          order=length(alpha_value),
                          LOAD=LOAD) %>%
                mutate(order=ifelse(grepl("LOAD|_AD|Alz|dement|Dement",hub),20,order))

p4 <- ggplot(mostly_connected_phenos %>% filter(order!=1), 
             aes(y=reorder(hub,-order), 
                 x=factor(alpha_value,levels=c("1","0.1","0.05","0.01","0.005","0.001",
                                               "0.0005","0.0001","5e-05","1e-05","5e-06",
                                               "1e-06","5e-07","1e-07","5e-08")), 
                 color= hex)) + 
                geom_point(aes(shape=LOAD),
                           size=3) +
                scale_color_manual(values=c("cyan"="cyan","darkgreen"="darkgreen","orange"="orange","saddlebrown"="saddlebrown","sienna3"="sienna3",
                                            "paleturquoise"="paleturquoise","darkolivegreen"="darkolivegreen","palevioletred3","royalblue"="royalblue","brown"="brown","darkred"="darkred",
                                            "magenta"="magenta","navajowhite2"="navajowhite2","green"="green","red"="red","pink"="pink","greenyellow"="greenyellow",
                                            "midnightblue"="midnightblue","salmon"="salmon","yellow"="yellow","salmon4"="salmon4","skyblue3"="skyblue3","darkturquoise"="darkturquoise",
                                            "steelblue"="steelblue","darkgrey"="darkgrey","skyblue"="skyblue","maroon"="maroon","lavenderblush3"="lavenderblush3",
                                            "tan"="tan","lightgreen"="lightgreen","black"="black","lightyellow"="lightyellow","floralwhite"="floralwhite","thistle1"="thistle1",
                                            "plum2"="plum2","turquoise"="turquoise","blue"="blue","ivory"="ivory","coral1"="coral1","darkseagreen4"="darkseagreen4",
                                            "plum1"="plum1","purple"="purple","violet"="violet","bisque4"="bisque4","yellowgreen"="yellowgreen","lightsteelblue1"="lightsteelblue1",
                                            "brown4"="brown4","darkorange"="darkorange","grey60"="grey60","lightcyan1"="lightcyan1","honeydew1"="honeydew1","darkorange2"="darkorange2",
                                            "mediumpurple3"="mediumpurple3","darkslateblue"="darkslateblue","thistle2"="thistle2","darkmagenta"="darkmagenta","lightpink4"="lightpink4","orangered4"="orangered4",
                                            "white"="lightgrey","lightcyan"="lightcyan3")) +
                geom_line(aes(group=hub),
                          color="grey28") +
                scale_shape_discrete(name=latex2exp::TeX("$\\PRS_{LOAD}$")) +
                theme_bw() +
                theme(axis.text.x = element_text(angle=90),
                      axis.text.y = element_text(size=c(7,
                                                        rep(7,150)),
                                                 angle=0,
                                                 face="bold"),
                      plot.title = element_text(hjust = 0.5))+
                guides(color = "none") +
                ggtitle(latex2exp::TeX("Module preservation in $\\Pi_{MHC-APOE-}")) +
                xlab(latex2exp::TeX("$\\alpha$-value threshold")) +
                ylab("") 

ggsave("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/module_preservation_ALL.jpg",
       p1,
       w=10,h=9.5, dpi=400)
ggsave("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/module_preservation_MHC-.jpg",
       p2,
       w=10,h=9.5, dpi=400)   
ggsave("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/module_preservation_APOE-.jpg",
       p3,
       w=10,h=9.5, dpi=400)   
ggsave("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/module_preservation_MHC-APOE-.jpg",
       p4,
       w=10,h=9.5, dpi=400)   
## 5.1) WGCNA module preservation testing/heatmap ----

modules <- c("weight ","asthma ",
             "high cholesterol ",
             "rheumatoid arthritis ")
selected_hubs <- list()
for(GenoType in c("No_APOE","No_MHC","No_MHC_APOE","With_MHC_APOE")){
                path_wgcna <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/WGCNA/",GenoType)
                path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/Resid_PRS/",GenoType)
                wgcna_dat_results <- readRDS(file=paste0(path_wgcna,"/wgcnaresults_all",".rds"))
                prs <- readRDS(paste0(path,"/Residual_results_all_p-vals.rds"))
                
                mostly_connected_phenos <- list()
                for(prs in names(wgcna_dat_results)){
                                print(prs)
                                # print(wgcna_dat_results[[prs]]$hubPRSbarplot)
                                hub <- wgcna_dat_results[[prs]]$hub
                                
                                # find AD scores
                                net <- wgcna_dat_results[[prs]]$net
                                # print(length(unique(net$colors)))
                                ad_col <- net$colors[grep("IGAP",names(net$colors))][1] #  vascular_implants|
                                hub$LOAD <- as.numeric(hub$module == ad_col)
                                mostly_connected_phenos[[prs]]<-hub
                                # print(table(net$colors))
                                # print(names(net$colors[which(net$colors==net$colors[grep("Alz|dement|IGAP|Dement|_AD",names(net$colors))][1])]))
                }
                
                mostly_connected_phenos <- bind_rows(mostly_connected_phenos, .id = "alpha_value")
                
                mostly_connected_phenos %<>%
                                mutate(hub=gsub("cabg","",
                                                gsub("medadj","",
                                                     tolower(str_sub(gsub("_"," ",sub("\\h0..*", "", hub)),4)))),
                                       hub=ifelse(grepl("covid",hub),"covid 19 ",
                                                  ifelse(grepl("smoking",hub),"smoking ",
                                                         ifelse(grepl("pr ",hub),str_sub(hub,4),
                                                                hub))),
                                       LOAD=ifelse(LOAD==1,"present","absent")) %>%
                                group_by(hub) %>%
                                summarise(alpha_value=alpha_value,
                                          hub=hub,
                                          hex=module,
                                          order=length(alpha_value),
                                          LOAD=LOAD) %>%
                                mutate(order=ifelse(grepl("IGAP|_AD|Alz|dement|Dement",hub),20,order))
                selected_hubs[[GenoType]] <- mostly_connected_phenos %>% filter(hub %in% modules)
}
dat <- bind_rows(selected_hubs, .id = "genotype") %>%
                mutate(genos = ifelse(genotype=="No_MHC","MHC-",
                                      ifelse(genotype=="No_APOE","APOE-",
                                             ifelse(genotype=="No_MHC_APOE","MHC-APOE-",
                                                    "ALL"))))

PC_lst <- list()
for(GenoType in unique(dat$genotype)){
                path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/Resid_PRS/",GenoType)
                prs <- readRDS(paste0(path,"/Residual_results_all_p-vals.rds"))
                
                path_wgcna <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/WGCNA/",GenoType)
                wgcna_dat_results <- readRDS(file=paste0(path_wgcna,"/wgcnaresults_all",".rds"))   
                tmp_dat <- dat %>% filter(genotype==GenoType)
                for(alpha in unique(tmp_dat$alpha_value)){
                                tmp_tmp_dat <- tmp_dat %>% filter(alpha_value==alpha)
                                for(col in unique(tmp_tmp_dat$hex)){
                                                name <- tmp_tmp_dat %>% filter(hex==col) %>% .$hub
                                                gen <- tmp_tmp_dat %>% filter(hex==col) %>% .$genos
                                                PC_lst[[paste0(gen,"_",alpha," (",col,")_",name)]] = wgcna_dat_results[[alpha]]$net$MEs[[paste0("ME",col)]]
                                }
                                
                }
}
                
module_preserve_dat <- bind_rows(PC_lst, .id = "id")

for(module in modules){
                # Weight module:
                aa <- module_preserve_dat %>% 
                                dplyr::select(grep(module,names(.))) %>% 
                                setNames(gsub(paste0("_",module),
                                              "",
                                              perl = TRUE,
                                              names(.)))
                # Creating a heatmap for all 4 (chosen) hub PRSs
                cormat <- round(cor(aa),2)
                # Reorder the correlation matrix
                cormat <- reorder_cormat(cormat)
                upper_tri <- get_upper_tri(cormat)
                # Melt the correlation matrix
                melted_cormat <- reshape2::melt(upper_tri, na.rm = TRUE)
                # Create a ggheatmap
                p <- ggplot(melted_cormat, aes(Var2, 
                                               Var1, 
                                               fill = value))+
                                geom_tile(color = "white")+
                                scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                                                     midpoint = 0, limit = c(-1,1), space = "Lab", 
                                                     name="Pearson\nCorrelation") +
                                theme_minimal()+ # minimal theme
                                theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                                                 hjust = 1),
                                      plot.title = element_text(hjust = 0.5))+
                                coord_fixed() + 
                                theme(
                                                axis.title.x = element_blank(),
                                                axis.title.y = element_blank(),
                                                panel.grid.major = element_blank(),
                                                panel.border = element_blank(),
                                                panel.background = element_blank(),
                                                axis.ticks = element_blank(),
                                                legend.justification = c(1, 0),
                                                legend.position = c(0.6, 0.7),
                                                legend.direction = "horizontal")+
                                guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                                             title.position = "top", title.hjust = 0.5)) +
                                ggtitle(paste0("Hub PRS of ",module)) +
                                geom_text(aes(Var2, Var1, label = value), 
                                          color = "black", size = 2.5)
                ggsave(paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/",
                              module,"_heatmap_preservation.jpg"),p,
                       w=10,h=10, dpi=400)
}

## 3.2) WGCNA denrogram examples
GenoType=c("No_APOE","No_MHC","No_MHC_APOE","With_MHC_APOE")[4]
path_wgcna <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/WGCNA/",GenoType)
path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/Resid_PRS/",GenoType)
wgcna_dat_results <- readRDS(file=paste0(path_wgcna,"/wgcnaresults_all",".rds"))
net <- wgcna_dat_results$`5e-06`$net

jpeg(paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/",
           "dendrogram_ALL_5e-06.jpg"),w=11,h=4,unit="in",res=400)
plotDendroAndColors(net$dendrograms[[1]],
                    net$colors,
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE,
                    hang = 0,
                    #addGuide = TRUE,
                    guideHang = 0.05,
                    main = latex2exp::TeX("Dendrogram and module colors from $\\Pi_{ALL,5e-06}$"))
dev.off()

# 6) WGCNA PC1 variation explained for each module ----

genotype_lst <- NULL
PVAL <- NULL
Module <- NULL
PC1 <- NULL
PC2 <- NULL
PC_top_80 <- NULL
N <- NULL
index <- 1
Study <- "ROSMAP"
for(GenoType in c("No_APOE","No_MHC","No_MHC_APOE","With_MHC_APOE")){ #"No_APOE","No_MHC","No_MHC_APOE",
                path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/Resid_PRS/",GenoType)
                prs <- readRDS(paste0(path,"/Residual_results_all_p-vals.rds"))
                
                path_wgcna <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/WGCNA/",GenoType)
                wgcna_dat_results <- readRDS(file=paste0(path_wgcna,"/wgcnaresults_all",".rds"))   
                for(pthres in names(prs)[order(as.numeric(names(prs)),decreasing = F)]){
                                for(color in names(table(wgcna_dat_results[[pthres]]$net$colors))){
                                                print(pthres)
                                                df <- prs[[pthres]]$residuals
                                                rownames(df) <- df$IID
                                                df$IID <- NULL
                                                prs_names <- names(wgcna_dat_results[[pthres]]$net$colors[wgcna_dat_results[[pthres]]$net$colors==color])
                                                matrix <- df %>% dplyr::select(prs_names)
                                                respca <- prcomp(matrix)
                                                pc1_cont <- sum(summary(respca)$importance[2,1])*100
                                                if(dim(matrix)[2]==1){
                                                                pc2_cont <- NA     
                                                } else{
                                                                pc2_cont <- sum(summary(respca)$importance[2,2])*100    
                                                }
                                                
                                                for(r in 1:dim(summary(respca)$importance)[2]){
                                                                if(sum(summary(respca)$importance[2,1:r])>0.80){
                                                                                top_r_80_explained <- r
                                                                                break
                                                                }
                                                }
                                                genotype_lst[index] <- GenoType
                                                PVAL[index] <- pthres
                                                Module[index] <- color
                                                PC1[index] <- pc1_cont
                                                PC2[index] <- pc2_cont
                                                PC_top_80[index] <- top_r_80_explained
                                                N[index] <- dim(matrix)[2]
                                                index <- index+1
                                                
                                }
                }
}               

wgcna_pca_dat <- as.data.frame(cbind(genotype_lst,PVAL,Module,PC1,PC2,PC_top_80,N))

fair_cols <- list("#38170B"="1",
                  "#BF1B0B"="0.1",
                  "#FFC465"="0.05",
                  "#66ADE5"="0.01",
                  "#252A52"="0.005",
                  "#38170B"="0.001",
                  "#BF1B0B"="0.0005",
                  "#FFC465"="0.0001", 
                  "#66ADE5"="5e-05",
                  "#252A52"="1e-05",
                  "#38170B"="5e-06",
                  "#BF1B0B"="1e-06",
                  "#FFC465"="5e-07",
                  "#66ADE5"="1e-07",
                  "#252A52"="5e-08")
wgcna_pca_dat %<>% 
                left_join(enframe(fair_cols) %>% unnest(value), 
                          by = c('PVAL' = 'value')) %>% 
                rename('color' = name) %>% 
                mutate(PC1=as.numeric(PC1),
                       genotype_lst=ifelse(genotype_lst=="No_APOE","APOE-",
                                           ifelse(genotype_lst=="No_MHC","MHC-",
                                                  ifelse(genotype_lst=="No_MHC_APOE","MHC-APOE-",
                                                         "ALL")))) %>% 
                filter(PVAL %in% c("1","0.1","0.05","0.01","0.005","0.001","0.0005","0.0001","5e-05","1e-05","5e-06","1e-06","5e-07","1e-07","5e-08"),
                       !Module %in% c("grey")) %>%
                mutate(ModuleSize=(as.numeric(N))/max(as.numeric(N)))  #color=adjust_transparency(color, alpha = N)


p1 <- ggplot(wgcna_pca_dat, 
             aes(x=factor(PVAL,levels=c("1","0.1","0.05","0.01","0.005","0.001",
                                        "0.0005","0.0001","5e-05","1e-05","5e-06",
                                        "1e-06","5e-07","1e-07","5e-08")), 
                 y=PC1,
                 color=color)) + 
                geom_beeswarm(aes(size=ModuleSize)) + #geom_boxplot() +
                # scale_fill_manual(values = brightness(fair_cols, 0.9)) +
                facet_wrap(~factor(genotype_lst,
                                   levels=c("ALL","MHC-","APOE-","MHC-APOE-")), scale="free_y",ncol=1) +
                theme_bw() +
                theme(plot.title = element_text(hjust = 0.5)) +
                guides(color = "none") +
                scale_size_continuous(name="Module Size",
                                      range  = c(0,5),
                                      limits = c(min(as.numeric(wgcna_pca_dat$ModuleSize)), max(as.numeric(wgcna_pca_dat$ModuleSize))), 
                                      breaks = c(min(as.numeric(wgcna_pca_dat$ModuleSize)),
                                                 mean(as.numeric(wgcna_pca_dat$ModuleSize)),
                                                 max(as.numeric(wgcna_pca_dat$ModuleSize))),
                                      labels = c(10,50,1200)) +
                # scale_colour_manual(name="ModuleSize",
                #                     values=c(min(as.numeric(wgcna_pca_dat$N)),
                #                              median(as.numeric(wgcna_pca_dat$N)),
                #                              max(as.numeric(wgcna_pca_dat$N)))) + 
                ggtitle('Total variation explained by first PC amongst all Modules') +
                labs(x = latex2exp::TeX("$\\alpha$-value threshold"), 
                     y = 'Percentage') +
                stat_summary(
                                fun.data = stat_box_data,
                                geom = "text",
                                hjust = 0.5,
                                vjust = 0.9
                )

ggsave("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/ePRS_explained.jpg",p1,
       w=13,h=8, dpi=700)

# 7) Module association results initial ----

genotype_lst <- NULL
PVAL <- NULL
Module <- NULL
hubs <- NULL
pval_assoc <- NULL
beta <- NULL
psigned_assoc <- NULL
pfdr_assoc <- NULL
pheno <- NULL
N <- NULL
index <- 1
Study <- "ROSMAP"
for(GenoType in c("No_APOE","No_MHC","No_MHC_APOE","With_MHC_APOE")){ #"No_APOE","No_MHC","No_MHC_APOE",
                path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/Resid_PRS/",GenoType)
                prs <- readRDS(paste0(path,"/Residual_results_all_p-vals.rds"))
                
                path_wgcna <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/WGCNA/",GenoType)
                wgcna_dat_results <- readRDS(file=paste0(path_wgcna,"/wgcnaresults_all",".rds"))   
                
                for(pthres in names(prs)[order(as.numeric(names(prs)),decreasing = F)]){
                                for(color in names(table(wgcna_dat_results[[pthres]]$net$colors))){
                                                if(color=="grey"){next}
                                                print(pthres)
                                                df <- prs[[pthres]]$residuals
                                                rownames(df) <- df$IID
                                                df$IID <- NULL
                                                prs_names <- names(wgcna_dat_results[[pthres]]$net$colors[wgcna_dat_results[[pthres]]$net$colors==color])
                                                matrix <- df %>% dplyr::select(prs_names)
                                                assoc <- wgcna_dat_results[[pthres]]$assoc_res %>%
                                                                filter(egene==paste0("ME",color))
                                                hub <- wgcna_dat_results[[pthres]]$hub
                                                for(phen in assoc$pheno){
                                                                tmp_dat <- assoc %>% filter(pheno==phen)
                                                                hub_name <- hub %>% filter(module==color) %>% 
                                                                                dplyr::select(hub)
                                                                
                                                                hubs[index] <- gsub("_"," ",gsub("_h0..*","",hub_name$hub)) 
                                                                genotype_lst[index] <- GenoType
                                                                PVAL[index] <- pthres
                                                                Module[index] <- color
                                                                pval_assoc[index] <- tmp_dat$p
                                                                psigned_assoc[index] <- tmp_dat$signedp
                                                                pfdr_assoc[index] <- tmp_dat$fdr
                                                                beta[index] <- tmp_dat$b
                                                                pheno[index] <- phen 
                                                                N[index] <- dim(matrix)[2]
                                                                index <- index+1
                                                }
                                                
                                                
                                                
                                                
                                                
                                                
                                                
                                }
                }
}               

# Adding if LOAD is present in each module
dat_most_connected <- list()
for(GenoType in c("No_APOE","No_MHC","No_MHC_APOE","With_MHC_APOE")){
                path_wgcna <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/WGCNA/",GenoType)
                mostly_connected_phenos <- list()
                for(prs in names(wgcna_dat_results)){
                                print(prs)
                                # print(wgcna_dat_results[[prs]]$hubPRSbarplot)
                                hub <- wgcna_dat_results[[prs]]$hub
                                
                                # find AD scores
                                net <- wgcna_dat_results[[prs]]$net
                                # print(length(unique(net$colors)))
                                ad_col <- net$colors[grep("LOAD",names(net$colors))][1] #  vascular_implants|
                                hub$LOAD <- as.numeric(hub$module == ad_col)
                                mostly_connected_phenos[[prs]]<-hub
                                # print(table(net$colors))
                                # print(names(net$colors[which(net$colors==net$colors[grep("Alz|dement|IGAP|Dement|_AD",names(net$colors))][1])]))
                }
                
                mostly_connected_phenos <- bind_rows(mostly_connected_phenos, .id = "PVAL")
                dat_most_connected[[GenoType]] <- mostly_connected_phenos
}
dat_most_connected <- bind_rows(dat_most_connected, .id = "geno")

wgcna_pca_dat <- data.frame(genotype_lst,
                            PVAL,
                            pheno,
                            Module,
                            hubs,
                            pval_assoc,
                            psigned_assoc,
                            pfdr_assoc,
                            beta,
                            N) %>%
                mutate(geno=genotype_lst,
                       genotype_lst=ifelse(genotype_lst=="No_APOE","APOE-",
                                           ifelse(genotype_lst=="No_MHC","MHC-",
                                                  ifelse(genotype_lst=="No_MHC_APOE","MHC-APOE-",
                                                         "ALL")))) %>% 
                filter(PVAL %in% c("1","0.1","0.05","0.01","0.005","0.001","0.0005","0.0001","5e-05","1e-05","5e-06","1e-06","5e-07","1e-07","5e-08"),
                       !Module %in% c("grey")) %>%
                mutate(ModuleSize=(as.numeric(N))/max(as.numeric(N)),
                       effectsize=ifelse(beta<0,"Negative","Positive"),
                       color=ifelse(pfdr_assoc<=0.05,Module,Module), # "lightgrey"
                       phenotype=pheno,
                       pheno=ifelse(pheno=="amyloid_sqrt","Total AB",
                                    ifelse(pheno=="cogdx","Final AD",
                                           ifelse(pheno=="cogn_global_random_slope","Global Slope",
                                                  ifelse(pheno=="cogn_globaln_lv","Global Last Visit",
                                                         "PHF tau")))),
                       hex = col2hex(color),
                       color_apoe_=ifelse(pval_assoc<0.05,Module,Module), # "lightgrey"
                       hubs = gsub(" PR "," ",hubs)) %>% 
                merge(.,data.frame(dat_most_connected), 
                       by.x=c("geno","PVAL","Module"),
                       by.y=c("geno","PVAL","module")) %>%
                mutate(LOAD=ifelse(LOAD==1,"present","absent"))

p1 <- ggplot(wgcna_pca_dat %>% filter(genotype_lst=="ALL"), 
             aes(x=factor(PVAL,levels=c("1","0.1","0.05","0.01","0.005","0.001",
                                        "0.0005","0.0001","5e-05","1e-05","5e-06",
                                        "1e-06","5e-07","1e-07","5e-08")), 
                 y=abs(psigned_assoc),
                 color=color,shape=LOAD)) + 
                geom_beeswarm(aes(size=ModuleSize)) + #geom_boxplot() +
                # scale_fill_manual(values = brightness(fair_cols, 0.9)) +
                facet_wrap(~pheno, scale="free_y",ncol=1) +
                theme_bw() +
                theme(plot.title = element_text(hjust = 0.5)) +
                guides(color = "none") +
                scale_color_manual(values=c("cyan"="cyan","darkgreen"="darkgreen","orange"="orange","saddlebrown"="saddlebrown","sienna3"="sienna3",
                                            "paleturquoise"="paleturquoise","darkolivegreen"="darkolivegreen","palevioletred3","royalblue"="royalblue","brown"="brown","darkred"="darkred",
                                            "magenta"="magenta","navajowhite2"="navajowhite2","green"="green","red"="red","pink"="pink","greenyellow"="greenyellow",
                                            "midnightblue"="midnightblue","salmon"="salmon","yellow"="yellow","salmon4"="salmon4","skyblue3"="skyblue3","darkturquoise"="darkturquoise",
                                            "steelblue"="steelblue","darkgrey"="darkgrey","skyblue"="skyblue","maroon"="maroon","lavenderblush3"="lavenderblush3",
                                            "tan"="tan","lightgreen"="lightgreen","black"="black","lightyellow"="lightyellow","floralwhite"="floralwhite","thistle1"="thistle1",
                                            "plum2"="plum2","turquoise"="turquoise","blue"="blue","ivory"="ivory","coral1"="coral1","darkseagreen4"="darkseagreen4",
                                            "plum1"="plum1","purple"="purple","violet"="violet","bisque4"="bisque4","yellowgreen"="yellowgreen","lightsteelblue1"="lightsteelblue1",
                                            "brown4"="brown4","darkorange"="darkorange","grey60"="grey60","lightcyan1"="lightcyan1","honeydew1"="honeydew1","darkorange2"="darkorange2",
                                            "mediumpurple3"="mediumpurple3","darkslateblue"="darkslateblue","thistle2"="thistle2","darkmagenta"="darkmagenta","lightpink4"="lightpink4","orangered4"="orangered4",
                                            "white"="lightgrey","lightcyan"="lightcyan3","lightgrey"="cornsilk3")) +
                # scale_shape(name = latex2exp::TeX("$\\beta$ sign")) +
                scale_size_continuous(name="Module Size",
                                      range  = c(2,11),
                                      limits = c(min(as.numeric(wgcna_pca_dat$ModuleSize)), max(as.numeric(wgcna_pca_dat$ModuleSize))), 
                                      breaks = c(min(as.numeric(wgcna_pca_dat$ModuleSize)),
                                                 mean(as.numeric(wgcna_pca_dat$ModuleSize)),
                                                 max(as.numeric(wgcna_pca_dat$ModuleSize))),
                                      labels = c(10,50,1200)) +
                # scale_colour_manual(name="ModuleSize",
                #                     values=c(min(as.numeric(wgcna_pca_dat$N)),
                #                              median(as.numeric(wgcna_pca_dat$N)),
                #                              max(as.numeric(wgcna_pca_dat$N)))) + 
                ggtitle(latex2exp::TeX("Phenotype ~ $\\ePRS_{\\Pi_{ALL}}\\ + Covariates")) +
                labs(x = latex2exp::TeX("$\\alpha$-value threshold"), 
                     y = latex2exp::TeX("-$\\log_{10}(p)$")) +
                geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "tomato3") +
                # geom_text(aes(x = "0.1", y = -log10(0.003), 
                #               label = "p-value = 0.05"),
                #           color="tomato3") +
                geom_label_repel(data=subset(wgcna_pca_dat %>% filter(genotype_lst=="ALL"), 
                                      pfdr_assoc < 0.05),
                                 fill = NA,
                          aes(label=paste(Module,"--",str_sub(tolower(hubs),4))),col="black",size=3,fontface="bold")

p2 <- ggplot(wgcna_pca_dat %>% filter(genotype_lst=="APOE-"), 
             aes(x=factor(PVAL,levels=c("1","0.1","0.05","0.01","0.005","0.001",
                                        "0.0005","0.0001","5e-05","1e-05","5e-06",
                                        "1e-06","5e-07","1e-07","5e-08")), 
                 y=abs(psigned_assoc),
                 color=color_apoe_,shape=LOAD)) + 
                geom_beeswarm(aes(size=ModuleSize)) + #geom_boxplot() +
                # scale_fill_manual(values = brightness(fair_cols, 0.9)) +
                facet_wrap(~pheno, scale="free_y",ncol=1) +
                theme_bw() +
                theme(plot.title = element_text(hjust = 0.5)) +
                guides(color = "none") +
                scale_color_manual(values=c("cyan"="cyan","darkgreen"="darkgreen","orange"="orange","saddlebrown"="saddlebrown","sienna3"="sienna3",
                                            "paleturquoise"="paleturquoise","darkolivegreen"="darkolivegreen","palevioletred3","royalblue"="royalblue","brown"="brown","darkred"="darkred",
                                            "magenta"="magenta","navajowhite2"="navajowhite2","green"="green","red"="red","pink"="pink","greenyellow"="greenyellow",
                                            "midnightblue"="midnightblue","salmon"="salmon","yellow"="yellow","salmon4"="salmon4","skyblue3"="skyblue3","darkturquoise"="darkturquoise",
                                            "steelblue"="steelblue","darkgrey"="darkgrey","skyblue"="skyblue","maroon"="maroon","lavenderblush3"="lavenderblush3",
                                            "tan"="tan","lightgreen"="lightgreen","black"="black","lightyellow"="lightyellow","floralwhite"="floralwhite","thistle1"="thistle1",
                                            "plum2"="plum2","turquoise"="turquoise","blue"="blue","ivory"="ivory","coral1"="coral1","darkseagreen4"="darkseagreen4",
                                            "plum1"="plum1","purple"="purple","violet"="violet","bisque4"="bisque4","yellowgreen"="yellowgreen","lightsteelblue1"="lightsteelblue1",
                                            "brown4"="brown4","darkorange"="darkorange","grey60"="grey60","lightcyan1"="lightcyan1","honeydew1"="honeydew1","darkorange2"="darkorange2",
                                            "mediumpurple3"="mediumpurple3","darkslateblue"="darkslateblue","thistle2"="thistle2","darkmagenta"="darkmagenta","lightpink4"="lightpink4","orangered4"="orangered4",
                                            "white"="lightgrey","lightcyan"="lightcyan3","lightgrey"="cornsilk3")) +
                # scale_shape(name = latex2exp::TeX("$\\beta$ sign")) +
                scale_size_continuous(name="Module Size",
                                      range  = c(2,11),
                                      limits = c(min(as.numeric(wgcna_pca_dat$ModuleSize)), max(as.numeric(wgcna_pca_dat$ModuleSize))), 
                                      breaks = c(min(as.numeric(wgcna_pca_dat$ModuleSize)),
                                                 mean(as.numeric(wgcna_pca_dat$ModuleSize)),
                                                 max(as.numeric(wgcna_pca_dat$ModuleSize))),
                                      labels = c(10,50,1200)) +
                ggtitle(latex2exp::TeX("Phenotype ~ $\\ePRS_{\\Pi_{APOE-}}\\ + Covariates")) +
                labs(x = latex2exp::TeX("$\\alpha$-value threshold"), 
                     y = latex2exp::TeX("-$\\log_{10}(p)$")) +
                geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "tomato3") +
                # geom_text(aes(x = "0.1", y = -log10(0.03), 
                #               label = "p-value = 0.05"),
                #           color="tomato3") +
                geom_label_repel(data=subset(wgcna_pca_dat %>% filter(genotype_lst=="APOE-") %>% 
                                                             group_by(genotype_lst,pheno) %>% 
                                                             mutate(pfdr_assoc = ifelse(pfdr_assoc==min(pfdr_assoc),1,0)), 
                                             pfdr_assoc == 1),
                                aes(label=paste(Module,"--",str_sub(tolower(hubs),4))),
                                col="black",
                                fill = NA,
                                size=3,fontface="bold",
                                max.overlaps = Inf)
p3 <- ggplot(wgcna_pca_dat %>% filter(genotype_lst=="MHC-"), 
             aes(x=factor(PVAL,levels=c("1","0.1","0.05","0.01","0.005","0.001",
                                        "0.0005","0.0001","5e-05","1e-05","5e-06",
                                        "1e-06","5e-07","1e-07","5e-08")), 
                 y=abs(psigned_assoc),
                 color=color,shape=LOAD)) + 
                geom_beeswarm(aes(size=ModuleSize)) + #geom_boxplot() +
                # scale_fill_manual(values = brightness(fair_cols, 0.9)) +
                facet_wrap(~pheno, scale="free_y",ncol=1) +
                theme_bw() +
                theme(plot.title = element_text(hjust = 0.5)) +
                guides(color = "none") +
                scale_color_manual(values=c("cyan"="cyan","darkgreen"="darkgreen","orange"="orange","saddlebrown"="saddlebrown","sienna3"="sienna3",
                                            "paleturquoise"="paleturquoise","darkolivegreen"="darkolivegreen","palevioletred3","royalblue"="royalblue","brown"="brown","darkred"="darkred",
                                            "magenta"="magenta","navajowhite2"="navajowhite2","green"="green","red"="red","pink"="pink","greenyellow"="greenyellow",
                                            "midnightblue"="midnightblue","salmon"="salmon","yellow"="yellow","salmon4"="salmon4","skyblue3"="skyblue3","darkturquoise"="darkturquoise",
                                            "steelblue"="steelblue","darkgrey"="darkgrey","skyblue"="skyblue","maroon"="maroon","lavenderblush3"="lavenderblush3",
                                            "tan"="tan","lightgreen"="lightgreen","black"="black","lightyellow"="lightyellow","floralwhite"="floralwhite","thistle1"="thistle1",
                                            "plum2"="plum2","turquoise"="turquoise","blue"="blue","ivory"="ivory","coral1"="coral1","darkseagreen4"="darkseagreen4",
                                            "plum1"="plum1","purple"="purple","violet"="violet","bisque4"="bisque4","yellowgreen"="yellowgreen","lightsteelblue1"="lightsteelblue1",
                                            "brown4"="brown4","darkorange"="darkorange","grey60"="grey60","lightcyan1"="lightcyan1","honeydew1"="honeydew1","darkorange2"="darkorange2",
                                            "mediumpurple3"="mediumpurple3","darkslateblue"="darkslateblue","thistle2"="thistle2","darkmagenta"="darkmagenta","lightpink4"="lightpink4","orangered4"="orangered4",
                                            "white"="lightgrey","lightcyan"="lightcyan3","lightgrey"="cornsilk3")) +
                # scale_shape(name = latex2exp::TeX("$\\beta$ sign")) +
                scale_size_continuous(name="Module Size",
                                      range  = c(2,11),
                                      limits = c(min(as.numeric(wgcna_pca_dat$ModuleSize)), max(as.numeric(wgcna_pca_dat$ModuleSize))), 
                                      breaks = c(min(as.numeric(wgcna_pca_dat$ModuleSize)),
                                                 mean(as.numeric(wgcna_pca_dat$ModuleSize)),
                                                 max(as.numeric(wgcna_pca_dat$ModuleSize))),
                                      labels = c(10,50,1200)) +
                ggtitle(latex2exp::TeX("Phenotype ~ $\\ePRS_{\\Pi_{MHC-}}\\ + Covariates")) +
                labs(x = latex2exp::TeX("$\\alpha$-value threshold"), 
                     y = latex2exp::TeX("-$\\log_{10}(p)$")) +
                geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "tomato3") +
                # geom_text(aes(x = "0.1", y = -log10(0.003), 
                #               label = "p-value = 0.05"),
                #           color="tomato3") +
                geom_label_repel(data=subset(wgcna_pca_dat %>% filter(genotype_lst=="MHC-"), 
                                            pfdr_assoc < 0.05),
                                 fill=NA,
                                aes(label=paste(Module,"--",str_sub(tolower(hubs),4))),col="black",size=3,fontface="bold",
                                max.overlaps=15)

p4 <- ggplot(wgcna_pca_dat %>% filter(genotype_lst=="MHC-APOE-"), 
             aes(x=factor(PVAL,levels=c("1","0.1","0.05","0.01","0.005","0.001",
                                        "0.0005","0.0001","5e-05","1e-05","5e-06",
                                        "1e-06","5e-07","1e-07","5e-08")), 
                 y=abs(psigned_assoc),
                 color=color_apoe_,shape=LOAD)) + 
                geom_beeswarm(aes(size=ModuleSize)) + #geom_boxplot() +
                # scale_fill_manual(values = brightness(fair_cols, 0.9)) +
                facet_wrap(~pheno, scale="free_y",ncol=1) +
                theme_bw() +
                theme(plot.title = element_text(hjust = 0.5)) +
                guides(color = "none") +
                scale_color_manual(values=c("cyan"="cyan","darkgreen"="darkgreen","orange"="orange","saddlebrown"="saddlebrown","sienna3"="sienna3",
                                            "paleturquoise"="paleturquoise","darkolivegreen"="darkolivegreen","palevioletred3","royalblue"="royalblue","brown"="brown","darkred"="darkred",
                                            "magenta"="magenta","navajowhite2"="navajowhite2","green"="green","red"="red","pink"="pink","greenyellow"="greenyellow",
                                            "midnightblue"="midnightblue","salmon"="salmon","yellow"="yellow","salmon4"="salmon4","skyblue3"="skyblue3","darkturquoise"="darkturquoise",
                                            "steelblue"="steelblue","darkgrey"="darkgrey","skyblue"="skyblue","maroon"="maroon","lavenderblush3"="lavenderblush3",
                                            "tan"="tan","lightgreen"="lightgreen","black"="black","lightyellow"="lightyellow","floralwhite"="floralwhite","thistle1"="thistle1",
                                            "plum2"="plum2","turquoise"="turquoise","blue"="blue","ivory"="ivory","coral1"="coral1","darkseagreen4"="darkseagreen4",
                                            "plum1"="plum1","purple"="purple","violet"="violet","bisque4"="bisque4","yellowgreen"="yellowgreen","lightsteelblue1"="lightsteelblue1",
                                            "brown4"="brown4","darkorange"="darkorange","grey60"="grey60","lightcyan1"="lightcyan1","honeydew1"="honeydew1","darkorange2"="darkorange2",
                                            "mediumpurple3"="mediumpurple3","darkslateblue"="darkslateblue","thistle2"="thistle2","darkmagenta"="darkmagenta","lightpink4"="lightpink4","orangered4"="orangered4",
                                            "white"="lightgrey","lightcyan"="lightcyan3","lightgrey"="cornsilk3")) +
                # scale_shape(name = latex2exp::TeX("$\\beta$ sign")) +
                scale_size_continuous(name="Module Size",
                                      range  = c(2,11),
                                      limits = c(min(as.numeric(wgcna_pca_dat$ModuleSize)), max(as.numeric(wgcna_pca_dat$ModuleSize))), 
                                      breaks = c(min(as.numeric(wgcna_pca_dat$ModuleSize)),
                                                 mean(as.numeric(wgcna_pca_dat$ModuleSize)),
                                                 max(as.numeric(wgcna_pca_dat$ModuleSize))),
                                      labels = c(10,50,1200)) +
                ggtitle(latex2exp::TeX("Phenotype ~ $\\ePRS_{\\Pi_{MHC-APOE-}}\\ + Covariates")) +
                labs(x = latex2exp::TeX("$\\alpha$-value threshold"), 
                     y = latex2exp::TeX("-$\\log_{10}(p)$")) +
                geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "tomato3") +
                # geom_text(aes(x = "0.1", y = -log10(0.03), 
                #               label = "p-value = 0.05"),
                #           color="tomato3") +
                geom_label_repel(data=subset(wgcna_pca_dat %>% filter(genotype_lst=="MHC-APOE-") %>% 
                                                             group_by(genotype_lst,pheno) %>% 
                                                             mutate(pfdr_assoc = ifelse(pfdr_assoc==min(pfdr_assoc),1,0)), 
                                             pfdr_assoc == 1),
                                 fill = NA,
                                aes(label=paste(Module,"--",str_sub(tolower(hubs),4))),col="black",size=3,fontface="bold")

ggsave("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/ePRS_assoc_ALL.jpg",p1,
       w=9,h=11,units = "in", dpi=700)
ggsave("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/ePRS_assoc_APOE-.jpg",p2,
       w=9,h=11, dpi=700)
ggsave("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/ePRS_assoc_MHC-.jpg",p3,
       w=9,h=11, dpi=700)
ggsave("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/ePRS_assoc_MHC-APOE-.jpg",p4,
       w=9,h=11, dpi=700)

## 7.1) Detailed Association data filter ----

# filtering for adjusted p-value < 0.05
aa <- wgcna_pca_dat %>% filter(pfdr_assoc<=0.05)
table(aa$genotype_lst)
# ALL  MHC- 
# 25   26 

### 7.1.1) If assoc results need an update, run this (takes long) ----
start = Sys.time()
# OLS_heatmap_dat_GEN_WGCNA(Study = "ROSMAP", GenoType = "No_APOE")
OLS_heatmap_dat_GEN_WGCNA(Study = "ROSMAP", GenoType = "No_MHC", dat=aa)
# OLS_heatmap_dat_GEN_WGCNA(Study = "ROSMAP", GenoType = "No_MHC_APOE")
OLS_heatmap_dat_GEN_WGCNA(Study = "ROSMAP", GenoType = "With_MHC_APOE", dat=aa)
end = Sys.time()
print(end-start)

### 7.1.2) assoc results with detailed table ----

GenoType <- "With_MHC_APOE"
path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/",GenoType)
dat1 <- read.csv(paste0(path,"/assocres_compare_hub_heatmap_dat6.csv"))
dat2 <- merge(aa %>% filter(genotype_lst=="ALL"),dat1,
              by.x=c("PVAL","phenotype"),
              by.y=c("alpha_level","pheno")) %>%
                mutate(r2_delta=r2val_full-r2val_base) %>%
                dplyr::select(geno, 
                              PVAL,
                              pheno,
                              Module,
                              hubs,
                              beta,
                              r2val_full, r2validated_full_CI_L, r2validated_base_CI_U,
                              r2_delta,
                              P_likelihood,
                              n.y) %>% 
                mutate(hubs=str_sub(hubs,4))
 
write.csv(dat2,paste0(path,"/assoc_results_detailed.csv"),row.names = F)               

GenoType <- "No_MHC"
path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/",GenoType)
dat1 <- read.csv(paste0(path,"/assocres_compare_hub_heatmap_dat.csv"))
dat2 <- merge(aa %>% filter(genotype_lst=="MHC-"),dat1,
              by.x=c("PVAL","phenotype"),
              by.y=c("alpha_level","pheno")) %>%
                mutate(r2_delta=r2val_full-r2val_base) %>%
                dplyr::select(geno, 
                              PVAL,
                              pheno,
                              Module,
                              hubs,
                              beta,
                              r2val_full, r2validated_full_CI_L, r2validated_base_CI_U,
                              r2_delta,
                              P_likelihood,
                              n.y) %>% 
                mutate(hubs=str_sub(hubs,4))

write.csv(dat2,paste0(path,"/assoc_results_detailed.csv"),row.names = F)

## 7.2) Modelling results with modules comparison plots ------

Genotype <- "With_MHC_APOE"
path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/","ROSMAP","/WGCNA/",Genotype,"/") #

dat <- read.csv(paste0(path,"assocres_compare_hub_heatmap_dat6.csv"))
# Selecting rows with q-value<0.05:
dat_gg <- dat %>% filter(modulevalues!="MEgrey") %>%
                filter(p.adjust(P_likelihood,method = "fdr")<0.05) %>%
                dplyr::select(alpha_level,
                              pheno,
                              coeff_eigenPS,
                              modulevalues,
                              r2val_base,
                              r2validated_base_CI_L,
                              r2validated_base_CI_U,
                              r2validated_full_CI_L,
                              r2validated_full_CI_U,
                              r2val_full,
                              P_likelihood,
                              n) %>%
                mutate(coeff_eigenPS=signif(coeff_eigenPS,3),
                       r2val_base=signif(r2val_base,3),
                       r2val_full=signif(r2val_full,3),
                       P_likelihood=signif(P_likelihood,3),
                       r2validated_base_CI_L=signif(r2validated_base_CI_L,3),
                       r2validated_base_CI_U=signif(r2validated_base_CI_U,3),
                       r2validated_full_CI_L=signif(r2validated_full_CI_L,3),
                       r2validated_full_CI_U=signif(r2validated_full_CI_U,3)) %>% 
                reshape2::melt(.,measure.vars=c("r2val_base","r2val_full")) %>%
                mutate(CI_U=ifelse(variable=="r2val_base",
                                   r2validated_base_CI_U,
                                   r2validated_full_CI_U),
                       CI_L=ifelse(variable=="r2val_base",
                                   r2validated_base_CI_L,
                                   r2validated_full_CI_L)) %>%
                dplyr::select(-r2validated_base_CI_L,
                              -r2validated_base_CI_U,
                              -r2validated_full_CI_L,
                              -r2validated_full_CI_U) %>%
                mutate(modulevalues = substring(modulevalues, 3,100),
                       Performance=ifelse(variable=="r2val_base",
                                          "Base Model",
                                          "Full Model"))

# if(dim(dat_gg)[1]==0){next}
indepvec.pathology <- c("cogdx","amyloid_sqrt","tangles_sqrt")
indepvec.cognition <- c("cogn_global_random_slope","cogn_globaln_lv") # Don't have this variable: "cogn_global_at_lastvisit"
pathnames <- c("Final AD","Total AB","PHF tau")
cognames <- c("Global slope","Global last visit")
varnameindex <- data.frame(name=c(pathnames,cognames),var=c(indepvec.pathology,indepvec.cognition))

dat_gg <- merge(varnameindex,dat_gg,by.x = "var",by.y="pheno")

pathnames <- c("Final AD","Total AB","PHF tau")
cognames <- c("Global slope","Global last visit") 

ggplot(dat_gg, 
       aes(x=factor(alpha_level,levels=c("1","0.1","0.05","0.01","0.005","0.001",
                                         "0.0005","0.0001","5e-05","1e-05","5e-06",
                                         "1e-06","5e-07","1e-07","5e-08")), 
           y=100*value, 
           fill=Performance)) + 
                geom_bar(stat="identity", position=position_dodge()) +
                geom_errorbar(aes(ymin=100*CI_L, ymax=100*CI_U), width=.2,
                              position=position_dodge(.9)) +
                scale_fill_brewer(palette="Paired") + 
                theme_bw() + 
                facet_wrap(~name+modulevalues,nrow = 4, scales = "free") +
                theme(axis.text.x=element_text(angle = -45, hjust = 0))+
                ylab(latex2exp::TeX("Model performance ($\\R^2,AUC$)")) + #  or $\\R^2$ $\\AUC$
                xlab(latex2exp::TeX("$\\alpha$-value threshold")) +
                ggtitle(latex2exp::TeX("Phenotype ~ $(\\ePRS_{Pi_{ALL}})+\\PRS_{LOAD}+Covariates$")) + 
                theme(plot.title = element_text(hjust = 0.5))
# geom_text(aes(label=paste(100*value,"%")), vjust=1.6, 
#           color="white",
#           position = position_dodge(0.9), 
#           size=1.6)
# data = dat_gg %>% 
# mutate(P_likelihood=ifelse(variable=="r2val_full",
# "Likelihood Ratio Test P-value:",
# P_likelihood)))

ggsave(filename = paste0(path,"/full_association_analysis_AUC_",Genotype,".jpg"), 
       width = 8, height = 10, dpi=500) #width=12

Genotype <- "No_MHC"
path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/","ROSMAP","/WGCNA/",Genotype,"/") #

dat <- read.csv(paste0(path,"assocres_compare_hub_heatmap_dat6.csv"))
# Selecting rows with q-value<0.05:
dat_gg <- dat %>% filter(modulevalues!="MEgrey") %>%
                filter(p.adjust(P_likelihood,method = "fdr")<0.05) %>%
                dplyr::select(alpha_level,
                              pheno,
                              coeff_eigenPS,
                              modulevalues,
                              r2val_base,
                              r2validated_base_CI_L,
                              r2validated_base_CI_U,
                              r2validated_full_CI_L,
                              r2validated_full_CI_U,
                              r2val_full,
                              P_likelihood,
                              n) %>%
                mutate(coeff_eigenPS=signif(coeff_eigenPS,3),
                       r2val_base=signif(r2val_base,3),
                       r2val_full=signif(r2val_full,3),
                       P_likelihood=signif(P_likelihood,3),
                       r2validated_base_CI_L=signif(r2validated_base_CI_L,3),
                       r2validated_base_CI_U=signif(r2validated_base_CI_U,3),
                       r2validated_full_CI_L=signif(r2validated_full_CI_L,3),
                       r2validated_full_CI_U=signif(r2validated_full_CI_U,3)) %>% 
                reshape2::melt(.,measure.vars=c("r2val_base","r2val_full")) %>%
                mutate(CI_U=ifelse(variable=="r2val_base",
                                   r2validated_base_CI_U,
                                   r2validated_full_CI_U),
                       CI_L=ifelse(variable=="r2val_base",
                                   r2validated_base_CI_L,
                                   r2validated_full_CI_L)) %>%
                dplyr::select(-r2validated_base_CI_L,
                              -r2validated_base_CI_U,
                              -r2validated_full_CI_L,
                              -r2validated_full_CI_U) %>%
                mutate(modulevalues = substring(modulevalues, 3,100),
                       Performance=ifelse(variable=="r2val_base",
                                          "Base Model",
                                          "Full Model"))

# if(dim(dat_gg)[1]==0){next}
indepvec.pathology <- c("cogdx","amyloid_sqrt","tangles_sqrt")
indepvec.cognition <- c("cogn_global_random_slope","cogn_globaln_lv") # Don't have this variable: "cogn_global_at_lastvisit"
pathnames <- c("Final AD","Total AB","PHF tau")
cognames <- c("Global slope","Global last visit")
varnameindex <- data.frame(name=c(pathnames,cognames),var=c(indepvec.pathology,indepvec.cognition))

dat_gg <- merge(varnameindex,dat_gg,by.x = "var",by.y="pheno")

pathnames <- c("Final AD","Total AB","PHF tau")
cognames <- c("Global slope","Global last visit") 

ggplot(dat_gg, 
       aes(x=factor(alpha_level,levels=c("1","0.1","0.05","0.01","0.005","0.001",
                                         "0.0005","0.0001","5e-05","1e-05","5e-06",
                                         "1e-06","5e-07","1e-07","5e-08")), 
           y=100*value, 
           fill=Performance)) + 
                geom_bar(stat="identity", position=position_dodge()) +
                geom_errorbar(aes(ymin=100*CI_L, ymax=100*CI_U), width=.2,
                              position=position_dodge(.9)) +
                scale_fill_brewer(palette="Paired") + 
                theme_bw() + 
                facet_wrap(~name+modulevalues,nrow = 4, scales = "free") +
                theme(axis.text.x=element_text(angle = -45, hjust = 0))+
                ylab(latex2exp::TeX("Model performance ($\\R^2,AUC$)")) + #  or $\\R^2$ $\\AUC$
                xlab(latex2exp::TeX("$\\alpha$-value threshold")) +
                ggtitle(latex2exp::TeX("Phenotype ~ $(\\ePRS_{Pi_{MHC-}})+\\PRS_{LOAD}+Covariates$")) + 
                theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = paste0(path,"/full_association_analysis_AUC_",Genotype,".jpg"), 
       width = 10, height = 8, dpi=500) #width=12
