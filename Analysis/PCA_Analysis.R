# Libraries ----
library(data.table)
library(openxlsx)
library(rms)
library(ggplot2)
library(ggrepel) # for overlapping labels in ggplot
library(corrplot)
library(factoextra) # for PCA interpretation
library(colorspace) # get nice colors
library(cluster)    # clustering algorithms
library(Hmisc)
library(dplyr) # data manipulation
library(tidyverse)  # data manipulation
library(summarytools)
library(ggrepel) # for overlapping labels in ggplot
library(rms)
library(broom) # tidy function


# 1) pca_results_GEN function ----
  # This function generates the info we need for PCA results:
pca_results_GEN <- function(Study=c("ROSMAP","ADNI"),
                            GenoType=c("No_APOE","No_MHC","No_MHC_APOE","With_MHC_APOE"),
                            trait_type = FALSE){
                
                path <- paste0("../Datasets/CLUMP_500_0.2/",Study,"/Resid_PRS/",GenoType)
                path_to_save <- paste0("../Datasets/CLUMP_500_0.2/",Study,"/PCA/",GenoType)
                
                # read in residuals of PRSs files
                lmatrix <- readRDS(paste0(path,"/Residual_results_all_p-vals.rds"))
                
                pca_results <- matrix(ncol = 7)
                for (prs_matrix in names(lmatrix)){
                                print(prs_matrix)
                                # removing IID columns
                                matrix <- lmatrix[[prs_matrix]]$residuals[,-dim(lmatrix[[prs_matrix]]$residuals)[2]]
                                
                                respca <- prcomp(matrix)
                                pc1_cont <- sum(summary(respca)$importance[2,1])*100
                                pc2_cont <- sum(summary(respca)$importance[2,2])*100
                                top_5_cont <- sum(summary(respca)$importance[2,1:5])*100
                                top_10_cont <- sum(summary(respca)$importance[2,1:10])*100
                                for(r in 1:dim(summary(respca)$importance)[2]){
                                                if(sum(summary(respca)$importance[2,1:r])>0.80){
                                                                top_r_80_explained <- r
                                                                break
                                                }
                                }
                                row <- c(prs_matrix,pc1_cont,pc2_cont,top_5_cont,top_10_cont,r,dim(summary(respca)$importance)[2])
                                pca_results <- rbind(pca_results,row)
                                
                                for(r in 5:dim(summary(respca)$importance)[2]){
                                                if(sum(summary(respca)$importance[2,1:r])>0.80){
                                                                top_r_80_explained <- r
                                                                break
                                                }
                                }
                                # ggsave(fviz_eig(respca,choice="variance",ncp = top_r_80_explained,main = paste0("Scree plot of PRSs calculated at P-value threshold of ",
                                #                                                                                 prs_matrix," for top ",top_r_80_explained, " PCs")),
                                #            file=paste0(path_to_save,"/scree_plot_",prs_matrix,".jpg"),
                                #        width = 250,
                                #        height = 100,
                                #        units = "mm")
                                
                }
                
                colnames(pca_results) <- c("P-value_thresh","PC1_var_explained(%)","PC2_var_explained(%)","PC1-PC5_var_explained(%)",
                                           "PC1-PC10_var_explained(%)","top_r_comps_explain_80%","total_components")
                pca_results <- pca_results[2:dim(pca_results)[1],]
                pca_results <- as.data.frame(pca_results)
                write.csv(pca_results,paste0(path_to_save,'/pca_summary.csv'), row.names = FALSE)
}

pca_results_GEN(Study = "ROSMAP",GenoType = "With_MHC_APOE")
pca_results_GEN(Study = "ROSMAP",GenoType = "No_MHC")
pca_results_GEN(Study = "ROSMAP",GenoType = "No_APOE")
pca_results_GEN(Study = "ROSMAP",GenoType = "No_MHC_APOE")


# 2) pca_association function ----
  # This function will find the p-values of each PC for each selected phenotype in ROSMAP, and generate plots for them
pca_association <- function(Study=c("ROSMAP","ADNI"),
                            GenoType=c("No_APOE","No_MHC","No_MHC_APOE","With_MHC_APOE"),
                            PCnum=25){
                path <- paste0("../Datasets/CLUMP_500_0.2/",Study,"/Resid_PRS/",GenoType)
                path_to_save <- paste0("../Datasets/CLUMP_500_0.2/",Study,"/PCA/",GenoType)
                
                # Reading ROS/MAP phenotype dataset
                ROSmaster <- readRDS("../Datasets/ROSMAP_Phenotype/ROSmaster.rds")
                # Reading Filtered PNUKBB manifest
                meta_pheno <- read.csv("./Pan_UKBB/Data/ukbb_manifest_filtered_phenos.csv")
                # Reading PCA of the ROS/MAP Genotype dataset
                geno_pcs <- read.table("../Datasets/PCA_Genotype/geno_qc.eigenvec_new.txt",header=F)
                names(geno_pcs) <- c("FID","IID",paste0("genoPC",seq(1:10)))
                
                # read in residuals of PRSs files
                lmatrix <- readRDS(paste0(path,"/Residual_results_all_p-vals.rds"))
                
                results.G <- list()
                for (prs_matrix in names(lmatrix)){
                                print(prs_matrix)
                                matrix <- lmatrix[[prs_matrix]]$residuals[,-dim(lmatrix[[prs_matrix]]$residuals)[2]]
                                
                                respca <- prcomp(matrix)
                                res.pcs <- as.data.frame(respca$x)
                                res.pcs$IID <- lmatrix[[prs_matrix]]$residuals$IID
                                pheno_pcs <- merge(res.pcs,ROSmaster,by="IID")
                                pheno_pcs2 <- merge(pheno_pcs,geno_pcs,by="IID")
                                ## identify relationships of PCs with phenotypes of interest
                                pheno.list <- c("amyloid_sqrt","tangles_sqrt","cogn_global_random_slope","age_death","parksc_lv_sqrt",
                                                "cogdx","gpath","educ","thyroid_ever", "stroke_ever","diabetes_sr_rx_ever",
                                                "cancer_ever","hypertension_ever")
                                
                                if (dim(summary(respca)$importance)[2] <= PCnum){
                                                PClist <- paste0("PC",seq(1:dim(summary(respca)$importance)[2]))
                                } else {
                                                PClist <- paste0("PC",seq(1:PCnum))
                                }
                                
                                index <- 1
                                pvalues <- NULL
                                bvalues <- NULL
                                nvalues <- NULL
                                phenovalues <- NULL
                                pcvalues <- NULL
                                for (pheno in pheno.list) {
                                                for (PC in PClist) {
                                                                form <- formula(paste(pheno,"~",PC))
                                                                mod <- lm(data=pheno_pcs2,form)
                                                                #mod1 <- glm(data=pheno_pcs2,form,family = binomial)
                                                                pvalues[index] <- anova(mod)[1,5]
                                                                bvalues[index] <- coef(mod)[2]
                                                                nvalues[index] <- dim(mod$model)[1]
                                                                pcvalues[index] <- PC
                                                                phenovalues[index] <- pheno
                                                                index <- index + 1
                                                }
                                }
                                assocres <- data.frame(pheno=phenovalues,
                                                       pc=pcvalues,
                                                       b=bvalues,
                                                       p=pvalues,
                                                       n=nvalues,
                                                       fdr=p.adjust(pvalues,method="fdr"))
                                
                                assocres$colour <- ifelse(assocres$b < 0, "Negative effect","Positive effect")
                                g <- ggplot(data=assocres, aes(x=pheno,y=-log10(p),group=pc))+
                                                geom_bar(stat="identity",position="dodge",aes(fill = colour))+
                                                geom_hline(aes(yintercept = -log10(0.05/nrow(assocres)),color="Bonferroni corrected P-value"),lty=2)+
                                                geom_hline(aes(yintercept = -log10(0.05),color="P-value < 0.05"),lty=2)+
                                                scale_linetype_manual(name = "limit", values = c(2, 2),
                                                                      guide = guide_legend(override.aes = list(color = c("green", "orange"))))+
                                                geom_text_repel(data=subset(assocres,p<0.05),aes(label=pc))+
                                                theme_minimal()+
                                                theme(axis.text.x=element_text(angle = -45, hjust = 0))+
                                                ggtitle(paste0('PRSs at P-value threshhold of',"\n",
                                                               prs_matrix))
                                # ggsave(paste0(path_to_save,"/Associations_top_",PCnum,"_PCs_at_", prs_matrix,".jpg"),
                                #        width = 54, height = 27, units = "cm")
                                results.G[[prs_matrix]] <- list(pc=respca,
                                                        resultsols=assocres,
                                                        plot=g)
                }
                saveRDS(results.G,file=paste0(path_to_save,"/PCA_results_all_p-vals",".rds"))
}
pca_association(Study = "ROSMAP",GenoType = "With_MHC_APOE")
pca_association(Study = "ROSMAP",GenoType = "No_MHC")
pca_association(Study = "ROSMAP",GenoType = "No_APOE")
pca_association(Study = "ROSMAP",GenoType = "No_MHC_APOE")

# 3) pca_heatmap_dat_GEN function ----
  # This function generates the data required to create the PCA heatmap explainer
pca_heatmap_dat_GEN <- function(Study=c("ROSMAP","ADNI"),
                                GenoType=c("No_APOE","No_MHC","No_MHC_APOE","With_MHC_APOE"),
                                PCNum = 25,
                                topcont = 20){
                
                path_pca <- paste0("../Datasets/CLUMP_500_0.2/",Study,"/PCA/",GenoType)
                # path_to_save_plot <- paste0("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Plots/PCA/",Study,"/",GenoType,"/Heatmaps")
                
                # Reading ROS/MAP phenotype dataset
                ROSmaster <- readRDS("../Datasets/ROSMAP_Phenotype/ROSmaster.rds")
                # Reading Filtered PNUKBB manifest
                meta_pheno <- read.csv("./Pan_UKBB/Data/ukbb_manifest_filtered_phenos.csv")
                # Reading PCA of the ROS/MAP Genotype dataset
                geno_pcs <- read.table("../Datasets/PCA_Genotype/geno_qc.eigenvec_new.txt",header=F)
                names(geno_pcs) <- c("FID","IID",paste0("genoPC",seq(1:10)))
                
                results.G <- readRDS(paste0(path_pca,"/PCA_results_all_p-vals",".rds"))
                
                pc_vector = NULL #----> Z
                pc_name = NULL #----> X
                alpha_val = NULL #-----> Y
                PC_contributors <- NULL
                top_pc_explain_80 <- NULL
                
                for (prs_matrix in names(results.G)){
                                respca <- results.G[[prs_matrix]]$pc
                                res.var <- get_pca_var(respca)
                                if(PCNum>=dim(respca$x)[2]){Num=dim(respca$x)[2]}else{Num<-PCNum}
                                top_k_80_explained=NULL
                                for(pc in 1:dim(respca$x)[2]){
                                                top_k_80_explained[pc] <- sum(summary(respca)$importance[2,pc])*100
                                                if(sum(top_k_80_explained)>=80){break}
                                }
                                for(pc in 1:Num){
                                                pc_vector <- append(pc_vector,sum(summary(respca)$importance[2,pc])*100)
                                                pc_name <- append(pc_name,paste0("PC",pc))
                                                alpha_val <- append(alpha_val,prs_matrix)
                                                top_pc_explain_80 <- append(top_pc_explain_80,length(top_k_80_explained))
                                                PC_contributors <- append(PC_contributors,paste0(paste0(as.character(names(res.var$contrib[,pc][order(res.var$contrib[,pc],decreasing = T)][1:topcont])),
                                                                                                        " -- ICD10 Code: ",
                                                                                                        meta_pheno[match(as.character(names(res.var$contrib[,pc][order(res.var$contrib[,pc],decreasing = T)][1:topcont])), meta_pheno$new_pheno_annot),]$ICD10,
                                                                                                        " -- Contribution: ",
                                                                                                        as.character(signif(res.var$contrib[,pc][order(res.var$contrib[,pc],decreasing = T)][1:topcont],digits = 3)),"%"), collapse="\n"))
                                }
                }
                
                pc_info_dat <- data.frame(pc_name,alpha_val,pc_vector,PC_contributors,top_pc_explain_80)
                pc_info_dat$alpha_val <- gsub("p_val_","",pc_info_dat$alpha_val)
                pc_info_dat$alpha_val <- gsub("_",".",pc_info_dat$alpha_val)
                
                write.csv(pc_info_dat,paste0(path_pca,"/pc_heatmap_dat.csv"),row.names = FALSE)
                
}
pca_heatmap_dat_GEN(PCNum = 2500,Study = "ROSMAP", GenoType = "With_MHC_APOE")
pca_heatmap_dat_GEN(PCNum = 2500,Study = "ROSMAP", GenoType = "No_MHC")
pca_heatmap_dat_GEN(PCNum = 2500,Study = "ROSMAP", GenoType = "No_APOE")
pca_heatmap_dat_GEN(PCNum = 2500,Study = "ROSMAP", GenoType = "No_MHC_APOE")


# 4) pca_assoc_dat_GEN function----
pca_assoc_dat_GEN <- function(Study=c("ROSMAP","ADNI"),
                              GenoType=c("No_APOE","No_MHC","No_MHC_APOE","With_MHC_APOE"),
                              PCNum = 25,
                              topcont = 25){
                path_pca <- paste0("../Datasets/CLUMP_500_0.2/",Study,"/PCA/",GenoType)
                path <- paste0("../Datasets/CLUMP_500_0.2/",Study,"/Resid_PRS/",GenoType)
                
                # Reading ROS/MAP phenotype dataset
                ROSmaster <- readRDS("../Datasets/ROSMAP_Phenotype/ROSmaster.rds")
                # Reading Filtered PNUKBB manifest
                meta_pheno <- read.csv("./Pan_UKBB/Data/ukbb_manifest_filtered_phenos.csv")
                # Reading PCA of the ROS/MAP Genotype dataset
                geno_pcs <- read.table("../Datasets/PCA_Genotype/geno_qc.eigenvec_new.txt",header=F)
                names(geno_pcs) <- c("FID","IID",paste0("genoPC",seq(1:10)))
                
                
                # read in residuals of PRSs files
                lmatrix <- readRDS(paste0(path,"/Residual_results_all_p-vals.rds"))
                
                results.G <- readRDS(paste0(path_pca,"/PCA_results_all_p-vals",".rds"))
                
                index <- 1
                pvalues <- NULL
                bvalues <- NULL
                nvalues <- NULL
                phenovalues <- NULL
                egenevalues <- NULL
                pcvalues <- NULL
                PC_contributors <- NULL
                alpha_val <- NULL
                
                #PC_contributors_value <- NULL
                for (prs_matrix in names(results.G)){
                                respca <- results.G[[prs_matrix]]$pc
                                res.var <- get_pca_var(respca)
                                res.pcs <- as.data.frame(respca$x)
                                res.pcs$IID <- lmatrix[[prs_matrix]]$residuals$IID
                                pheno_pcs <- merge(res.pcs,ROSmaster,by="IID")
                                #pheno_pcs2 <- merge(pheno_pcs,geno_pcs,by="IID")
                                if (dim(summary(respca)$importance)[2] <= PCNum){
                                                PClist <- paste0("PC",seq(1:dim(summary(respca)$importance)[2]))
                                } else {
                                                PClist <- paste0("PC",seq(1:PCNum))
                                }
                                
                                # Changing Cogdx variable:
                                # variable = cogdx
                                # Q: Overall cognitive diagnostic category?
                                #                 Coding Detail
                                # 1 NCI, No cognitive impairment (No impaired domains)
                                # 2 MCI, Mild cognitive impairment (One impaired domain) and NO other cause of CI
                                # 3 MCI, Mild cognitive impairment (One impaired domain) AND another cause of CI
                                # 4 AD, Alzheimer's disease and NO other cause of CI (NINCDS PROB AD)
                                # *5 AD, Alzheimer's disease AND another cause of CI (NINCDS POSS AD)
                                # *6 Other dementia. Other primary cause of dementia
                                
                                pheno_pcs$cogdx[which((pheno_pcs$cogdx==2)|(pheno_pcs$cogdx==3)|(pheno_pcs$cogdx==5)|(pheno_pcs$cogdx==6))] <- NA
                                pheno_pcs$cogdx[which((pheno_pcs$cogdx==4))] <- 2
                                
                                #### define outcomes and covariates
                                covars.pathology <- c("msex","pmi","age_death")
                                covars.cognition <- c("educ","msex","age_death")
                                covars.death <- c("msex","age_bl")
                                
                                # indepvec.pathology <- c("plaq_n_sqrt","plaq_d_sqrt","amyloid_sqrt","tangles_sqrt","nft_sqrt",
                                #                         "ci_num2_gct","ci_num2_mct","arteriol_scler","caa_4gp","cvda_4gp2",
                                #                         "dlbdx","hspath_any","tdp_stage4","parkdx","pathoAD","vm3123",
                                #                         "pput3123","it3123","mf3123")
                                
                                indepvec.pathology <- c("plaq_n_sqrt","plaq_d_sqrt","cogdx","amyloid_sqrt","tangles_sqrt",
                                                        "pathoAD","nft_sqrt","dlbdx")
                                
                                indepvec.cognition <- c("cogn_global_random_slope") # Don't have this variable: "cogn_global_at_lastvisit"
                                
                                indepvec.death <- c("age_death")
                                
                                pathnames <- c("Neuritic plaques","Diffuse plaques","Final AD","Total AB","PHF tau",
                                               "Patho AD","NFT","Lewy body stage")
                                
                                # pathnames <- c("Neuritic plaques","Diffuse plaques","Total AB","PHF tau","NFT","Gross cerebral infarcts",
                                #                "Micro cerebral infarcts","Arteriolosclerosis","Cerebral AA","Cerebral atherosclerosis",
                                #                "Lewy body stage","Hippocampal sclerosis","TDP-43","PD Dx","Patho AD","PAM VM Caudate",
                                #                "PAM post. putamen","PAM IT","PAM MF")
                                
                                cognames <- c("Global slope") #,"Global last visit"
                                
                                deathnames <- c("Age of death")
                                
                                varnameindex <- data.frame(name=c(pathnames,cognames,deathnames),var=c(indepvec.pathology,indepvec.cognition,indepvec.death))
                                
                                nocovars <- FALSE
                                
                                print("Running Associations")
                                
                                for(vartype in c("path","cog","deat")) {
                                                if (vartype=="path") {
                                                                indepvec <- indepvec.pathology
                                                                covars <- covars.pathology
                                                } else if (vartype=="cog") {
                                                                indepvec <- indepvec.cognition
                                                                covars <- covars.cognition
                                                } else {
                                                                indepvec <- indepvec.death
                                                                covars <- covars.death
                                                }
                                                
                                                if (length(indepvec)>1){
                                                                tablist <- apply(pheno_pcs[,indepvec],2,table) %>%
                                                                                lapply(length) %>%
                                                                                as.numeric()
                                                } else {
                                                                tablist <- table(pheno_pcs[,indepvec])  %>%
                                                                                lapply(length) %>%
                                                                                as.numeric() %>% sum()      
                                                }
                                                
                                                cont.index <- tablist > 2 #TRUE if numeric indepvec
                                                phenindex <- 1
                                                for (pheno in indepvec) {
                                                                if (cont.index[phenindex]==T) {
                                                                                for (PC in PClist) {
                                                                                                if (nocovars==TRUE) {
                                                                                                                form <- formula(paste(pheno,"~",PC))
                                                                                                } else {
                                                                                                                form <- formula(paste(pheno,"~",PC,"+",paste0(covars,collapse = " + ")))
                                                                                                }
                                                                                                pcnum <- as.integer(sub('PC', '', PC))
                                                                                                PC_contributors <- append(PC_contributors,paste0(paste0(as.character(names(res.var$contrib[,pcnum][order(res.var$contrib[,pcnum],decreasing = T)][1:topcont])),
                                                                                                                                                        " -- ICD10 Code: ",
                                                                                                                                                        meta_pheno[match(as.character(names(res.var$contrib[,pcnum][order(res.var$contrib[,pcnum],decreasing = T)][1:topcont])), meta_pheno$new_pheno_annot),]$ICD10,
                                                                                                                                                        " -- Contribution: ",
                                                                                                                                                        as.character(signif(res.var$contrib[,pcnum][order(res.var$contrib[,pcnum],decreasing = T)][1:topcont],digits = 3)),"%"), collapse="\n"))
                                                                                                
                                                                                                mod <- lm(data=pheno_pcs, form)
                                                                                                pvalues[index] <- tidy(mod)$p.value[2]
                                                                                                bvalues[index] <- coef(mod)[2]
                                                                                                nvalues[index] <- dim(mod$model)[1]
                                                                                                pcvalues[index] <- PC
                                                                                                phenovalues[index] <- pheno
                                                                                                alpha_val[index] <- prs_matrix
                                                                                                index <- index+1
                                                                                }
                                                                } else {
                                                                                for (PC in PClist) {
                                                                                                if (nocovars==TRUE) {
                                                                                                                form <- formula(paste(pheno,"~",PC))
                                                                                                } else {
                                                                                                                form <- formula(paste(pheno,"~",PC,"+",paste0(covars,collapse = " + ")))
                                                                                                }
                                                                                                pcnum <- as.integer(sub('PC', '', PC))
                                                                                                PC_contributors <- append(PC_contributors,paste0(paste0(as.character(names(res.var$contrib[,pcnum][order(res.var$contrib[,pcnum],decreasing = T)][1:topcont])),
                                                                                                                                                        " -- ICD10 Code: ",
                                                                                                                                                        meta_pheno[match(as.character(names(res.var$contrib[,pcnum][order(res.var$contrib[,pcnum],decreasing = T)][1:topcont])), meta_pheno$new_pheno_annot),]$ICD10,
                                                                                                                                                        " -- Contribution: ",
                                                                                                                                                        as.character(signif(res.var$contrib[,pcnum][order(res.var$contrib[,pcnum],decreasing = T)][1:topcont],digits = 3)),"%"), collapse="\n"))
                                                                                                
                                                                                                mod <- lrm(data=pheno_pcs, form)
                                                                                                pvalues[index] <- anova(mod)[1,3]
                                                                                                bvalues[index] <- coef(mod)[2]
                                                                                                nvalues[index] <- mod$stats[1]
                                                                                                pcvalues[index] <- PC
                                                                                                phenovalues[index] <- pheno
                                                                                                alpha_val[index] <- prs_matrix
                                                                                                index <- index+1
                                                                                }
                                                                }
                                                                phenindex <- phenindex + 1
                                                }
                                }
                                
                }
                assocres <- data.frame(alpha_level=alpha_val,
                                       pheno=phenovalues,
                                       pc=pcvalues,
                                       b=bvalues,
                                       p=pvalues,
                                       n=nvalues,
                                       signedp=-log10(pvalues)*sign(bvalues),
                                       PC_contributors=PC_contributors)
                write.csv(assocres,paste0(path_pca,"/assocres_pc_heatmap_dat.csv"),row.names = FALSE)
}
pca_assoc_dat_GEN(Study = "ROSMAP", GenoType = "With_MHC_APOE")
pca_assoc_dat_GEN(Study = "ROSMAP", GenoType = "No_MHC")
pca_assoc_dat_GEN(Study = "ROSMAP", GenoType = "No_APOE")
pca_assoc_dat_GEN(Study = "ROSMAP", GenoType = "No_MHC_APOE")


# (Don't run) pca_contribution function ----
# This function will plot contributions of top n PCs for each selected PRS matrix
pca_contribution <- function(Study=c("ROSMAP","ADNI"),
                             GenoType=c("No_APOE","No_MHC","No_MHC_APOE","With_MHC_APOE"),
                             label_subset = levels(as.factor(meta_pheno$category_manual)),
                             trait_type = FALSE,
                             PCnum=25){
                path_pca <- paste0("../Datasets/CLUMP_500_0.2/",Study,"/PCA/",GenoType)
                path_to_save <- paste0("../Datasets/CLUMP_500_0.2/",Study,"/PCA/",GenoType)
                
                results.G <- readRDS(paste0(path_pca,"/PCA_results_all_p-vals",".rds"))
                for (p in names(results.G)){
                                print(p)
                                respca <- results.G[[p]]$pc
                                res.var <- get_pca_var(respca)
                                
                                if (dim(summary(respca)$importance)[2] <= PCnum){
                                                N <- dim(summary(respca)$importance)[2]
                                } else {
                                                N <- PCnum
                                }
                                
                                for (pcnum in 1:N){
                                                var <-sum(summary(respca)$importance[2,pcnum])*100
                                                ### get contributing PRS for given PC
                                                pc_top20 <- as.data.frame(cbind(res.var$cos2[,pcnum][order(res.var$cos2[,pcnum],decreasing = T)][1:10],res.var$contrib[,pcnum][order(res.var$contrib[,pcnum],decreasing = T)][1:10]))
                                                pc_top20$phenotypes <- rownames(pc_top20)
                                                names(pc_top20) <- c("Cos2","Contribution","phenotypes")
                                                # represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
                                                
                                                # contains the contributions (in percentage) of the variables to the principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
                                                g <- ggplot(pc_top20, aes(x=reorder(phenotypes, Contribution), y=Contribution)) +
                                                                geom_bar(stat="identity",position="dodge")+
                                                                theme_minimal()+
                                                                theme(axis.text.x=element_text(angle = -90, hjust = 0))+
                                                                ggtitle(paste0('Contributions of top 10 variables','\n', 'to PC',pcnum,
                                                                               ' (',var,'% variation explained',')'))+
                                                                ylab("Contributions of each phenotype (%)")+
                                                                xlab("")+
                                                                coord_flip()
                                                #print(g)
                                                ggsave(paste0(path_to_save,"/top_",PCnum,"_contribution_at_", p ,"_threshold_PC_",pcnum,".jpg"))
                                }
                }
}
pca_contribution(Study = "ROSMAP",GenoType = "With_MHC_APOE")
pca_contribution(Study = "ROSMAP",GenoType = "No_MHC")
pca_contribution(Study = "ROSMAP",GenoType = "No_APOE")
pca_contribution(Study = "ROSMAP",GenoType = "No_MHC_APOE")

