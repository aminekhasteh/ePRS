# Libraries ----
library(data.table)
library(rms)
library(limma)
library(ggplot2)
library(diptest)
library(lsr)
library(dplyr)
library(harrietr) # dist to long format
library(HiClimR) # fast cor
library(WGCNA)

# Data cleaning function ----
  # This function applies initial QC on the PRS matrices for PC and WGCN analysis
primary_resid_data_GEN <- function(Study=c("ROSMAP","ADNI"),
                                   GenoType=c("No_APOE","No_MHC","No_MHC_APOE","With_MHC_APOE"),
                                   use_snp_count=TRUE,
                                   multimodal=TRUE,
                                   snp_count_n=5,
                                   cor_thresh =0.8) {
                # Reading ROS/MAP phenotype dataset
                ROSmaster <- readRDS("../Datasets/ROSMAP_Phenotype/ROSmaster.rds")
                # Reading batch effect data:
                batch_effect <- fread("../Datasets/batch_index.txt",col.names = c("IID","batch"))
                ROSmaster <- merge(ROSmaster,batch_effect,by.x="IID",by.y="IID")
                # Reading Filtered PNUKBB manifest
                meta_pheno <- read.csv("../ePRS/Pan_UKBB/Data/ukbb_manifest_filtered_phenos.csv")
                # Reading PCA of the ROS/MAP Genotype dataset
                geno_pcs <- read.table("../Datasets/PCA_Genotype/geno_qc.eigenvec_new.txt",header=F)
                names(geno_pcs) <- c("FID","IID",paste0("genoPC",seq(1:10)))

                # ----------------------------------------------------------- # Reading PRS matrices # ----------------------------------------------------------- #                                                                                       #|

                path <- paste0("../Datasets/CLUMP_500_0.2/",Study,"/PRS/",GenoType)
                path_to_save <- paste0("../Datasets/CLUMP_500_0.2/",Study,"/Resid_PRS/",GenoType)
                
                if(isTRUE(use_snp_count)){
                                snp_count <- as.data.frame(fread(paste0(path,'/SNP_count/snp_count.txt'),header=T))
                                snp_count$category <- c("1","0.1","0.05","0.01",
                                                        "0.005","0.001","0.0005","0.0001",
                                                        "5e-05","1e-05","5e-06","1e-06",
                                                        "5e-07","1e-07","5e-08")
                }
                #------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #
                # read in PRS files
                prs_filenames <- list.files(path,full.names = T)
                prs_list <- list()
                for (prsfile in prs_filenames[which(grepl(".txt", prs_filenames))]) {
                                prs_list[[tail(strsplit(prsfile,split="/")[[1]],n=1)]] <- as.data.frame(fread(prsfile,header=T,sep=","))
                }

                ## Loop to generate the residuals and the PCA datasets:
                dsets_torun <- names(prs_list)
                allres <- list()
                details <- matrix(ncol = 6)
                matrix0 <- prs_list$`p_val_5e-08.txt` # setting matrix0 to lowest p-value threshold
                for (thresh in dsets_torun[which(grepl("p_val",dsets_torun))][order(as.numeric(gsub("p_val_","",
                                                                                                    gsub(".txt","",
                                                                                                         dsets_torun[which(grepl("p_val",dsets_torun))]))),decreasing = F)]) {
                                pval <- gsub(".txt","",thresh)
                                pval <- gsub("p_val_","",pval)
                                
                                # PRSice removes any phenotype if at a specific p-value threshold the PRS does not add any additional information compared to the one at lower p-value threshold
                                # We want to still add those PRSs regardless of if they don't add additional information to the higher p-value treshold
                                matrix1 <- prs_list[[thresh]]
                                names <- names(matrix0)[which(names(matrix0) %nin% names(matrix1))]
                                matrix1 <- as.data.frame(cbind(matrix0[names],matrix1))
                                matrix0 <- prs_list[[thresh]]
                                prs_list[[thresh]] <- matrix1
                                # need to replace prs_list[[thresh]] by matrix1 
                                
                                print(paste("Starting run, p-value threshold=",pval))
                                print(paste("number of subjects =",nrow(prs_list[[thresh]])))
                                print(paste("number of PRS =",ncol(prs_list[[thresh]])-2))
                                #names(prs_list[[thresh]]) <- gsub(",","",names(prs_list[[thresh]]),fixed=T)
                                
                                formod <- merge(prs_list[[thresh]],geno_pcs,by=c("IID","FID"))
                                formod <- merge(formod,ROSmaster,by=c("IID")) # 7 user removed here...
                                formod <- formod[,which(names(formod)!="V1")]
                                formod <- formod[which(formod$IID != "ROS10424475"),] # removing this individual since they have a different race.
                                formod$batch <- as.factor(formod$batch)
                                
                                # this could be removed from loop, but in case any subjects are removed/added to either geno_pcs or PRS files, the design matrix will need to be recalculated, so left in to prevent error
                                designvars <- c(names(geno_pcs)[-c(1:2)])  
                                designtext <- paste("model.matrix( ~ ",paste("formod$",designvars,sep="",collapse=" + ")," + formod$batch)",sep="")
                                design <- eval(parse(text=designtext))
                                
                                # OLS regression using machinery build for gene expression in limma
                                print(paste("Fitting model with genetic PCs"))
                                fit <- lmFit(t(formod[,3:ncol(prs_list[[thresh]])]),design)
                                
                                print(paste("Extracting model residuals"))
                                res <- t(residuals.MArrayLM(fit,y = t(formod[,3:ncol(prs_list[[thresh]])])))
                                
                                res <- scale(res,scale=TRUE)

                                #Finding columns with NAs
                                na_cols <- colnames(res)[ apply(res, 2, anyNA) ]
                                print(paste0('number of columns with NAs ',length(na_cols)))
                                res <- res[,which(!colnames(res) %in% na_cols)]
                                
                                #SNP count flag:
                                if (isTRUE(use_snp_count)){
                                                pheno_failed <- names(snp_count[,which(snp_count[which(snp_count$category == pval),] < snp_count_n)])[which(!names(snp_count[,which(snp_count[which(snp_count$category == pval),] < snp_count_n)]) %in% "category")]
                                                print(paste("Number of phenotypes with less than",snp_count_n,"SNPs used in each phenotype is",length(pheno_failed))) # removing 1 : category column
                                                
                                                if (length(pheno_failed)>0) {
                                                                low_snp_pheno <- sum(colnames(res) %in% pheno_failed)
                                                                res <- res[,colnames(res)[which(!colnames(res) %in% pheno_failed)]]
                                                } else {
                                                                low_snp_pheno <- 0
                                                                lowsnpcount.names <- NA
                                                                res <- res
                                                }
                                                
                                } else{low_snp_pheno<-0}
                                
                                # Change the phenotype labels
                                i1 <- match(colnames(res), meta_pheno$phenocode_annotate_lst)
                                i2 <- !is.na(i1) # to take care of non matches which are NA
                                colnames(res)[i2] <- meta_pheno$new_pheno_annot[i1[i2]]
                                
                                # Removing phenotypes we selected manually
                                phenos_to_remove <- which(colnames(res) %in% meta_pheno$new_pheno_annot[which(meta_pheno$to_remove==1)])
                                
                                if (length(phenos_to_remove)>0) {
                                                phenos_to_remove.phenos.names <- meta_pheno$new_pheno_annot[which(meta_pheno$to_remove==1)][which(meta_pheno$new_pheno_annot[which(meta_pheno$to_remove==1)] %in% colnames(res))]
                                                res.clean <- res[,-phenos_to_remove]
                                } else {
                                                phenos_to_remove.phenos.names <- NA
                                                res.clean <- res
                                }
                                
                                PRS_dim <- dim(res.clean)[2]
                                print(paste0('removing ',length(phenos_to_remove),' acute diseases, except for acute mental or nueroligical phenotypes, and other administrative phenotypes, selected manually.'))
                                
                                # flag multimodal scores - Using dip test
                                
                                if (isTRUE(multimodal)){
                                                print(paste("Flagging scores with multimodal distributions"))
                                                dts <- apply(res.clean,2,dip.test)
                                                dtsp <- lapply(dts,function(x) { x$p.value })
                                                multimodal <- which(dtsp < 0.05/ncol(res.clean))
                                                
                                                print(paste("Removing",length(multimodal),"multimodal scores"))
                                                
                                                if (length(multimodal)>0) {
                                                                multimodal.names <- colnames(res)[multimodal]
                                                                res.clean <- res.clean[,-multimodal]
                                                } else {
                                                                multimodal.names <- NA
                                                                res.clean <- res.clean
                                                }
                                } else {
                                                res.clean <- res.clean
                                                multimodal <- NULL
                                                multimodal.names <- NA
                                }
                                
                                # Use correlation test to remove duplicated or almost duplicated phenotypes
                                
                                matrix <- as.data.frame(res.clean)
                                
                                #cols.cor <- cor(matrix, use = "pairwise.complete.obs", method = "pearson")
                                #cols.cor1 <- fastCor(matrix)
                                cols.cor <- WGCNA::cor(matrix, use = "pairwise.complete.obs", method = "pearson")
                                matrix_dist <- as.matrix(as.dist(cols.cor))
                                
                                #matrix_dist <- as.matrix(dist(t(matrix)))
                                dist_df <- melt_dist(matrix_dist)
                                #dist_df_passed <- dist_df[which(log(dist_df$dist) >= quantile(log(dist_df$dist),0.025)),]
                                check <- dist_df[which(abs(dist_df$dist) > cor_thresh ),]
                                to_remove = NULL
                                if(dim(check)[1] > 0){
                                                for (i in 1:dim(check)[1]){
                                                                corr <- cor.test(as.numeric(unlist(matrix[check[i,1]])),
                                                                                 as.numeric(unlist(matrix[check[i,2]])))
                                                                if((corr$p.value < 0.05/dim(matrix)[2])){
                                                                                if ((grepl( "IC", check[i,1], fixed = TRUE)) & 
                                                                                    (grepl( "PH", check[i,2], fixed = TRUE))){
                                                                                                to_remove <- append(to_remove,check[i,2])
                                                                                } else if ((grepl( "IC", check[i,2], fixed = TRUE)) & 
                                                                                           (grepl( "PH", check[i,1], fixed = TRUE))){
                                                                                                to_remove <- append(to_remove,check[i,1])
                                                                                } else {
                                                                                                to_remove <- append(to_remove,sample(c(check[i,1], check[i,2]), 1))
                                                                                }
                                                                }
                                                }
                                }
                                
                                print(paste0(length(unique(to_remove))," duplicated phenotypes were found"))
                                
                                if (length(unique(to_remove))>0) {
                                                res.clean <- res.clean[ , -which(colnames(res.clean) %in% unique(to_remove))]
                                                
                                } else {
                                                res.clean <- res.clean
                                }

                                # Changing res.clean --> dataframe
                                res.clean <- as.data.frame(res.clean)
                                res.clean$IID <- formod$IID
                                
                                # Save files in list
                                print(paste("Saving data"))
                                
                                allres[[pval]] <- list(residuals=res.clean,
                                               multimodal=multimodal.names,
                                               highlycorrelated=check)
                                
                                row <- c(pval,PRS_dim,low_snp_pheno,length(multimodal),length(to_remove),dim(res.clean)[2]-1)
                                details <- rbind(details,row)

                }
                details <- details[2:dim(details)[1],]
                details <- as.data.frame(details)
                names(details)<-c('P-value threshhold','Number of PRSs',paste0("Number of phenotypes with less than ", snp_count_n," SNPs"), 
                                  "Number of phenotypes with multimodal distribution", 
                                  paste0("Number of duplicated phenotypes with |corr|<",cor_thresh), "Number of PRSs after QC")
                write.csv(details,paste0(path_to_save,'/details.csv'),row.names = FALSE)
                saveRDS(allres,file=paste0(path_to_save,"/Residual_results_all_p-vals.rds"))
                
}

primary_resid_data_GEN(Study="ROSMAP", GenoType="No_MHC", use_snp_count = TRUE, snp_count_n = 5, cor_thresh = 0.95, multimodal=FALSE)
primary_resid_data_GEN(Study="ROSMAP", GenoType="No_APOE", use_snp_count = TRUE, snp_count_n = 5, cor_thresh = 0.95, multimodal=FALSE)
primary_resid_data_GEN(Study="ROSMAP", GenoType="No_MHC_APOE", use_snp_count = TRUE, snp_count_n = 5, cor_thresh = 0.95, multimodal=FALSE)
primary_resid_data_GEN(Study="ROSMAP", GenoType="With_MHC_APOE", use_snp_count = TRUE, snp_count_n = 5, cor_thresh = 0.95, multimodal=FALSE)


