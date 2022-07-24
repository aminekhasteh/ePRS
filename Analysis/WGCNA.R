# Libraries ----
library(WGCNA)
library(rms)
library(limma)
library(flashClust)
library(ggplot2)
library(cowplot)
library(gplots)
library(ggdendro)
library(ggraph)
library(dplyr)
library(broom)
library(parallel)
library(doParallel)
library(latex2exp)
library(ggbeeswarm)
library(shades)
library(colorspace)

# 1) Soft thresholding measurement ----
  # This function generates initial power analysis to determine softpowers for each alpha-value threshold
WGCNA_sftpowers_GEN <- function(Study=c("ROSMAP","ADNI"),
                                GenoType=c("No_APOE","No_MHC",
                                           "No_MHC_APOE","With_MHC_APOE")){
                path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/Resid_PRS/",GenoType)
                path_to_save <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/WGCNA/",GenoType)
                prs <- readRDS(paste0(path,"/Residual_results_all_p-vals.rds"))
                all_power_analysis <- list()
                for (pthres in names(prs)[order(as.numeric(names(prs)),decreasing = F)]){
                                print(pthres)
                                df <- prs[[pthres]]$residuals
                                rownames(df) <- df$IID
                                df$IID <- NULL
                                # sftpower <- sftpowers[sftindex]
                                # cluster subjects and plot dendrogram
                                sampleTree1 <- flashClust(dist(df),method = "ward")
                                ggdendrogram(sampleTree1) + labs(title=paste0("n=",nrow(df)))
                                ### cluster PRS
                                PRSTree1 <- flashClust(dist(t(df)),method = "ward")
                                dend1 <- ggdendrogram(PRSTree1)+
                                                labs(title=paste0("n=",ncol(df)))
                                # determine optimal power visually - here set at 11
                                print("Choosing Soft Power Threshold")
                                powers = c(1:30)
                                sft = pickSoftThreshold(df,powerVector=powers,
                                                        verbose =5,
                                                        networkType="unsigned",
                                                        blockSize = ncol(df))
                                if(is.na(sft$powerEstimate)) { sft$powerEstimate <- 1 }
                                print(paste0("Threshold= ",sft$powerEstimate))
                                # p1 <- ggplot(data=sft$fitIndices,aes(y=SFT.R.sq,x=Power))+
                                #                 geom_hline(yintercept=c(0.85,0.9),lty=2,col="red")+
                                #                 geom_line(size=1.5,alpha = 0.2, stat = "smooth", method = "loess",col="blue")+
                                #                 geom_text(aes(label=Power),size=3)+
                                #                 geom_vline(xintercept=16,lty=1,col="darkgreen")+
                                #                 labs(y="Signed R2")+
                                #                 theme_minimal()+
                                #                 ggtitle("Scale Independence")
                                # 
                                # p2 <- ggplot(data=sft$fitIndices,aes(y=mean.k.,x=Power))+
                                #                 geom_line(size=1.5,alpha = 0.2, stat = "smooth", method = "loess",col="blue")+
                                #                 geom_text(aes(label=Power),size=3)+
                                #                 geom_vline(xintercept=16,lty=1,col="darkgreen")+
                                #                 labs(y="Mean connectivity")+
                                #                 theme_minimal()+
                                #                 ggtitle("Mean connectivity")
                                # p <- plot_grid(p1, p2, labels=c('A', 'B'))
                                # title <- ggdraw() + draw_label(TeX(paste0("$\\alpha$-value threshold of ",pthres)), fontface='bold')
                                # pdf(paste0(path_to_save,"/powerstudy_",pthres,".pdf"),w=10,h=8)
                                # print(plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)))
                                # dev.off()
                                all_power_analysis[[pthres]] <- sft
                }
                saveRDS(all_power_analysis,paste0(path_to_save,"/powerstudy_results.rds"))
}

WGCNA_sftpowers_GEN(Study="ROSMAP",GenoType="With_MHC_APOE")
WGCNA_sftpowers_GEN(Study="ROSMAP",GenoType="No_APOE")
WGCNA_sftpowers_GEN(Study="ROSMAP",GenoType="No_MHC")
WGCNA_sftpowers_GEN(Study="ROSMAP",GenoType="No_MHC_APOE")

# 2) Network generation ----
  # This function creates detailed hub analysis from the WGCNA package, also creates initial association results
WGCNA_dat_plot_GEN <- function(sftpowers,
                               Study,
                               GenoType){
                path <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/Resid_PRS/",GenoType)
                path_to_save <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/",Study,"/WGCNA/",GenoType)
                prs <- readRDS(paste0(path,"/Residual_results_all_p-vals.rds"))
                sftindex <- 1
                wgcna_results_lst <-  list()
                for (pthres in names(prs)[order(as.numeric(names(prs)),decreasing = F)]) {
                                df <- prs[[pthres]]$residuals
                                rownames(df) <- df$IID
                                df$IID <- NULL
                                sftpower <- sftpowers[sftindex]
                                print("Running WGCNA")
                                net <- blockwiseModules(datExpr = df,
                                                        power=sftpower,
                                                        TOMType = "unsigned",
                                                        minModuleSize = 10,
                                                        minKMEtoStay = 0.01,
                                                        maxBlockSize = 20000,
                                                        deepSplit = 4,
                                                        detectCutHeight = 0.999,
                                                        saveTOMs = F,
                                                        pamStage = T,
                                                        pamRespectsDendro = T,
                                                        mergeCutHeight = 0.1,
                                                        networkType = "unsigned",
                                                        verbose = 10,
                                                        minCoreKME = 0,
                                                        minModuleKMESize = 1)
                                
                                # gplt_dendrogram <- plotDendroAndColors(net$dendrograms[[1]],
                                #                                        net$colors,
                                #                                        "Dynamic Tree Cut",
                                #                                        dendroLabels = FALSE,
                                #                                        hang = 0,
                                #                                        #addGuide = TRUE,
                                #                                        guideHang = 0.05,
                                #                                        main = "Gene dendrogram and module colors")
                                # print("Plotting Dendrogram")
                                # pdf(paste0(path_to_save,"/dendrogram_",pthres,".pdf"),w=10,h=7)
                                # print(gplt_dendrogram)
                                # dev.off()
                                
                                print("Getting Hub PRS")
                                #omitColors = ifelse(sum(table(net$colors) < 2) > 0,
                                #                    names(which(table(net$colors) < 2)),
                                #                    NA) # Omitting colors with less than 2 number of PRSs, otherwise this function won't work.
                                omitColors = names(table(net$colors))[which(as.numeric(table(net$colors)) < 2)]
                                if(identical(names(table(net$colors))[which(as.numeric(table(net$colors)) < 2)], character(0))){
                                                omitColors=F
                                }
                                hubprs <- chooseTopHubInEachModule(datExpr = df,
                                                                   colorh = net$colors,
                                                                   power = sftpower,
                                                                   omitColors = omitColors,
                                                                   type = "unsigned")
                                
                                hubdf <- data.frame(module=names(hubprs),hub=hubprs)
                                hubdf$n <- as.numeric(table(net$colors))[which(as.numeric(table(net$colors)) >= 2)]
                                hubdf$hex <- col2hex(hubdf$module)
                                hubdf$value <- 1
                                hubdf <- subset(hubdf, module %nin% "grey")
                                
                                hexindex <- hubdf$hex
                                names(hexindex) <- hubdf$module
                                
                                hubdf$hub <- factor(hubdf$hub,levels=hubdf$hub[order(hubdf$n,decreasing = T)])
                                write.csv(hubdf,paste0(path_to_save,"/",pthres,"_hubdf.csv"),row.names = F)
                                
                                gplt_hubPRS <- ggplot(data=hubdf,aes(x=n,y=hub,fill=module))+
                                                geom_bar(stat="identity",show.legend = F)+
                                                geom_text(aes(label=module),nudge_x = 5)+
                                                scale_fill_manual(values=hexindex)+
                                                labs(y="Hub gene",x="Number of PRS")+
                                                theme_classic()
                                
                                # pdf(paste0(path_to_save,"/hubPRSbarplot_",pthres,".pdf"),w=10,h=10)
                                # print(gplt_hubPRS)
                                # dev.off()
                                
                                # gplt_moduleclusters <- plot(hclust(dist(t(net$MEs))))
                                
                                # dendrogram of modules
                                # pdf(paste0(path_to_save,"/moduleclusters_",pthres,".pdf"),w=10,h=8)
                                # print()
                                # dev.off()
                                
                                ##### association with AD phenotypes
                                ###################################################
                                ##### module - trait associations, monocytes & DLPFC
                                ROSmaster <- readRDS("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/ROSMAP_Phenotype/ROSmaster.rds")
                                # Reading Filtered PNUKBB manifest
                                meta_pheno <- read.csv("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/code/Pan_UKBB/ukbb_manifest_filtered_phenos.csv")
                                # Changing Cogdx variable:
                                ROSmaster$cogdx[which((ROSmaster$cogdx==2)|(ROSmaster$cogdx==3)|(ROSmaster$cogdx==5)|(ROSmaster$cogdx==6))] <- NA
                                ROSmaster$cogdx[which((ROSmaster$cogdx==4))] <- 2
                                
                                MEs <- net$MEs
                                MEs$IID <- rownames(net$MEs)
                                md <- merge(MEs,ROSmaster,by="IID",all.x=T)
                                ####################################
                                ####################################
                                #### define outcomes and covariates
                                ####################################
                                covars.pathology <- c("msex","pmi","age_death")
                                covars.cognition <- c("educ","msex") # ,"age_death"
                                covars.death <- c("msex","age_bl")
                                
                                # indepvec.pathology <- c("plaq_n_sqrt","plaq_d_sqrt","amyloid_sqrt","tangles_sqrt","nft_sqrt",
                                #                         "ci_num2_gct","ci_num2_mct","arteriol_scler","caa_4gp","cvda_4gp2",
                                #                         "dlbdx","hspath_any","tdp_stage4","parkdx","pathoAD","vm3123",
                                #                         "pput3123","it3123","mf3123")
                                
                                # indepvec.pathology <- c("plaq_n_sqrt","plaq_d_sqrt","cogdx","amyloid_sqrt","tangles_sqrt",
                                #                         "pathoAD","nft_sqrt","dlbdx")
                                # indepvec.cognition <- c("cogn_global_random_slope","cogn_globaln_lv") # Don't have this variable: "cogn_global_at_lastvisit"
                                # indepvec.death <- c("age_death")
                                # pathnames <- c("Neuritic plaques","Diffuse plaques","Final AD","Total AB","PHF tau",
                                #                "Patho AD","NFT","Lewy body stage")
                                # # pathnames <- c("Neuritic plaques","Diffuse plaques","Total AB","PHF tau","NFT","Gross cerebral infarcts",
                                # #                "Micro cerebral infarcts","Arteriolosclerosis","Cerebral AA","Cerebral atherosclerosis",
                                # #                "Lewy body stage","Hippocampal sclerosis","TDP-43","PD Dx","Patho AD","PAM VM Caudate",
                                # #                "PAM post. putamen","PAM IT","PAM MF")
                                # cognames <- c("Global slope","Global last visit")
                                # deathnames <- c("Age of death")
                                
                                
                                indepvec.pathology <- c("cogdx","amyloid_sqrt","tangles_sqrt")
                                indepvec.cognition <- c("cogn_global_random_slope","cogn_globaln_lv") # Don't have this variable: "cogn_global_at_lastvisit"
                                pathnames <- c("Final AD","Total AB","PHF tau")
                                cognames <- c("Global slope","Global last visit") 
                                
                                
                                varnameindex <- data.frame(name=c(pathnames,cognames),var=c(indepvec.pathology,indepvec.cognition))
                                
                                nocovars <- FALSE
                                #######################################
                                print("Running Associations")
                                index <- 1
                                pvalues <- NULL
                                bvalues <- NULL
                                nvalues <- NULL
                                phenovalues <- NULL
                                egenevalues <- NULL
                                egen_contributors <- NULL
                                hub_contributors <- NULL
                                for(vartype in c("path","cog")) {
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
                                                egenelist <- grep("IID",names(MEs),value=T,invert = T)
                                                
                                                if (length(indepvec)>1){
                                                                tablist <- apply(md[,indepvec],2,table) %>%
                                                                                lapply(length) %>%
                                                                                as.numeric()
                                                } else {
                                                                tablist <- table(md[,indepvec])  %>%
                                                                                lapply(length) %>%
                                                                                as.numeric() %>% sum()      
                                                }
                                                
                                                cont.index <- tablist > 2 #TRUE if numeric indepvec
                                                phenindex <- 1
                                                for (pheno in indepvec) {
                                                                if (cont.index[phenindex]==T) {
                                                                                for (egene in egenelist) {
                                                                                                if (nocovars==TRUE) {
                                                                                                                form <- formula(paste(pheno,"~",egene))
                                                                                                } else {
                                                                                                                form <- formula(paste(pheno,"~",egene,"+",paste0(covars,collapse = " + ")))
                                                                                                }
                                                                                                mod <- lm(data=md, form)
                                                                                                hub_contributors[index] <- paste0(paste0(names(net$colors[which(net$colors == substring(egene, 3))]),
                                                                                                                                         " -- ICD10 Code: ",
                                                                                                                                         meta_pheno[match(names(net$colors[which(net$colors == substring(egene, 3))]), meta_pheno$new_pheno_annot),]$ICD10,
                                                                                                                                         collapse="\n"))
                                                                                                pvalues[index] <- tidy(mod)$p.value[2]
                                                                                                bvalues[index] <- coef(mod)[2]
                                                                                                nvalues[index] <- dim(mod$model)[1]
                                                                                                egenevalues[index] <- egene
                                                                                                phenovalues[index] <- pheno
                                                                                                egen_contributors[index] <- paste(names(net$colors[which(net$colors == egene)]),collapse = "\n")
                                                                                                index <- index+1
                                                                                }
                                                                } else {
                                                                                for (egene in egenelist) {
                                                                                                if (nocovars==TRUE) {
                                                                                                                form <- formula(paste(pheno,"~",egene))
                                                                                                } else {
                                                                                                                form <- formula(paste(pheno,"~",egene,"+",paste0(covars,collapse = " + ")))
                                                                                                }
                                                                                                mod <- lrm(data=md, form)
                                                                                                hub_contributors[index] <- paste0(paste0(names(net$colors[which(net$colors == substring(egene, 3))]),
                                                                                                                                         " -- ICD10 Code: ",
                                                                                                                                         meta_pheno[match(names(net$colors[which(net$colors == substring(egene, 3))]), meta_pheno$new_pheno_annot),]$ICD10,
                                                                                                                                         collapse="\n"))
                                                                                                pvalues[index] <- anova(mod)[1,3]
                                                                                                bvalues[index] <- coef(mod)[2]
                                                                                                nvalues[index] <- mod$stats[1]
                                                                                                egenevalues[index] <- egene
                                                                                                phenovalues[index] <- pheno
                                                                                                egen_contributors[index] <- paste(names(net$colors[which(net$colors == egene)]),collapse = "\n")
                                                                                                index <- index+1
                                                                                }
                                                                }
                                                                phenindex <- phenindex + 1
                                                }
                                }
                                
                                results <- data.frame(pheno=phenovalues,
                                                      egene=egenevalues,
                                                      egene_contributors = hub_contributors,
                                                      b=bvalues,
                                                      p=pvalues,
                                                      n=nvalues,
                                                      allfdr=p.adjust(pvalues, method="fdr"),
                                                      signedp=-log10(pvalues)*sign(bvalues))
                                
                                #results <- subset(results, pheno %in% c("pathoAD","amyloid_sqrt","tangles_sqrt","cogn_global_at_lastvisit","cogn_global_slope"))
                                #results$fdr <- p.adjust(results$p)
                                results$fdr <- results$allfdr
                                
                                ######################################
                                ###### Heatplot showing associations
                                print("Plotting Heatmap")
                                
                                ressub <- results
                                ressub$pheno <- varnameindex$name[match(ressub$pheno,varnameindex$var)]
                                ressub$egene <- gsub(".*ME","",ressub$egene)
                                
                                moddefs <- data.frame(module=names(table(net$colors)),
                                                      ngenes=as.character(table(net$colors)))
                                
                                ressub$modn <- moddefs$ngenes[match(ressub$egene,moddefs$module)]
                                ressub$egene <- paste0(ressub$egene," (",ressub$modn,")")
                                
                                # bonfT <- 0.05/nrow(ressub)
                                rescast <- reshape2::dcast(ressub, pheno ~ egene,value.var = "signedp")
                                
                                rownames(rescast) <- rescast$pheno
                                rescast$pheno <- NULL
                                rescast <- as.matrix(rescast)
                                
                                egenedend <- as.dendrogram(hclust(dist(t(rescast))))
                                egene.order <- order.dendrogram(egenedend)
                                
                                ressub$egene.f <- factor(x = ressub$egene,
                                                         levels = colnames(rescast)[egene.order], 
                                                         ordered = TRUE)
                                
                                ressub$pheno.f <- factor(x = ressub$pheno,
                                                         levels = c(pathnames,cognames), 
                                                         ordered = TRUE)
                                write.csv(ressub,paste0(path_to_save,"/",pthres,"_ressub.csv"),row.names = F)
                                
                                gplt_assoc <- ggplot(data=ressub,aes(y=pheno.f,x=egene.f,
                                                                     fill=-log10(p)*sign(b)))+
                                                geom_tile(show.legend = T)+
                                                scale_fill_gradient2(low="cornflowerblue",high="brown1",mid = "white")+
                                                geom_text(data=subset(ressub, p < 0.05 & fdr > 0.05),aes(label=signif(p,2)),col="black",size=2)+
                                                geom_text(data=subset(ressub, fdr < 0.05),aes(label=signif(p,2)),col="black",size=3,fontface="bold")+
                                                labs(x="Phenotype",y="PRS Module",title="PRS module-trait associations")+
                                                theme_minimal()+
                                                theme(axis.text.x=element_text(angle = -45, hjust = 0))
                                # pdf(paste0(path_to_save,"/associationheatmap_",pthres,".pdf"),w=16,h=10)
                                # print(gplt_assoc)
                                # dev.off()
                                
                                sftindex <- sftindex + 1
                                
                                wgcna_results_lst[[pthres]] <- list(net = net,
                                                                    assoc_res = results,
                                                                    associationheatmap = gplt_assoc,
                                                                    hubPRSbarplot = gplt_hubPRS,
                                                                    hub = hubdf)
                }
                saveRDS(wgcna_results_lst,file=paste0(path_to_save,"/wgcnaresults_all",".rds"))
}

WGCNA_dat_plot_GEN(sftpowers = c(3,2,17,1,1,1,1,1,1,1,1,2,1,1,1),Study="ROSMAP",GenoType="With_MHC_APOE")
WGCNA_dat_plot_GEN(sftpowers = c(8,7,4,1,1,1,1,1,1,1,1,1,1,1,1),Study="ROSMAP",GenoType="No_MHC")
WGCNA_dat_plot_GEN(sftpowers = c(4,2,3,1,1,1,1,1,1,1,1,2,1,1,1),Study="ROSMAP",GenoType="No_APOE")
WGCNA_dat_plot_GEN(sftpowers = c(10,4,1,1,1,1,1,1,1,1,1,1,1,1,1),Study="ROSMAP",GenoType="No_MHC_APOE")
# 3) Module information ----
  # Getting module information (number of modules for each phenotype)
PRSs <- NULL
LEN_MODULE <- NULL
Pi <- NULL
index<-1
path_wgcna <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/","ROSMAP","/WGCNA/","No_MHC")
wgcna_dat_results <- readRDS(file=paste0(path_wgcna,"/wgcnaresults_all",".rds"))
for(prs in names(wgcna_dat_results)){
                for(GenoType in c("No_APOE","No_MHC","No_MHC_APOE","With_MHC_APOE")){
                                path_wgcna <- paste0("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/","ROSMAP","/WGCNA/",GenoType)
                                wgcna_dat_results <- readRDS(file=paste0(path_wgcna,"/wgcnaresults_all",".rds"))   
                                colors <- wgcna_dat_results[[prs]]$net$colors
                                Pi[index] <- GenoType
                                PRSs[index] <- prs
                                LEN_MODULE[index] <- length(table(wgcna_dat_results[[prs]]$net$colors))
                                index <- index+1
                }
}

module_count <- as.data.frame(cbind(Pi,PRSs,LEN_MODULE))
module_count <- reshape2::dcast(module_count,PRSs~Pi,value.var='LEN_MODULE')
write.csv(module_count,'/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/module_counts.csv',row.names = F)
                     