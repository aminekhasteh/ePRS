# Intro ----
# To respond to one of the reviewers of the HMG, we are going to apply the following analysis:
  # First, we will group all 15 alpha-value thresholds of each trait
  # Second, we apply PCA to each trait (2044 ROSMAP participants and up to 15 columns representing the PRS of the given phenotype at different alpha-value threshold)
  # Third, we take the first PC of each trait as a representative for the given phenotype
  # Lastly, we apply the WGCNA pipeline and apply the rest of our original analysis (Models 1-3)

# Libraries ----
GenoTypes <- c("No_APOE","No_MHC","No_MHC_APOE","With_MHC_APOE")
library(tidyverse)
library(ggplot2)
library(flashClust)
library(WGCNA)
library(ggdendro)
library(gplots)
library(ggbeeswarm)
library(Hmisc)
library(rms)
library(broom)
library(ggrepel)
library(reshape2)

# 1. Grouping all alpha-value thresholds per phenotype ----
path <- '../Datasets/CLUMP_500_0.2/ROSMAP/Resid_PRS/'

meta_data_lst <- list()
for(genotype in GenoTypes){
  prs_dat <- readRDS(paste0(path,genotype,
                            "/Residual_results_all_p-vals.rds"))
  data_lst_tmp <- list()
  for(alpha in names(prs_dat)){
    # getting the main PRS matrices
    tmp_dat <- prs_dat[[alpha]]$residuals
    names(tmp_dat) <- c(paste0(names(tmp_dat)[1:length(names(tmp_dat))-1],
                             "__",alpha), "IID")
    data_lst_tmp[[alpha]] <- tmp_dat
  }
  prs_dat <- data_lst_tmp %>% reduce(full_join, by = "IID")
  meta_data_lst[[genotype]] <- prs_dat
}

# 2. Applying the PC Analysis to each phenotype matrix ----
pc1_PRS_lst <- list()
for(genotype in names(meta_data_lst)){
  prs_dat <- meta_data_lst[[genotype]]
  tmp_lst <- list()
  for(pheno in unique(gsub("\\__.*","",names(prs_dat)))){
    if(pheno=="IID"){next}
    tmp_dat <- prs_dat %>% select(starts_with(pheno))
    # applying PCA:
    pc <- prcomp(tmp_dat)
    tmp_dat <- data.frame(pc$x[,1],prs_dat$IID) %>%
      mutate(pc_explained=sum(summary(pc)$importance[2,1])*100)
    names(tmp_dat) <- c(pheno,"IID",paste0("var_explained_",pheno))
    tmp_lst[[pheno]] <- tmp_dat
  }
  # tmp_dat <- bind_rows(tmp_lst,.id="Phenotype")
  pc1_prs_dat <- tmp_lst %>% reduce(full_join, by = "IID")
  pc1_PRS_lst[[genotype]] <- pc1_prs_dat
}

PC_meta_data <- bind_rows(pc1_PRS_lst,.id = "GenoType")
table(PC_meta_data$GenoType)

# Here we get the information on total variation explained by each phenotype matrix
PC_var_explained <- PC_meta_data %>%
  select(GenoType,starts_with("var_explained_")) %>%
  distinct() %>%
  t() %>%
  as.data.frame()
names <- PC_var_explained[1,]
PC_var_explained <- PC_var_explained[2:dim(PC_var_explained)[1],]
names(PC_var_explained) <- names
PC_var_explained$Phenotypes <- gsub("var_explained_","",row.names(PC_var_explained))
row.names(PC_var_explained) <- NULL

# Getting the PC1 of each phenotype
PC_PRS_dat <- PC_meta_data %>%
  select(!starts_with("var_explained_"))

# Saving
write.csv(PC_var_explained,
          "../Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/ALL/PC_var_explained.csv",
          row.names=FALSE)
write.csv(PC_PRS_dat,
          "../Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/ALL/PC_PRS_dat.csv",
          row.names=FALSE)

# 3. Applying WGCNA ----

all_power_analysis <- list()
for (genotype in names(table(PC_PRS_dat$GenoType))){
  print(genotype)
  df <- PC_PRS_dat %>% 
    filter(GenoType==genotype) %>%
    select_if(~sum(!is.na(.)) > 0) %>%
    select(-GenoType)
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
  all_power_analysis[[genotype]] <- sft
}
names(all_power_analysis) <- names(table(PC_PRS_dat$GenoType))
saveRDS(all_power_analysis,paste0("../Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/ALL/powerstudy_results.rds"))

#all_power_analysis <- readRDS("../Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/ALL/powerstudy_results.rds")
wgcna_results_lst <-  list()
for (genotype in names(table(PC_PRS_dat$GenoType))) {
  print(genotype)
  df <- PC_PRS_dat %>% 
    filter(GenoType==genotype) %>%
    select_if(~sum(!is.na(.)) > 0) %>%
    select(-GenoType)
  rownames(df) <- df$IID
  df$IID <- NULL
  sftpower <- all_power_analysis[[genotype]]$powerEstimate
  if (sftpower >20){sftpower=6}
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
  write.csv(hubdf,paste0("../Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/ALL/",
                         genotype,"_hubdf.csv"),row.names = F)
  
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
  ##### module - trait associations, monocytes & DLPFC
  ROSmaster <- readRDS("../Datasets/ROSMAP_Phenotype/ROSmaster.rds")
  # Reading Filtered PNUKBB manifest
  meta_pheno <- read.csv("./Pan_UKBB/Data/ukbb_manifest_filtered_phenos.csv")
  # Changing Cogdx variable:
  ROSmaster$cogdx[which((ROSmaster$cogdx==2)|(ROSmaster$cogdx==3)|(ROSmaster$cogdx==5)|(ROSmaster$cogdx==6))] <- NA
  ROSmaster$cogdx[which((ROSmaster$cogdx==4))] <- 2
  
  MEs <- net$MEs
  MEs$IID <- rownames(net$MEs)
  md <- merge(MEs,ROSmaster,by="IID",all.x=T)
  # define outcomes and covariates
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
  
  # Heatplot showing associations
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
  write.csv(ressub,paste0("../Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/ALL/",genotype,"_ressub.csv"),row.names = F)
  
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
  
  # sftindex <- sftindex + 1
  
  wgcna_results_lst[[genotype]] <- list(net = net,
                                      assoc_res = results,
                                      associationheatmap = gplt_assoc,
                                      hubPRSbarplot = gplt_hubPRS,
                                      hub = hubdf)
}

saveRDS(wgcna_results_lst,paste0("../Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/ALL/wgcnaresults_all.rds"))

# 4. Plotting the WGCNA results ----

genotype_lst <- NULL
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
prs <- read.csv("../Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/ALL/PC_PRS_dat.csv")
wgcna_dat_results <- readRDS("../Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/ALL/wgcnaresults_all.rds")  
for(genotype in c("No_APOE","No_MHC","No_MHC_APOE","With_MHC_APOE")){
  print(genotype)
  df <- prs %>%
    filter(GenoType == genotype) %>%
    select(-GenoType) %>%
    rename(`PR_5-alpha_reductase_inhibitor_BPH_benign_prostatic_hyperplasia_h0.05708`=`PR_5.alpha_reductase_inhibitor_BPH_benign_prostatic_hyperplasia_h0.05708`,
           `CO_Bone-densitometry_of_heel_right_h0.13979`=`CO_Bone.densitometry_of_heel_right_h0.13979`,
           `CO_Bone-densitometry_of_heel_left_h0.13831`=`CO_Bone.densitometry_of_heel_left_h0.13831`) %>%
    select_if(~sum(!is.na(.)) > 0) 
  rownames(df) <- df$IID
  df$IID <- NULL
  for(color in names(table(wgcna_dat_results[[genotype]]$net$colors))){
    if(color=="grey"){next}
    prs_names <- names(wgcna_dat_results[[genotype]]$net$colors[wgcna_dat_results[[genotype]]$net$colors==color])
    matrix <- df %>% dplyr::select(all_of(prs_names))
    assoc <- wgcna_dat_results[[genotype]]$assoc_res %>%
      filter(egene==paste0("ME",color))
    hub <- wgcna_dat_results[[genotype]]$hub
    for(phen in assoc$pheno){
      tmp_dat <- assoc %>% filter(pheno==phen)
      hub_name <- hub %>% filter(module==color) %>% 
        dplyr::select(hub)
      
      hubs[index] <- gsub("_"," ",gsub("_h0..*","",hub_name$hub)) 
      genotype_lst[index] <- genotype
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

# Adding if LOAD is present in each module
mostly_connected_phenos <- list()
for(genotype in c("No_APOE","No_MHC","No_MHC_APOE","With_MHC_APOE")){
  print(genotype)
  hub <- wgcna_dat_results[[genotype]]$hub
  # find AD scores
  net <- wgcna_dat_results[[genotype]]$net
  ad_col <- net$colors[grep("LOAD",names(net$colors))][1] #  vascular_implants|
  hub$LOAD <- as.numeric(hub$module == ad_col)
  mostly_connected_phenos[[genotype]] <- hub
}
dat_most_connected <- bind_rows(mostly_connected_phenos, .id = "GenoType")

wgcna_pca_dat <- data.frame(genotype_lst,
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
  filter(!Module %in% c("grey")) %>%
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
        by.x=c("geno","Module"),
        by.y=c("GenoType","module")) %>%
  mutate(LOAD=ifelse(LOAD==1,"present","absent"))

p1 <- ggplot(wgcna_pca_dat, 
             aes(x=factor(genotype_lst), 
                 y=abs(psigned_assoc),
                 color=color,shape=LOAD)) + 
  geom_beeswarm(aes(size=ModuleSize)) + 
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
  ggtitle(latex2exp::TeX("Phenotype ~ $\\ePRS_{\\PCA}\\ + Covariates")) +
  labs(x = latex2exp::TeX("Genome Area Inclusion"), 
       y = latex2exp::TeX("-$\\log_{10}(p)$")) +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "tomato3") +
  geom_label_repel(data=subset(wgcna_pca_dat %>% 
                                 group_by(pheno) %>%
                                 mutate(pfdr_assoc = ifelse(pfdr_assoc==min(pfdr_assoc),1,0)), 
                               pfdr_assoc == 1),
                   fill = NA,
                   aes(label=paste(Module,"--",str_sub(tolower(hubs),4))),col="black",size=3,fontface="bold")

ggsave("../Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/ALL/ePRS_assoc.jpg",p1,
       w=10,h=11, dpi=150)
ggsave("../Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/ALL/ePRS_assoc1.jpg",p1,
       w=8,h=8, dpi=150)

# 5. Taking a look at the variations explained by PC1 ----

PC_var_explained <- read.csv("../Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/ALL/PC_var_explained.csv")
pc_var_explained_long <- melt(PC_var_explained)
# row.names(PC_var_explained) <- PC_var_explained$Phenotypes
# PC_var_explained$Phenotypes <- NULL

p2 <- ggplot(pc_var_explained_long, 
          aes(x=factor(variable), 
              y=value)) + 
     geom_quasirandom(color="grey78") +
     theme_bw() +
     theme(plot.title = element_text(hjust = 0.5)) +
     ggtitle(latex2exp::TeX("Total Variation Explained By PRS-PCA For Each Phenotype")) +
     labs(x = latex2exp::TeX("Genome Area Inclusion"), 
          y = latex2exp::TeX("Total Vairaion Explained"))

ggsave("../Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/ALL/PRS-PCA_var_explained.jpg",p2,
       w=9,h=10, dpi=150)
ggsave("../Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/ALL/PRS-PCA_var_explained1.jpg",p2,
       w=8,h=8, dpi=150)

pc_var_explained_long %>% 
  group_by(variable) %>% 
  summarise(mean=mean(value,na.rm=T),
            median=median(value,na.rm=T),
            sd=sd(value,na.rm=T))

# 6. Taking a look at the best performing LOAD PRS across all alpha-values ----
path <- '../Datasets/CLUMP_500_0.2/ROSMAP/Resid_PRS/'
GenoTypes <- c("No_APOE","No_MHC","No_MHC_APOE","With_MHC_APOE")
prs_load_lst <- list()
for(genotype in GenoTypes){
  prs_dat <- readRDS(paste0(path,genotype,
                            "/Residual_results_all_p-vals.rds"))
  data_lst_tmp <- list()
  for(alpha in names(prs_dat)){
    # getting the main PRS matrices
    tmp_dat <- prs_dat[[alpha]]$residuals %>% select(LOAD,IID)
    names(tmp_dat) <- c(paste0(names(tmp_dat)[1:length(names(tmp_dat))-1],
                               "__",alpha), "IID")
    data_lst_tmp[[alpha]] <- tmp_dat
  }
  prs_dat <- data_lst_tmp %>% reduce(full_join, by = "IID")
  prs_load_lst[[genotype]] <- prs_dat
}
prs_load_dat <- bind_rows(prs_load_lst,.id="genotype")

names(prs_load_dat) <- c("genotype","LOAD__5e08","IID","LOAD__1e07","LOAD__5e07",
                         "LOAD__1e06","LOAD__5e06","LOAD__1e05",
                         "LOAD__5e05","LOAD__0_0001","LOAD__0_0005","LOAD__0_001",
                         "LOAD__0_005","LOAD__0_01","LOAD__0_05","LOAD__0_1","LOAD__1")

ROSmaster <- readRDS("../Datasets/ROSMAP_Phenotype/ROSmaster.rds")
# Changing Cogdx variable:
ROSmaster$cogdx[which((ROSmaster$cogdx==2)|(ROSmaster$cogdx==3)|(ROSmaster$cogdx==5)|(ROSmaster$cogdx==6))] <- NA
ROSmaster$cogdx[which((ROSmaster$cogdx==4))] <- 2

# define outcomes and covariates
covars.pathology <- c("msex","pmi","age_death")
covars.cognition <- c("educ","msex") # ,"age_death"
covars.death <- c("msex","age_bl")

indepvec.pathology <- c("cogdx","amyloid_sqrt","tangles_sqrt")
indepvec.cognition <- c("cogn_global_random_slope","cogn_globaln_lv") # Don't have this variable: "cogn_global_at_lastvisit"
pathnames <- c("Final AD","Total AB","PHF tau")
cognames <- c("Global slope","Global last visit") 

varnameindex <- data.frame(name=c(pathnames,cognames),var=c(indepvec.pathology,indepvec.cognition))

nocovars <- FALSE
print("Running Associations")
index <- 1
pvalues <- NULL
bvalues <- NULL
nvalues <- NULL
phenovalues <- NULL
egenevalues <- NULL
genos <- NULL

for(geno in GenoTypes){
  prs_load_dat_geno <- prs_load_dat %>% filter(genotype==geno) %>% select(-genotype)
  md <- merge(prs_load_dat_geno,ROSmaster,by="IID",all.x=T)
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
    
    egenelist <- grep("IID",names(prs_load_dat_geno),value=T,invert = T)
    
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
            form <- as.formula(paste(pheno,"~",egene,"+",paste0(covars,collapse = " + ")))
          }
          mod <- lm(data=md, form)
          pvalues[index] <- tidy(mod)$p.value[2]
          bvalues[index] <- coef(mod)[2]
          nvalues[index] <- dim(mod$model)[1]
          egenevalues[index] <- egene
          phenovalues[index] <- pheno
          genos[index] <- geno
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
          pvalues[index] <- anova(mod)[1,3]
          bvalues[index] <- coef(mod)[2]
          nvalues[index] <- mod$stats[1]
          egenevalues[index] <- egene
          phenovalues[index] <- pheno
          genos[index] <- geno
          index <- index+1
        }
      }
      phenindex <- phenindex + 1
    }
  }
}


results <- data.frame(pheno=phenovalues,
                      egene=egenevalues,
                      genotype=genos,
                      b=bvalues,
                      p=pvalues,
                      n=nvalues,
                      allfdr=p.adjust(pvalues, method="fdr"),
                      signedp=-log10(pvalues)*sign(bvalues)) %>%
  separate(egene,c("LOAD","alpha_value"),sep="__",) %>%
  mutate(alpha_value=gsub('\\_', '.', alpha_value),
         alpha_value=gsub('e', 'e-', alpha_value)) %>%
  merge(.,varnameindex,by.x="pheno",by.y="var",all=T) %>%
  mutate(genotype=ifelse(genotype=="No_APOE","APOE-",
                         ifelse(genotype=="No_MHC","MHC-",
                                ifelse(genotype=="No_MHC_APOE","MHC/APOE-","ALL"))))


p3 <- ggplot(results, 
       aes(x=factor(alpha_value,levels=c("1","0.1","0.05","0.01","0.005","0.001",
                                         "0.0005","0.0001","5e-05","1e-05","5e-06",
                                         "1e-06","5e-07","1e-07","5e-08")), 
           y=signedp, 
           fill=genotype)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Paired") + 
  theme_bw() + 
  facet_wrap(~name,nrow = 4, scales = "free") +
  theme(axis.text.x=element_text(angle = -45, hjust = 0))+
  ylab(latex2exp::TeX("-$\\log_{10}(p)$ $\\times$ sign($\\beta$)")) + #  or $\\R^2$ $\\AUC$
  xlab(latex2exp::TeX("$\\alpha$-value threshold")) +
  ggtitle(latex2exp::TeX("Phenotype ~ $\\PRS_{LOAD,\\alpha}$ + Covariates$")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  guides(fill=guide_legend(title=""))

ggsave("../Datasets/CLUMP_500_0.2/ROSMAP/best_perfoming_PRS_LOAD.jpg",p3,
       w=12,h=8, dpi=150)

# 7. PRSLOAD further detail ----

LOAD_var_count <- list()
for(genotype in GenoTypes[-1]){
  snp_dat <- as.data.frame(fread(paste0("../Datasets/CLUMP_500_0.2/ROSMAP/PRS/",genotype,
                                        "/SNP_count/snp_count.txt"),header=T))
  snp_dat$category <- c("1","0.1","0.05","0.01",
                        "0.005","0.001","0.0005","0.0001",
                        "5e-05","1e-05","5e-06","1e-06",
                        "5e-07","1e-07","5e-08")
  LOAD_var_count[[genotype]] <- snp_dat %>% 
    select(category,LOAD_NEW)
}
LOAD_var_count_dat <- bind_rows(LOAD_var_count,.id="genotype")
write.csv(LOAD_var_count_dat,"../Datasets/CLUMP_500_0.2/ROSMAP/SNP_count_LOAD_NEW.csv")

# 8. PC Analysis of the whole PRS matrix and their association -----
path <- "../Datasets/CLUMP_500_0.2/ROSMAP/PCA/"

indepvec.pathology <- c("cogdx","amyloid_sqrt","tangles_sqrt")
indepvec.cognition <- c("cogn_global_random_slope","cogn_globaln_lv") # Don't have this variable: "cogn_global_at_lastvisit"
pathnames <- c("Final AD","Total AB","PHF tau")
cognames <- c("Global slope","Global last visit") 

varnameindex <- data.frame(name=c(pathnames,cognames),var=c(indepvec.pathology,indepvec.cognition))


for(genotype in GenoTypes){
  dat_pc_assoc <- read.csv(paste0(path,genotype,"/assocres_pc_heatmap_dat.csv"))
  dat_pc_assoc <- dat_pc_assoc %>% 
    mutate(allfdr=p.adjust(p, method="fdr"),
           color=ifelse(allfdr<=0.05,"cornflowerblue","grey78"),
           top_pc_cont = sub("\\--.*", "", PC_contributors),
           LOAD = ifelse(grepl("LOAD",PC_contributors),
                                 "Present","Absent")) %>%
    mutate(alpha_level = as.character(alpha_level),
           alpha_level = ifelse(alpha_level=="5e-04","0.0005",
                                ifelse(alpha_level=="1e-04","0.0001",alpha_level))) %>%
    merge(.,varnameindex,by.x="pheno",by.y="var")
  p <- ggplot(dat_pc_assoc, 
              aes(x=factor(alpha_level,
                           levels=c("1","0.1","0.05","0.01","0.005","0.001",
                                    "0.0005","0.0001","5e-05","1e-05","5e-06",
                                    "1e-06","5e-07","1e-07","5e-08")), 
                  y=abs(signedp),
                  color=color,
                  shape=LOAD)) + 
    geom_quasirandom(aes(color=color)) +
    facet_wrap(~name, scale="free_y",ncol=2) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    guides(color = "none") +
    labs(x = latex2exp::TeX("$\\alpha$-value threshold"), 
         y = latex2exp::TeX("-$\\log_{10}(p)$")) +
    geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "tomato3") +
    geom_label_repel(data=subset(dat_pc_assoc, 
                                 allfdr<=0.05),
                     aes(label=gsub("PR","",paste0(str_sub(gsub("_"," ",sub("\\h0..*", "", PC_contributors)),4)," (",
                                      str_extract(PC_contributors,"\\d+(\\.\\d+){0,1}%"),")"))),
                     col="black",size=3,
                     fontface="bold",
                     fill=NA,
                     max.overlaps=Inf)
  if(genotype=="No_APOE"){
    p <- p+ ggtitle(latex2exp::TeX("Phenotype ~ $\\PC_{APOE-}\\ + Covariates"))
  } else if(genotype=="No_MHC"){
    p <- p+ ggtitle(latex2exp::TeX("Phenotype ~ $\\PC_{MHC-}\\ + Covariates"))
  } else if(genotype=="No_MHC_APOE"){
    p <- p+ ggtitle(latex2exp::TeX("Phenotype ~ $\\PC_{MHC/APOE-}\\ + Covariates"))
  } else{
    p <- p+ ggtitle(latex2exp::TeX("Phenotype ~ $\\PC_{ALL}\\ + Covariates"))
  }
  ggsave(paste0("../Datasets/CLUMP_500_0.2/ROSMAP/PCA/pc_assoc_",genotype,".jpg"),p,
         w=14,h=8, dpi=150)
  
  
}


# 9. Correlation of PRSLOAD and PCA-PRS LOAD ----
# First, we need the PCA-PRS of LOAD
PC_PRS_dat <- read.csv("../Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/ALL/PC_PRS_dat.csv")
PCA_PRS_LOAD_dat <- PC_PRS_dat %>% 
  select(IID,GenoType,starts_with("LOAD")) %>%
  rename(LOAD_PCA=LOAD)

# Second, we need the PRSLOADs
path <- '../Datasets/CLUMP_500_0.2/ROSMAP/Resid_PRS/'
meta_data_lst <- list()
for(genotype in GenoTypes){
  prs_dat <- readRDS(paste0(path,genotype,
                            "/Residual_results_all_p-vals.rds"))
  data_lst_tmp <- list()
  for(alpha in names(prs_dat)){
    # getting the main PRS matrices
    tmp_dat <- prs_dat[[alpha]]$residuals %>% select(LOAD,IID)
    names(tmp_dat) <- c(paste0(names(tmp_dat)[1:length(names(tmp_dat))-1],
                               "__",alpha), "IID")
    data_lst_tmp[[alpha]] <- tmp_dat
  }
  prs_dat <- data_lst_tmp %>% reduce(full_join, by = "IID")
  meta_data_lst[[genotype]] <- prs_dat
}
PRS_LOAD_dat <- bind_rows(meta_data_lst,.id = "GenoType")

# combining the PRS-PCA_LOADs with PRS_LOADs
LOAD_dat <- merge(PCA_PRS_LOAD_dat,PRS_LOAD_dat,
                  by=c("IID","GenoType"))

# Now, we produce three heatmaps

alpha_values <- NULL
corr <- NULL
p_values <- NULL
geno <- NULL
index <- 1
for(genotype in GenoTypes){
  tmp_dat <- LOAD_dat %>% filter(GenoType==genotype) %>% select(-GenoType)
  row.names(tmp_dat) <- tmp_dat$IID
  tmp_dat$IID <- NULL
  LOAD_PC <- tmp_dat$LOAD_PCA
  tmp_dat <- tmp_dat %>% select(-LOAD_PCA)
  for(prs in names(tmp_dat)){
    cortest <- cor.test(LOAD_PC,tmp_dat[[prs]])
    alpha_values[index] <- prs
    corr[index] <- cortest$estimate
    p_values[index] <- cortest$p.value
    geno[index] <- genotype
    index <- index+1
  }
}

cor_test_results <- data.frame(geno,alpha_values,corr,p_values) %>%
  mutate(alpha_values = gsub("LOAD__","",alpha_values),
         geno=ifelse(geno=="No_APOE","APOE-",
                     ifelse(geno=="No_MHC","MHC-",
                            ifelse(geno=="No_MHC_APOE","MHC-APOE-",
                                   "ALL"))),
         cor_sign = ifelse(corr<0,"Negative","Positive")) 
# Plotting:
g <- ggplot(cor_test_results, 
       aes(x = factor(geno), 
           y = factor(alpha_values,levels=c("1","0.1","0.05","0.01","0.005","0.001",
                                            "0.0005","0.0001","5e-05","1e-05","5e-06",
                                            "1e-06","5e-07","1e-07","5e-08")), 
           fill = factor(cor_sign),
           size = abs(corr))) +
  geom_point(shape = 21, stroke = 0) +
  geom_hline(yintercept = seq(.5, 16.5, 1), size = .2) +
  scale_x_discrete(position = "bottom") +
  scale_radius(range = c(1, 15)) +
  # scale_fill_gradient(low = "orange", high = "blue") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        plot.title = element_text(hjust = 0.5)) +
  ylab(latex2exp::TeX("$\\alpha$-value threshold")) + 
  xlab(latex2exp::TeX("Genome Area Inclusion")) +
  guides(size = guide_legend(override.aes = list(fill = NA, color = 
                                                   "black", stroke = .25), 
                             label.position = "bottom",
                             title.position = "top", 
                             order = 1)) +
  labs(size = latex2exp::TeX("Area = |$\\cor(PCAPRS_{LOAD},PRS_{LOAD})|"), 
       fill = "Correlation:",
       title = "Correlation of PRS-PCA of LOAD and PRS of LOAD")

ggsave("../Datasets/CLUMP_500_0.2/ROSMAP/PRS_LOAD_PCAPRS_LOAD_cor.jpg",g,
       w=7,h=8, dpi=200)
 