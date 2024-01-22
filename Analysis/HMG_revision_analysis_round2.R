# Intro ----
# To respond to one of the reviewers of the HMG, we are going to apply the following analysis:
  # Q2)  
    # Perform bootstrapping or resampling on the genetic covariance matrix (Π'Π/n) to assess phenotype/module robustness across different samples.
    # Estimate phenotype membership probabilities within each module through repeated sampling.
    # Implement and document this analysis for at least one matrix, considering multiple alpha thresholds.
    # Include findings and implications for the robustness of eigen-PRS methodology in the manuscript revision.
  # Q4)
    # Recalculate R-squared values in Table 2 to account for variance in cognitive decline rate estimates:
      # Compute mean variance of cognitive decline slopes ('c').
      # Identify total residual variance ('b').
      # Apply correction factor 1/(1-c/b) to current R-squared values.
    # Update Table 2 with corrected R-squared values.
    # Briefly discuss the adjusted R-squared impact in the manuscript.


# Libraries ----
GenoTypes <- c("No_APOE","No_MHC","No_MHC_APOE","With_MHC_APOE")
library(tidyverse)
library(ggplot2)
library(gplots) # col2hex
library(WGCNA)
library(flashClust)
library(ggdendro)
library(Hmisc)
library(reshape2)
library(purrr)
library(ggbeeswarm)
library(ggrepel)
library(lme4) # For linear mixed modelling
library(broom.mixed)

# 1. Estimating probability of phenotype's beloning to each module ----
path <- '../Datasets/CLUMP_500_0.2/ROSMAP/Resid_PRS/'
genotype <- GenoTypes[4]
prs_dat <- readRDS(paste0(path,genotype,"/Residual_results_all_p-vals.rds"))

all_wgcna_alpha <- list()
for(alpha in names(prs_dat)[2]){
  dat <- prs_dat[[alpha]]$residuals
  row.names(dat) <- dat$IID
  dat$IID <- NULL
  
  # Now we apply WGCNA with bootstrapping
  all_wgcna <- list()
  for(rand in 1:1000){
    set.seed(rand)
    df <- dat[, sample(ncol(dat), size=as.integer(ncol(dat)*0.95))] # selecting a random subset
    df <- cov(df) # getting the covariance matrix
    powers = c(1:30)
    sft = pickSoftThreshold(df,powerVector=powers,
                            verbose =5,
                            networkType="unsigned",
                            blockSize = ncol(df))
    if(is.na(sft$powerEstimate)) { sft$powerEstimate <- 6 } # the default value by WGCNA
    print(paste0("Threshold= ",sft$powerEstimate))
    print("Running WGCNA")
    sftpower <- sft$powerEstimate
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
    hubdf <- subset(hubdf, module %nin% "grey")
    hexindex <- hubdf$hex
    names(hexindex) <- hubdf$module
    hubdf$hub <- factor(hubdf$hub,levels=hubdf$hub[order(hubdf$n,decreasing = T)])
    # getting all colors
    dat_colour <- data.frame(net$colors)
    dat_colour$phenotype <- row.names(dat_colour)
    row.names(dat_colour) <- NULL
    colnames(dat_colour) <- c("module","phenotype")
    dat_modules <- merge(dat_colour,hubdf,by="module")
    all_wgcna[[rand]] <- dat_modules
  }
  all_wgcna_results <- bind_rows(all_wgcna,.id = "permutation")
  all_wgcna_alpha[[alpha]] <- all_wgcna_results
}

all_wgcna_alpha_results <- bind_rows(all_wgcna_alpha,.id = "alpha_values")
write.csv(all_wgcna_alpha_results,"../Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/ALL/all_wgcna_alpha_results_1e07.csv")
all_wgcna_alpha_results <- read.csv("../Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/ALL/all_wgcna_alpha_results.csv")
# Let's get some statistics from our new dataset:
  # First data set: A cross tabular count between Modules and phenotypes
  # Second data set: A cross tabular count between Modules and phenoype hub

cross_tab_dat1 <- all_wgcna_alpha_results %>%
  filter(alpha_values=="1e-07") %>%
  select(module,phenotype) %>%
  dcast(phenotype~module) %>%
  mutate(newSum = select_if(., is.numeric) %>% 
           reduce(`+`)) %>% 
  mutate_if(is.numeric, list(~ ./newSum)) %>% 
  select(-newSum) #%>%
  # rowwise() %>%
  # mutate(Probability = max(c_across(-phenotype), na.rm = TRUE), # Find the max value
  #        Module = names(.)[-1][which.max(c_across(-phenotype))]) %>% # Find the column name of the max value
  # select(phenotype, Module, Probability) # Select the relevant columns


cross_tab_dat2 <- all_wgcna_alpha_results %>%
  filter(alpha_values=="1e-07") %>%
  select(module,hub,permutation) %>%
  distinct() %>%
  select(-permutation) %>%
  dcast(hub~module) %>%
  mutate(newSum = select_if(., is.numeric) %>% 
           reduce(`+`)) %>% 
  mutate_if(is.numeric, list(~ ./newSum)) %>% 
  select(-newSum) #%>%
  # rowwise() %>%
  # mutate(Probability = max(c_across(-hub), na.rm = TRUE), # Find the max value
  #        Module = names(.)[-1][which.max(c_across(-hub))]) %>% # Find the column name of the max value
  # select(hub, Module, Probability) # Select the relevant columns

# Plotting our results

# Heatmap
dat <- cross_tab_dat1 %>%
  mutate(phenotype=ifelse(phenotype=="LOAD","LOAD",
                          str_to_title(gsub("_"," ",gsub("_h0..*","",
                                                         str_sub(tolower(phenotype),4))))))
row.names(dat) <- dat$phenotype
dat$phenotype <- NULL
heatmap(as.matrix(dat))

# Setting up data for our ggplot
gg_dat <- melt(cross_tab_dat1) %>%
  mutate(category=ifelse(grepl("CO",phenotype),
                         "Continous",
                         ifelse(grepl("PR",phenotype),
                                "Prescription",
                                ifelse(grepl("BI",phenotype),
                                       "Biomarker",
                                       ifelse(grepl("CA",phenotype),
                                              "Categorical",
                                              "ICD")))),
         # Sort phenotype levels alphabetically and reassign them to the phenotype column
         phenotype=factor(phenotype, levels = sort(unique(phenotype),decreasing=T)),
         variable = factor(variable, 
                           levels = c("turquoise","blue","brown",
                                      "black","yellow","pink",
                                      "red","green"),
                           labels = c("turquoise","red","brown",
                                      "black","yellow","pink",
                                      "blue","green"))) %>% # This order was fround from the heatmap
  rename(Module=variable,`Probability Estimate`=value)


gg_grouping <- cross_tab_dat1 %>%
  mutate(groups = cut(seq_len(nrow(.)), breaks = 5, labels = FALSE)) %>%
  select(groups,phenotype)

# merging the groups and the gg_data
gg_dat <- merge(gg_dat,gg_grouping,by="phenotype") %>%
  mutate(phenotype=ifelse(phenotype=="LOAD","LOAD",
                          str_to_title(gsub("pr "," ",gsub("_"," ",gsub("_h0..*","",
                                                                        str_sub(tolower(phenotype),
                                                                                4)))))),
         phenotype=gsub("Nec","NEC",phenotype),
         phenotype=gsub("Dvt","DVT",phenotype),
         phenotype=gsub("Cabg","CABG",phenotype),
         phenotype=gsub("Ptca","PTCA",phenotype),
         phenotype=gsub("Igf1","IGF-1",phenotype),
         phenotype=gsub("Ldl","LDL",phenotype),
         phenotype=gsub("Crit","CRIT",phenotype),
         phenotype=gsub("Ecg","ECG",phenotype),
         phenotype=gsub("Nos","NOS",phenotype),
         phenotype=gsub("Shbg","Sex Hormone-Binding Globulin",phenotype),
         phenotype=gsub("Smoking 1","Smoking",phenotype),
         phenotype=gsub("Sexual Factors 2","Sexual Factors",phenotype),
         phenotype=gsub("Numeric Memory 1","Numeric Memory",phenotype),
         # Sort phenotype levels alphabetically and reassign them to the phenotype column
         phenotype=factor(phenotype, levels = sort(unique(phenotype),decreasing=T)))

# Create the plot with faceting for phenotype and consistent color mapping for variables
p1 <- ggplot(gg_dat %>%
               filter(groups==1),
             aes(x = Module, 
                 y = phenotype, 
                 color = Module)) + 
  geom_line(aes(group=phenotype),
            color="grey") +
  geom_point(aes(size=`Probability Estimate`),shape=16) +
  scale_color_manual(values=c("cyan"="cyan","darkgreen"="darkgreen","orange"="orange","saddlebrown"="saddlebrown","sienna3"="sienna3",
                              "paleturquoise"="paleturquoise","darkolivegreen"="darkolivegreen","palevioletred3","royalblue"="royalblue","brown"="brown","darkred"="darkred",
                              "magenta"="magenta","navajowhite2"="navajowhite2","green"="green","red"="red","pink"="pink","greenyellow"="greenyellow",
                              "midnightblue"="midnightblue","salmon"="salmon","yellow"="yellow","salmon4"="salmon4","skyblue3"="skyblue3","darkturquoise"="darkturquoise",
                              "steelblue"="steelblue","darkgrey"="darkgrey","skyblue"="skyblue","maroon"="maroon","lavenderblush3"="lavenderblush3",
                              "tan"="tan","lightgreen"="lightgreen","black"="black","lightyellow"="lightyellow","floralwhite"="floralwhite","thistle1"="thistle1",
                              "plum2"="plum2","turquoise"="turquoise","blue"="blue","ivory"="ivory","coral1"="coral1","darkseagreen4"="darkseagreen4",
                              "plum1"="plum1","purple"="purple","violet"="violet","bisque4"="bisque4","yellowgreen"="yellowgreen","lightsteelblue1"="lightsteelblue1",
                              "brown4"="brown4","darkorange"="darkorange","grey60"="grey60","lightcyan1"="lightcyan1","honeydew1"="honeydew1","darkorange2"="darkorange2",
                              "mediumpurple3"="mediumpurple3","darkslateblue"="darkslateblue","thistle2"="thistle2","darkmagenta"="darkmagenta","lightpink4"="lightpink4",
                              "orangered4"="orangered4",
                              "white"="lightgrey","lightcyan"="lightcyan3","lightgrey"="cornsilk3")) +  # Map colors to variables
  #facet_grid(phenotype ~ ., scales = "free_y", space = "free_y") +  # Faceting for phenotype with separate scales
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +  # Rotate x-axis labels if needed
  guides(color = "none") +
  ggtitle(latex2exp::TeX("Estimated Probability Of PRSs Belonging To Modules From $\\Pi_{ALL,1$\\times$10^{-7}$}$")) +
  xlab(latex2exp::TeX("")) +
  ylab("")

p2 <- ggplot(gg_dat %>%
               filter(groups==2),
             aes(x = Module, 
                 y = phenotype, 
                 color = Module)) + 
  geom_line(aes(group=phenotype),
            color="grey") +
  geom_point(aes(size=`Probability Estimate`),shape=16) +
  scale_color_manual(values=c("cyan"="cyan","darkgreen"="darkgreen","orange"="orange","saddlebrown"="saddlebrown","sienna3"="sienna3",
                              "paleturquoise"="paleturquoise","darkolivegreen"="darkolivegreen","palevioletred3","royalblue"="royalblue","brown"="brown","darkred"="darkred",
                              "magenta"="magenta","navajowhite2"="navajowhite2","green"="green","red"="red","pink"="pink","greenyellow"="greenyellow",
                              "midnightblue"="midnightblue","salmon"="salmon","yellow"="yellow","salmon4"="salmon4","skyblue3"="skyblue3","darkturquoise"="darkturquoise",
                              "steelblue"="steelblue","darkgrey"="darkgrey","skyblue"="skyblue","maroon"="maroon","lavenderblush3"="lavenderblush3",
                              "tan"="tan","lightgreen"="lightgreen","black"="black","lightyellow"="lightyellow","floralwhite"="floralwhite","thistle1"="thistle1",
                              "plum2"="plum2","turquoise"="turquoise","blue"="blue","ivory"="ivory","coral1"="coral1","darkseagreen4"="darkseagreen4",
                              "plum1"="plum1","purple"="purple","violet"="violet","bisque4"="bisque4","yellowgreen"="yellowgreen","lightsteelblue1"="lightsteelblue1",
                              "brown4"="brown4","darkorange"="darkorange","grey60"="grey60","lightcyan1"="lightcyan1","honeydew1"="honeydew1","darkorange2"="darkorange2",
                              "mediumpurple3"="mediumpurple3","darkslateblue"="darkslateblue","thistle2"="thistle2","darkmagenta"="darkmagenta","lightpink4"="lightpink4",
                              "orangered4"="orangered4",
                              "white"="lightgrey","lightcyan"="lightcyan3","lightgrey"="cornsilk3")) +  # Map colors to variables
  #facet_grid(phenotype ~ ., scales = "free_y", space = "free_y") +  # Faceting for phenotype with separate scales
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +  # Rotate x-axis labels if needed
  guides(color = "none") +
  ggtitle(latex2exp::TeX("Estimated Probability Of PRSs Belonging To Modules From $\\Pi_{ALL,1$\\times$10^{-7}$}$")) +
  xlab(latex2exp::TeX("")) +
  ylab("")

p3 <- ggplot(gg_dat %>%
               filter(groups==3),
             aes(x = Module, 
                 y = phenotype, 
                 color = Module)) + 
  geom_line(aes(group=phenotype),
            color="grey") +
  geom_point(aes(size=`Probability Estimate`),shape=16) +
  scale_color_manual(values=c("cyan"="cyan","darkgreen"="darkgreen","orange"="orange","saddlebrown"="saddlebrown","sienna3"="sienna3",
                              "paleturquoise"="paleturquoise","darkolivegreen"="darkolivegreen","palevioletred3","royalblue"="royalblue","brown"="brown","darkred"="darkred",
                              "magenta"="magenta","navajowhite2"="navajowhite2","green"="green","red"="red","pink"="pink","greenyellow"="greenyellow",
                              "midnightblue"="midnightblue","salmon"="salmon","yellow"="yellow","salmon4"="salmon4","skyblue3"="skyblue3","darkturquoise"="darkturquoise",
                              "steelblue"="steelblue","darkgrey"="darkgrey","skyblue"="skyblue","maroon"="maroon","lavenderblush3"="lavenderblush3",
                              "tan"="tan","lightgreen"="lightgreen","black"="black","lightyellow"="lightyellow","floralwhite"="floralwhite","thistle1"="thistle1",
                              "plum2"="plum2","turquoise"="turquoise","blue"="blue","ivory"="ivory","coral1"="coral1","darkseagreen4"="darkseagreen4",
                              "plum1"="plum1","purple"="purple","violet"="violet","bisque4"="bisque4","yellowgreen"="yellowgreen","lightsteelblue1"="lightsteelblue1",
                              "brown4"="brown4","darkorange"="darkorange","grey60"="grey60","lightcyan1"="lightcyan1","honeydew1"="honeydew1","darkorange2"="darkorange2",
                              "mediumpurple3"="mediumpurple3","darkslateblue"="darkslateblue","thistle2"="thistle2","darkmagenta"="darkmagenta","lightpink4"="lightpink4",
                              "orangered4"="orangered4",
                              "white"="lightgrey","lightcyan"="lightcyan3","lightgrey"="cornsilk3")) +  # Map colors to variables
  #facet_grid(phenotype ~ ., scales = "free_y", space = "free_y") +  # Faceting for phenotype with separate scales
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +  # Rotate x-axis labels if needed
  guides(color = "none") +
  ggtitle(latex2exp::TeX("Estimated Probability Of PRSs Belonging To Modules From $\\Pi_{ALL,1$\\times$10^{-7}$}$")) +
  xlab(latex2exp::TeX("")) +
  ylab("")

p4 <- ggplot(gg_dat %>%
               filter(groups==4),
             aes(x = Module, 
                 y = phenotype, 
                 color = Module)) + 
  geom_line(aes(group=phenotype),
            color="grey") +
  geom_point(aes(size=`Probability Estimate`),shape=16) +
  scale_color_manual(values=c("cyan"="cyan","darkgreen"="darkgreen","orange"="orange","saddlebrown"="saddlebrown","sienna3"="sienna3",
                              "paleturquoise"="paleturquoise","darkolivegreen"="darkolivegreen","palevioletred3","royalblue"="royalblue","brown"="brown","darkred"="darkred",
                              "magenta"="magenta","navajowhite2"="navajowhite2","green"="green","red"="red","pink"="pink","greenyellow"="greenyellow",
                              "midnightblue"="midnightblue","salmon"="salmon","yellow"="yellow","salmon4"="salmon4","skyblue3"="skyblue3","darkturquoise"="darkturquoise",
                              "steelblue"="steelblue","darkgrey"="darkgrey","skyblue"="skyblue","maroon"="maroon","lavenderblush3"="lavenderblush3",
                              "tan"="tan","lightgreen"="lightgreen","black"="black","lightyellow"="lightyellow","floralwhite"="floralwhite","thistle1"="thistle1",
                              "plum2"="plum2","turquoise"="turquoise","blue"="blue","ivory"="ivory","coral1"="coral1","darkseagreen4"="darkseagreen4",
                              "plum1"="plum1","purple"="purple","violet"="violet","bisque4"="bisque4","yellowgreen"="yellowgreen","lightsteelblue1"="lightsteelblue1",
                              "brown4"="brown4","darkorange"="darkorange","grey60"="grey60","lightcyan1"="lightcyan1","honeydew1"="honeydew1","darkorange2"="darkorange2",
                              "mediumpurple3"="mediumpurple3","darkslateblue"="darkslateblue","thistle2"="thistle2","darkmagenta"="darkmagenta","lightpink4"="lightpink4",
                              "orangered4"="orangered4",
                              "white"="lightgrey","lightcyan"="lightcyan3","lightgrey"="cornsilk3")) +  # Map colors to variables
  #facet_grid(phenotype ~ ., scales = "free_y", space = "free_y") +  # Faceting for phenotype with separate scales
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +  # Rotate x-axis labels if needed
  guides(color = "none") +
  ggtitle(latex2exp::TeX("Estimated Probability Of PRSs Belonging To Modules From $\\Pi_{ALL,1$\\times$10^{-7}$}$")) +
  xlab(latex2exp::TeX("")) +
  ylab("")

p5 <- ggplot(gg_dat %>%
               filter(groups==5),
             aes(x = Module, 
                 y = phenotype, 
                 color = Module)) + 
  geom_line(aes(group=phenotype),
            color="grey") +
  geom_point(aes(size=`Probability Estimate`),shape=16) +
  scale_color_manual(values=c("cyan"="cyan","darkgreen"="darkgreen","orange"="orange","saddlebrown"="saddlebrown","sienna3"="sienna3",
                              "paleturquoise"="paleturquoise","darkolivegreen"="darkolivegreen","palevioletred3","royalblue"="royalblue","brown"="brown","darkred"="darkred",
                              "magenta"="magenta","navajowhite2"="navajowhite2","green"="green","red"="red","pink"="pink","greenyellow"="greenyellow",
                              "midnightblue"="midnightblue","salmon"="salmon","yellow"="yellow","salmon4"="salmon4","skyblue3"="skyblue3","darkturquoise"="darkturquoise",
                              "steelblue"="steelblue","darkgrey"="darkgrey","skyblue"="skyblue","maroon"="maroon","lavenderblush3"="lavenderblush3",
                              "tan"="tan","lightgreen"="lightgreen","black"="black","lightyellow"="lightyellow","floralwhite"="floralwhite","thistle1"="thistle1",
                              "plum2"="plum2","turquoise"="turquoise","blue"="blue","ivory"="ivory","coral1"="coral1","darkseagreen4"="darkseagreen4",
                              "plum1"="plum1","purple"="purple","violet"="violet","bisque4"="bisque4","yellowgreen"="yellowgreen","lightsteelblue1"="lightsteelblue1",
                              "brown4"="brown4","darkorange"="darkorange","grey60"="grey60","lightcyan1"="lightcyan1","honeydew1"="honeydew1","darkorange2"="darkorange2",
                              "mediumpurple3"="mediumpurple3","darkslateblue"="darkslateblue","thistle2"="thistle2","darkmagenta"="darkmagenta","lightpink4"="lightpink4",
                              "orangered4"="orangered4",
                              "white"="lightgrey","lightcyan"="lightcyan3","lightgrey"="cornsilk3")) +  # Map colors to variables
  #facet_grid(phenotype ~ ., scales = "free_y", space = "free_y") +  # Faceting for phenotype with separate scales
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +  # Rotate x-axis labels if needed
  guides(color = "none") +
  ggtitle(latex2exp::TeX("Estimated Probability Of PRSs Belonging To Modules From $\\Pi_{ALL,1$\\times$10^{-7}$}$")) +
  xlab(latex2exp::TeX("")) +
  ylab("")

ggsave("../Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/ALL/pheno_p_estimate_1e07_1.jpg",
       p1,
       w=11,h=10, dpi=250)
ggsave("../Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/ALL/pheno_p_estimate_1e07_2.jpg",
       p2,
       w=11,h=10, dpi=250)
ggsave("../Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/ALL/pheno_p_estimate_1e07_3.jpg",
       p3,
       w=11,h=10, dpi=250)
ggsave("../Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/ALL/pheno_p_estimate_1e07_4.jpg",
       p4,
       w=11,h=10, dpi=250)
ggsave("../Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/ALL/pheno_p_estimate_1e07_5.jpg",
       p5,
       w=11,h=10, dpi=250)


# 2. Adjusting the R-squared values for models of predicting cognitive decline ----
# We need to use lmer model on the longitiduinal data of the cognitive decline in ROSMAP particpants
path <- '../Datasets/CLUMP_500_0.2/ROSMAP/Resid_PRS/'
genotype <- GenoTypes[4]
ROSmaster <- readRDS("../Datasets/ROSMAP_Phenotype/ROSmaster.rds")
ROSmaster_long <- read.csv("../Datasets/ROSMAP_Phenotype/dataset_978_long_01-13-2022.csv")
prs_dat <- readRDS(paste0(path,genotype,"/Residual_results_all_p-vals.rds"))
prs_dat <- prs_dat[[1]]$residuals
IIDs <- prs_dat$IID

# merging the two rosmaster dataset and select needed variables
dat_main <- merge(ROSmaster,
                   ROSmaster_long,
                   by=c("projid","study"),
                   all.y=TRUE) %>%
  filter(IID %in% IIDs) %>%
  select(projid,FID,IID,
         study,fu_year,
         msex, educ,age_at_visit,
         cogn_global,
         cogn_global_random_slope,
         cogn_globaln_lv) %>%
  drop_na() %>%
  mutate(IID=as.factor(IID))

# First let's run the lmer model

m <- lmer(cogn_global ~ age_at_visit*fu_year + 
            (age_at_visit + educ) + 
            (1+ fu_year|IID),dat_main)
summary(m)
coef(summary(m))

# Get variance estimates for random effects
var_estimates <- VarCorr(m)

# Print the variance estimates
print(var_estimates)


# Extract random effects
random_effects <- ranef(m, condVar = TRUE)

# Obtain the random slopes for `fu_year` by IID
random_slopes <- as.data.frame(random_effects$IID)
random_slopes$IID <- row.names(random_slopes)
row.names(random_slopes) <- NULL
names(random_slopes) <- c("model_intercept","model_random_slope","IID")

# If you want to incorporate these random slopes into your original data frame
# (assuming the original data frame is 'dat_main' and has a column 'IID')
dat_main_with_slopes <- dat_main %>%
  left_join(random_slopes, by = "IID") %>% 
  mutate(cogn_resids = resid(m)) # adding residuals

# Plotting
# Let's take a quick look at:  

  # 1) model_random_slopes vs cogn_global_random_slope
aa <- dat_main_with_slopes %>% 
  select(IID,cogn_global_random_slope,
         model_random_slope) %>%
  distinct()
cor.test(aa$cogn_global_random_slope,aa$model_random_slope)

# Here x_pos and y_pos represent desired coordinates for the annotation
x_pos <- max(aa$cogn_global_random_slope, na.rm = TRUE) * 0.8
y_pos <- min(aa$model_random_slope, na.rm = TRUE) * 0.8

# Perform the correlation test
cor_test <- cor.test(aa$cogn_global_random_slope, aa$model_random_slope)

# Create a text label with the correlation coefficient and p-value
cor_label <- sprintf("Correlation = %.2f\np-value = %.3f", 
                     cor_test$estimate, cor_test$p.value)

# Create the ggplot
ggplot(aa, aes(x = cogn_global_random_slope, y = model_random_slope)) +
  geom_point() + 
  geom_smooth(method = "lm", color = "blue") +
  annotate("text", x = x_pos, y = y_pos, label = cor_label, 
           hjust = 1, vjust = 1, size = 5, color = "red4") +
  theme_bw() +
  xlab("Random Slope of Global Cognition (Old)") +
  ylab("Random Slope of Global Cognition (New)")

  # 2) cogn_global (last visit) vs cogn_globaln_lv
aa1 <- dat_main_with_slopes %>% 
  group_by(IID) %>%
  filter(fu_year==max(fu_year)) %>%
  select(IID,cogn_global,
         cogn_globaln_lv)
cor.test(aa1$cogn_globaln_lv,aa1$cogn_global)

# Here x_pos and y_pos represent desired coordinates for the annotation
x_pos <- max(aa1$cogn_globaln_lv, na.rm = TRUE) * 0.8
y_pos <- min(aa1$cogn_global, na.rm = TRUE) * 0.8

# Perform the correlation test
cor_test <- cor.test(aa1$cogn_globaln_lv, aa1$cogn_global)

# Create a text label with the correlation coefficient and p-value
cor_label <- sprintf("Correlation = %.2f\np-value = %.3f", 
                     cor_test$estimate, cor_test$p.value)

# Create the ggplot
ggplot(aa1, aes(x = cogn_globaln_lv, y = cogn_global)) +
  geom_point() + 
  geom_smooth(method = "lm", color = "blue") +
  annotate("text", x = x_pos, y = y_pos, label = cor_label, 
           hjust = 1, vjust = 1, size = 5, color = "red4") +
  theme_bw() +
  xlab("Last Visit Global Cognition (Old)") +
  ylab("Last Visit Global Cognition (New)")

  # 3) The cognitive decline estimated slopes:
# Create a sequence of fu_year values for the predictions
fu_year_seq <- seq(min(dat_main_with_slopes$fu_year, na.rm = TRUE), 
                   max(dat_main_with_slopes$fu_year, na.rm = TRUE), 
                   length.out = 100)

# Generate predictions for each IID over the fu_year sequence
by_IID_predictions <- dat_main_with_slopes %>%
  select(IID, model_intercept, model_random_slope) %>%
  distinct() %>%
  group_by(IID) %>%
  do(data.frame(fu_year = fu_year_seq, 
                predicted_cogn_global = .$model_intercept + fu_year_seq * .$model_random_slope)) %>%
  ungroup()

# Plot the observed data points and the predicted lines for each IID
ggplot(data = dat_main_with_slopes, 
       aes(x = fu_year, 
           y = cogn_global, 
           color=IID)) + 
  geom_line(data = by_IID_predictions,
            aes(x = fu_year,
                y = predicted_cogn_global,
                group = IID,
                color=IID),
            color = 'grey',
            alpha = 0.5) +
  geom_point(alpha = 0.5) +
  # geom_line() +
  theme_minimal() +
  guides(color = "none")

# Now, let's extract the statistics to adjust our R^2 values:

# Extract the variance for the random slopes (cognitive decline slopes)
# and the residual variance
var_cov_matrix <- VarCorr(m)
c <- attr(var_cov_matrix, "sc")^2 * var_cov_matrix$IID[["fu_year", "fu_year"]]  # variance of fu_year slope
# Extract the residual standard deviation and square it to get the residual variance
residual_std_dev <- sigma(m)
b <- residual_std_dev^2

