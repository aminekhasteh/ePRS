# Intro ----
# To respond to one of the reviewers of the HMG, we are going to apply the following analysis:
  # Q2)  
    # Perform bootstrapping or resampling on the genetic covariance matrix (Π'Π/n) to assess phenotype/module robustness across different samples.
    # Estimate phenotype membership probabilities within each module through repeated sampling.
    # Implement and document this analysis for at least one matrix, considering multiple alpha thresholds.
    # Include findings and implications for the robustness of eigen-PRS methodology in the manuscript revision.

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
