#install.packages("plotly")
# Libraries
library(ggplot2)
library(plotly)
library(shiny)
library(shinyjs)
library(shinyWidgets)
library(shinydashboard)
library(hrbrthemes)
library(RColorBrewer)
library(tidyverse)
library(latex2exp)
library(factoextra)
library(WGCNA)


# APP #

#header <- dashboardHeader(disable = TRUE)

sidebar <- dashboardSidebar(
                sidebarMenu(
                                selectInput("Genotype",
                                            "Genotype Region",
                                            choices=c("No_APOE","No_MHC",
                                                      "No_MHC_APOE","With_MHC_APOE"), #c("Full Genotype", "MHC region removed", "APOE region removed", "MHC and APOE region removed")
                                            selected = NULL,
                                            multiple = FALSE,
                                            selectize = TRUE,
                                            width = NULL,
                                            size = NULL
                                ),
                                menuItem("Results", icon = icon("list-alt",lib = "glyphicon"), tabName = "Results"),
                                menuItem("Source code", icon = icon("file-code-o"), 
                                         href = "https://github.com/aminekhasteh/Master_Thesis")
                )
                
)

body <- dashboardBody(
                tabItems(
                                # First tab item
                                tabItem(tabName = "Overview",
                                        h2("Association of whole-person polygenic component scores with Alzheimer's disease pathology"),
                                        "Common neurodegenerative disorders, such as Alzheimer’s Disease (AD), 
                                        are highly heterogeneous and genetically correlated with many other 
                                        age-related traits. Growing online compendia of GWAS summary statistics 
                                        now enable phenome-wide assessments of SNP-based polygenic scores (PRS). 
                                        We harnessed such summary statistics to describe the “whole person” 
                                        genetic risk landscape in an elderly population and test if meta-PRS 
                                        scores derived from over 2,000 traits can improve risk prediction of 
                                        AD and AD-related phenotypes. The Religious Order Study and Memory 
                                        Ageing Project (ROS/MAP) are two ongoing longitudinal studies of ageing 
                                        involving annual cognitive testing and post-mortem neuropathological 
                                        assessment. We imputed genotypes for 2,052 unrelated ROS/MAP participants 
                                        of European ancestry using the TopMed Imputation panel. We calculated 
                                        PRSs for up to 2,312 heritable (h2 > 0.05) and sex non-specific phenotypes, 
                                        selected from the Pan-UKBB consortium, at 15 different P-value thresholds 
                                        using PRSice v2.3.3. PRSs were corrected for fine population structure and 
                                        omitted if they were constructed from fewer than five SNPs. Principal 
                                        component analysis was applied to PRSs for all traits at each P-value 
                                        threshold, with the resulting latent components representing “meta-PRS”. 
                                        We observed that the top meta-PRSs explained up to 5.7% of the variability 
                                        in all PRSs at lower p-value thresholds (e.g. 5x10-8) and loaded most 
                                        strongly onto cardiovascular and metabolic phenotypes. At higher P-value 
                                        thresholds, prescriptions for pain killers were the primary loadings of 
                                        top meta-PRSs. In direct phenotypic association analyses, meta-PRSs 
                                        representing all-cause dementia, cholesterol levels, ocular health, 
                                        and cardiovascular disease were most strongly associated with rates 
                                        of cognitive decline and levels of amyloid and tau neuropathology. 
                                        The clustering of “whole person” genetic risk in this population offers 
                                        a map of heterogeneity in genetic risk profiles for AD-related traits. 
                                        It reveals methodological influences on broad-scope polygenic correlation 
                                        in the elderly population."
                                ),
                                tabItem(tabName = "Results",
                                        fluidRow(
                                                        tabBox(
                                                                        tabPanel("PRS Distribution Check",
                                                                                 box(plotOutput('plot0',width = "900px",height = "600px")),
                                                                                 box(
                                                                                                 width = 12,
                                                                                                 title = "Distribution of Resdiuals of the PRSs",
                                                                                                 sliderInput("slider1", "Seed Number:", 1, 10000, 50)
                                                                                 )),
                                                                        
                                                                        tabPanel("Hierarchical Clustering",
                                                                                 box(plotOutput('plot1',width = "1000px",height = "800px")),
                                                                                 box(
                                                                                                 column(
                                                                                                                 width = 12,
                                                                                                                 selectInput("Pval","P-value Threshold",
                                                                                                                             choices=c("0.0001","0.0005","0.001",
                                                                                                                                       "0.005","0.01","0.05","0.1",
                                                                                                                                       "1","1e-05","1e-06","1e-07",
                                                                                                                                       "5e-05","5e-06","5e-07","5e-08"),
                                                                                                                             selected = NULL,
                                                                                                                             multiple = FALSE,
                                                                                                                             selectize = TRUE
                                                                                                 )),
                                                                                                 column(
                                                                                                                 width = 12,
                                                                                                                 selectInput("method","Linkage Method",
                                                                                                                             choices=c("ward", "single", 
                                                                                                                                       "complete", "average", 
                                                                                                                                       "mcquitty", "median", 
                                                                                                                                       "centroid"),
                                                                                                                             selected = NULL,
                                                                                                                             multiple = FALSE,
                                                                                                                             selectize = TRUE
                                                                                                                 ))
                                                                                 )),
                                                                        
                                                                        tabPanel("PCA",
                                                                                 box(
                                                                                                 plotlyOutput('plot2',width = "1000px",height = "800px"))
                                                                                 ),
                                                                        
                                                                        tabPanel("PCA Distribution Check",
                                                                                box(plotOutput('plot3',width = "900px",height = "600px")),
                                                                                box(
                                                                                                width = 12,
                                                                                                title = "Distribution of Loadings from the PCA results",
                                                                                                sliderInput("slider2", "Seed Number:", 1, 10000, 50)
                                                                                )),
                                                                        tabPanel("PCA Association", 
                                                                                 box(plotlyOutput('plot4',width = "1000px",height = "700px")),  #plotOutput("plot1", height = 250)
                                                                                 
                                                                                 box(
                                                                                                 width=12,
                                                                                                 selectInput(
                                                                                                                 "Phenotypes",
                                                                                                                 "Phenotypes",
                                                                                                                 choices=c("amyloid_sqrt","tangles_sqrt",
                                                                                                                           "cogn_global_random_slope",
                                                                                                                           "age_death","parksc_lv_sqrt",
                                                                                                                           "cogdx","gpath",
                                                                                                                           "educ","thyroid_ever",
                                                                                                                           "stroke_ever","diabetes_sr_rx_ever",
                                                                                                                           "cancer_ever","hypertension_ever"),
                                                                                                                 selected = NULL,
                                                                                                                 multiple = FALSE,
                                                                                                                 selectize = TRUE,
                                                                                                                 size = NULL
                                                                                                 )
                                                                                                 

                                                                                 )),
                                                                        tabPanel("Modelling with top PCs", 
                                                                                 box(plotlyOutput('plot5',width = "1000px",height = "700px")),
                                                                                 
                                                                                 box(
                                                                                                 width=12,
                                                                                                 selectInput(
                                                                                                                 "Phenotypes1",
                                                                                                                 "Phenotypes",
                                                                                                                 choices=c("amyloid_sqrt","tangles_sqrt",
                                                                                                                           "cogn_global_random_slope",
                                                                                                                           "age_death","gpath","cogdx"),
                                                                                                                 selected = NULL,
                                                                                                                 multiple = FALSE,
                                                                                                                 selectize = TRUE,
                                                                                                                 size = NULL
                                                                                                 )
                                                                                                 
                                                                                                 
                                                                                 ),
                                                                                 box(
                                                                                                 width=12,
                                                                                                 "Base Model: Phenotype ~ PRS_IGAP_AD + Gender + Baseline Age",
                                                                                                 "PC Model: Phenotype ~ PCi + Gender + Baseline Age",
                                                                                                 "Full Model: Phenotype ~ PRS_IGAP_AD + PCi + Gender + Baseline Age",
                                                                                 )),
                                                                        tabPanel("WGCNA", 
                                                                                 fluidRow(
                                                                                                 column(8, plotOutput("plot6",width = "1000px",height = "400px")),
                                                                                                 column(12, plotOutput("plot7",width = "1000px",height = "400px")),
                                                                                                 column(12, plotOutput("plot8",width = "1000px",height = "400px")),
                                                                                                 column(12, plotOutput("plot9",width = "1000px",height = "400px")),
                                                                                                 box(
                                                                                                                 width=12,
                                                                                                                 selectInput("Pval1",
                                                                                                                             "P-value Threshold",
                                                                                                                             choices=c("0.0001","0.0005","0.001",
                                                                                                                                       "0.005","0.01","0.05","0.1",
                                                                                                                                       "1","1e-05","1e-06","1e-07",
                                                                                                                                       "5e-05","5e-06","5e-07","5e-08"),
                                                                                                                             selected = NULL,
                                                                                                                             multiple = FALSE,
                                                                                                                             selectize = TRUE,
                                                                                                                             size = NULL
                                                                                                 ))
                                                                                 ))
                                                        )
                                                        
                                        )
                                )
                )
)

# Put them together into a dashboardPage
ui <- dashboardPage(
                dashboardHeader(disable = TRUE),
                sidebar,
                body
)

server <- function(input, output) {
                dat_prs_dist <- reactive({
                                path <- paste0("../Datasets/CLUMP_500_0.2/","ROSMAP","/Resid_PRS/",input$Genotype)
                                dat <- readRDS(paste0(path,"/Residual_results_all_p-vals.rds"))
                                set.seed(input$slider1)
                                i <- sample(1:nrow(as.matrix(dat[[input$Pval]]$residuals)),1)
                                hist(as.matrix(dat[[input$Pval]]$residuals[,i]), 
                                     main=paste0("Distribution of Residuals of PRS at \nalpha-value threshold of ",
                                                 input$Pval,
                                                 "\nfor ",colnames(dat[[input$Pval]]$residuals)[i]),
                                     xlab="",col=sample(1:40,1))
                })
                dat_hc <- reactive({
                                path <- paste0("../Datasets/CLUMP_500_0.2/","ROSMAP","/Resid_PRS/",input$Genotype) #
                                dat <- readRDS(paste0(path,"/Residual_results_all_p-vals.rds"))
                                dat[[input$Pval]]$residuals$IID <- NULL
                                par(mar=c(5,5,7,5),cex=1,font=3)
                                par(oma=c(15,3,9,3))
                                heatmap(as.matrix(dat[[input$Pval]]$residuals), hclustfun=function(d) hclust(d, method=input$method))
                                mtext(paste0('Heatmap of PRSs at phenotypes level',"\n",
                                             dim(as.matrix(dat[[input$Pval]]$residuals))[1],' individuals and ',
                                             dim(as.matrix(dat[[input$Pval]]$residuals))[2],' phenotypes',"\n",
                                             'alpha value threshold of ',
                                             input$Pval), side = 3, line = 1, cex = 1, outer=TRUE)
                                
                })
                
                dat_pca <- reactive({
                                path <- paste0("../Datasets/CLUMP_500_0.2/","ROSMAP","/PCA/",input$Genotype)
                                pc_info_dat <- read.csv(paste0(path,"/pc_heatmap_dat.csv"))
                                g <- ggplot(pc_info_dat, aes(reorder(alpha_val,as.numeric(alpha_val)), reorder(pc_name, -as.numeric(pc_vector)),
                                                             fill= as.numeric(pc_vector),
                                                             text=paste0(pc_name," explaining ",pc_vector,"% of the total variation",'\nTop contributors: \n',
                                                                         PC_contributors))) +
                                                geom_tile() +
                                                scale_fill_gradient(low="white", high="red",n.breaks=10) +
                                                theme_ipsum() +
                                                labs(x = "P-value threshold",y="Princiapl Components",fill = "Variation explained")
                                
                                g
                })
                
                dat_pca_dist <- reactive({
                                path <- paste0("../Datasets/CLUMP_500_0.2/","ROSMAP","/PCA/",input$Genotype)
                                dat <- readRDS(paste0(path,"/PCA_results_all_p-vals.rds"))
                                respca <- dat[[input$Pval]]$pc
                                res.var <- get_pca_var(respca)
                                set.seed(input$slider2)
                                i <- sample(1:ncol(respca$x),1)
                                hist(respca$x[,i],
                                     main=paste0("Distribution of PCs for PRS at \nalpha-value threshold of ",
                                                 input$Pval,
                                                 "\nfor ",colnames(respca$x)[i]),
                                     xlab="",col=sample(1:40,1))
                })
                
                dat_pca_assoc <- reactive({
                                path <- paste0("../Datasets/CLUMP_500_0.2/","ROSMAP","/PCA/",input$Genotype)
                                assocres <- read.csv(paste0(path,"/assocres_pc_heatmap_dat.csv"))
                                assocres_pheno <- assocres %>% filter(pheno==input$Phenotypes)
                                assocres_pheno$index <- as.numeric(rownames(assocres_pheno))
                                assocres_pheno$corrected_p <- -log10(p.adjust(assocres_pheno$p,method="fdr"))
                                for(i in seq(1,nrow(assocres_pheno))){
                                                if(assocres_pheno$corrected_p[i]>-log10(0.05)){
                                                                assocres_pheno$psig[i]= "**"
                                                } else if (assocres_pheno$p[i] < 0.05){
                                                                assocres_pheno$psig[i]= "*"
                                                } else {
                                                                assocres_pheno$psig[i]= ""
                                                                
                                                }
                                }

                                # Color Brewer palette
                                assocres_pheno$colour <- ifelse(assocres_pheno$p < 0, "Negative effect","Positive effect")
                                g <- ggplot(assocres_pheno, aes(reorder(alpha_level,as.numeric(alpha_level)),
                                                                reorder(pc,index), fill= as.numeric(corrected_p),
                                                                text=paste0(pc," : uncorrected p-value = ",signif(p,digits = 3),'\n',
                                                                            'Beta coefficient: ',signif(b,digits = 3),'\nTop contributors: \n',
                                                                            'Degrees of freedom: ', n,'\n',
                                                                            PC_contributors))) +
                                                geom_tile() +
                                                geom_text(aes(label = psig)) +
                                                scale_fill_gradient(high="red", low="white",n.breaks=10) +
                                                theme_ipsum() +
                                                labs(x = "SNP selection criteria threshold",y="Princiapl Components",fill = "FDR adjusted p-value (-Log10)",
                                                     title=paste(input$Phenotypes,"~","PCi"))
                                g
                                
                })
                
                dat_pca_compare_assoc <- reactive({
                                path <- paste0("../Datasets/CLUMP_500_0.2/","ROSMAP","/PCA/",input$Genotype)
                                assocres <- read.csv(paste0(path,"/assocres_compare_pc_heatmap_dat.csv"))
                                assocres_pheno <- assocres %>% filter(pheno==input$Phenotypes1) #"cogdx" input$Phenotypes1
                                assocres_pheno$index <- as.numeric(rownames(assocres_pheno))
                                
                                assocres_pheno$corrected_p <- -log10(p.adjust(assocres_pheno$P_likelihood,method="fdr"))
                                for(i in seq(1,nrow(assocres_pheno))){
                                                if(assocres_pheno$corrected_p[i]>-log10(0.05)){
                                                                assocres_pheno$psig[i]= "**"
                                                } else if (assocres_pheno$P_likelihood[i] < 0.05){
                                                                assocres_pheno$psig[i]= "*"
                                                } else {
                                                                assocres_pheno$psig[i]= ""
                                                                
                                                }
                                }

                                # for(i in seq(1,nrow(assocres_pheno))){
                                #                 if(assocres_pheno$r2val_full[i] > assocres_pheno$r2val_base[i]){
                                #                                 assocres_pheno$star[i]= "*"
                                #                 } else {
                                #                                 assocres_pheno$star[i]= ""
                                #                 }
                                # }
                                mean_IGAP_r2 <- NULL
                                for (i in unique(assocres_pheno$alpha_level)){
                                                mean_IGAP_r2[which(assocres_pheno$alpha_level == i)] <- mean(assocres_pheno[which(assocres_pheno$alpha_level == i),]$r2val_base)
                                }
                                
                                # Color Brewer palette
                                if (input$Phenotypes1 == "cogdx"){
                                                g <- ggplot(assocres_pheno, aes(reorder(alpha_level,as.numeric(alpha_level)),
                                                                                reorder(pc,index), fill= as.numeric(r2val_full),
                                                                                text=paste0(pc,'\n',
                                                                                            "Uncorrected p-value = ",signif(P_likelihood,digits = 3),'\n',
                                                                                            'Boot strapped AUC (Base model) = ',signif(mean_IGAP_r2,digits = 3),'\n',
                                                                                            'Boot strapped AUC (PC model) = ',signif(r2val,digits = 3),'\n',
                                                                                            'Boot strapped AUC (Full model) = ',signif(r2val_full,digits = 3),'\n',
                                                                                            'Coefficient of PC = ',signif(coeff_PC,digits = 3), '\n',
                                                                                            'Coefficient of IGAP PRS = ',signif(coeff_IGAP,digits = 3), '\n',
                                                                                            'Degrees of freedom: ', n,
                                                                                            '\nTop contributors: \n',
                                                                                            PC_contributors))) +
                                                                geom_tile() + 
                                                                # geom_point(aes(colour = as.numeric(r2val), 
                                                                #                size = as.numeric(r2val)),
                                                                #            high="red",low="blue") + 
                                                                scale_size(range = c(1, 5)) +
                                                                geom_text(aes(label = psig)) +
                                                                scale_fill_gradient(high="white", low="orange", n.breaks=15) +
                                                                theme_ipsum() +
                                                                labs(x = "SNP selection criteria threshold",
                                                                     y="Princiapl Components",
                                                                     fill = "AUC value of the full model",
                                                                     title= paste0("Y: ",input$Phenotypes1)) + 
                                                                theme(plot.title = element_text(size=5)) #+ 
                                                # annotate("text", x = 1, y = 4.7, size = 3.2,label = str(assocres_pheno$r2val_base[1])) +
                                                # annotate("text", x = 2, y = 26, size = 3.2,label = str(assocres_pheno$r2val_base[2])) +
                                                # annotate("text", x = 3, y = 4.7, size=  3.2, label = str(assocres_pheno$r2val_base[3])) +
                                                # annotate("text", x = 4, y = 4.7, size = 3.2, label = str(assocres_pheno$r2val_base[4]))
                                                
                                }else {
                                                g <- ggplot(assocres_pheno, aes(reorder(alpha_level,as.numeric(alpha_level)),
                                                                                reorder(pc,index), fill= as.numeric(r2val_full),
                                                                                text=paste0(pc,'\n',
                                                                                            "Uncorrected p-value = ",signif(P_likelihood,digits = 3),'\n',
                                                                                            'Boot strapped R2 (Base model) = ',signif(mean_IGAP_r2,digits = 3),'\n',
                                                                                            'Boot strapped R2 (PC model) = ',signif(r2val,digits = 3),'\n',
                                                                                            'Boot strapped R2 (Full model) = ',signif(r2val_full,digits = 3),'\n',
                                                                                            'Coefficient of PC = ',signif(coeff_PC,digits = 3), '\n',
                                                                                            'Coefficient of IGAP PRS = ',signif(coeff_IGAP,digits = 3), '\n',
                                                                                            'Degrees of freedom: ', n,
                                                                                            '\nTop contributors: \n',
                                                                                            PC_contributors))) +
                                                                geom_tile() + 
                                                                # geom_point(aes(colour = as.numeric(r2val), 
                                                                #                size = as.numeric(r2val)),
                                                                #            high="red",low="blue") + 
                                                                scale_size(range = c(1, 5)) +
                                                                geom_text(aes(label = psig)) +
                                                                scale_fill_gradient(high="white", low="orange", n.breaks=15) +
                                                                theme_ipsum() +
                                                                labs(x = "SNP selection criteria threshold",
                                                                     y="Princiapl Components",
                                                                     fill = "Variation explained by the full model",
                                                                     title= paste0("Y: ",input$phenotypes1)) + 
                                                                theme(plot.title = element_text(size=5)) #+ 
                                                # annotate("text", x = 1, y = 4.7, size = 3.2,label = str(assocres_pheno$r2val_base[1])) +
                                                # annotate("text", x = 2, y = 26, size = 3.2,label = str(assocres_pheno$r2val_base[2])) +
                                                # annotate("text", x = 3, y = 4.7, size=  3.2, label = str(assocres_pheno$r2val_base[3])) +
                                                # annotate("text", x = 4, y = 4.7, size = 3.2, label = str(assocres_pheno$r2val_base[4]))
                                }
                                g
                                
                })
                
                dat_wgcna_compare_assoc <- reactive({
                                path <- paste0("../Datasets/CLUMP_500_0.2/","ROSMAP","/WGCNA/",input$Genotype)
                                #path <- "../Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/With_MHC_APOE/"
                                dat <- readRDS(paste0(path,"/wgcnaresults_all.rds"))
                                g <- dat[[input$Pval1]]$associationheatmap
                                g
                                
                })
                
                dat_wgcna_hubs <- reactive({
                                path <- paste0("../Datasets/CLUMP_500_0.2/","ROSMAP","/WGCNA/",input$Genotype)
                                #path <- "../Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/With_MHC_APOE/"
                                dat <- readRDS(paste0(path,"/wgcnaresults_all.rds"))
                                g <- dat[[input$Pval1]]$hubPRSbarplot
                                g
                                
                })
                
                dat_wgcna_clustering <- reactive({
                                path <- paste0("../Datasets/CLUMP_500_0.2/","ROSMAP","/WGCNA/",input$Genotype)
                                #path <- "../Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/With_MHC_APOE/"
                                dat <- readRDS(paste0(path,"/wgcnaresults_all.rds"))
                                net <- dat[[input$Pval1]]$net
                                plotDendroAndColors(net$dendrograms[[1]],
                                                    net$colors,
                                                    "Dynamic Tree Cut",
                                                    dendroLabels = FALSE,
                                                    hang = 0,
                                                    #addGuide = TRUE,
                                                    guideHang = 0.05,
                                                    main = "PRS dendrogram and module colors")
                                
                })
                
                dat_wgcna_hubs_clustering <- reactive({
                                path <- paste0("../Datasets/CLUMP_500_0.2/","ROSMAP","/WGCNA/",input$Genotype)
                                #path <- "../Thesis_Project/Datasets/CLUMP_500_0.2/ROSMAP/WGCNA/With_MHC_APOE/"
                                dat <- readRDS(paste0(path,"/wgcnaresults_all.rds"))
                                net <- dat[[input$Pval1]]$net
                                plot(hclust(dist(t(net$MEs))),xlab="",ylab="Colours")
                                
                })
                
                
                output$plot0<-renderPlot({dat_prs_dist()})
                output$plot1<-renderPlot({dat_hc()})
                output$plot2<-renderPlotly({
                                ggplotly(dat_pca(), tooltip="text")
                })
                output$plot3<-renderPlot({dat_pca_dist()})
                output$plot4<-renderPlotly({
                                ggplotly(dat_pca_assoc(), tooltip="text")})
                output$plot5<-renderPlotly({
                                ggplotly(dat_pca_compare_assoc(), tooltip="text")})
                output$plot6<-renderPlot({dat_wgcna_compare_assoc()})
                output$plot7<-renderPlot({dat_wgcna_hubs()})
                output$plot8<-renderPlot({dat_wgcna_clustering()})
                output$plot9<-renderPlot({dat_wgcna_hubs_clustering()})
                
}




shinyApp(ui, server)



