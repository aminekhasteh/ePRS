# Reading PNUKBB manifest
panuk_manifest <- read.csv("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/code/Pan_UKBB/ukbb_manifest_EUR_h2_05_both_sex_selected_pheno_new_added_annotated.csv")


# Changing the n value at the end of the description:
# panuk_manifest$new_pheno_annot <- paste0(sub("\\._*", '', panuk_manifest$new_pheno_annot) , "_" , panuk_manifest$n_cases_EUR)
# 
# 
# main_data$new_pheno_annot  <- gsub("]", "_", main_data$new_pheno_annot)
# main_data$new_pheno_annot  <- gsub("[[]", "_", main_data$new_pheno_annot)
# 
# write.csv(panuk_manifest,"/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/code/Pan_UKBB/ukbb_manifest_EUR_h2_05_both_sex_selected_pheno_new_added_annotated.csv")


a <- read.csv("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/code/Pan_UKBB/UKBIOBANK/PHECODE_FINAL.csv",sep=";")
a1 <- read.csv("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/code/Pan_UKBB/UKBIOBANK/CATEGORICAL_MEDICATION_FINAL.csv",sep = ";")
a2 <- read.csv("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/code/Pan_UKBB/UKBIOBANK/CATEGORICAL_MEDICAL_FINAL.csv")
a3 <- read.csv("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/code/Pan_UKBB/UKBIOBANK/PRESCRIPTIONS_FINAL.csv")
a4 <- read.csv("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/code/Pan_UKBB/UKBIOBANK/phecode_icd10.csv")

sum(panuk_manifest$aws_link %in% a$aws_link)

aa <- as.data.frame(rbind(a[c("aws_link","AB_ICD10")],a1[c("aws_link","AB_ICD10")],a2[c("aws_link","AB_ICD10")],a3[c("aws_link","AB_ICD10")]))


aaa <- merge(aa,panuk_manifest,by="aws_link")
aaa <-aaa[,-1]
# Selecting only ICD10:

aaaa <- panuk_manifest[which(panuk_manifest$trait_type=="icd10"),]

ICD10_code <- rep("NO_ICD10_code",nrow(panuk_manifest))
ICD10_code[which(panuk_manifest$aws_link %in% aaa$aws_link) ] <- aaa$AB_ICD10
panuk_manifest$ICD10_code <- ICD10_code


panuk_manifest$ICD10_code[which(panuk_manifest$trait_type=="icd10")] <- panuk_manifest$phenocode[which(panuk_manifest$trait_type=="icd10")]
sum(is.na(panuk_manifest$ICD10_code))

##############################


panuk_manifest$new_annot <- NA
for ( i in seq(1,nrow(panuk_manifest))){
                panuk_manifest$new_annot[i] <- paste0(paste(as.vector(str_split(panuk_manifest$new_pheno_annot[i],"_")[[1]][-(length(str_split(panuk_manifest$new_pheno_annot[i],"_")[[1]])-1):-(length(str_split(panuk_manifest$new_pheno_annot[i],"_")[[1]]))]),
                                                            collapse = "_"),"_",panuk_manifest$ICD10_code[i])
                
}

write.csv(panuk_manifest,"/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/code/Pan_UKBB/ukbb_manifest_EUR_h2_05_both_sex_selected_pheno_new_added_annotated1111.csv")
####################
panuk_dat <- read.csv("/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/code/Pan_UKBB/ukbb_manifest_EUR_h2_05_both_sex_selected_pheno_new_added_annotated1111.csv")
panuk_dat$new_annot[which(is.na(panuk_dat$ICD10_code))] <- sub("_[^_]+$", "", panuk_dat$new_annot[which(is.na(panuk_dat$ICD10_code))])

for ( i in seq(1,nrow(panuk_dat))){
                panuk_dat$new_annot[i] <- paste0(gsub("h0","h0.",unlist(str_split(panuk_dat$new_annot[i],"_"))),collapse = "_")
}
write.csv(panuk_dat,"/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/code/Pan_UKBB/ukbb_manifest_EUR_h2_05_both_sex_selected_pheno_new_added_annotated1111.csv")
aa <- panuk_dat[which(is.na(panuk_dat$ICD10_code)),]
write.csv(aa,"/Users/amin/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/code/Pan_UKBB/ukbb_manifest_EUR_h2_05_both_sex_selected_pheno_new_added_annotated2222.csv")



panuk_dat <- read.csv("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/code/Pan_UKBB/ukbb_manifest_EUR_h2_05_both_sex_selected_pheno_new_added_annotated1111.csv")
panuk_dat$new_annot <- sub("_[^_]+$", "", panuk_dat$new_annot)

pheno_remove <- read.table("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/code/Pan_UKBB/phenos_to_remove.txt")
names(pheno_remove) <- "pheno"

sum(panuk_dat$new_pheno_annot %in% pheno_remove$pheno)

write.csv(panuk_dat,"/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/code/Pan_UKBB/ukbb_manifest_EUR_h2_05_both_sex_selected_pheno_new_added_annotated111111111.csv")

l = NULL
for (i in seq(1,nrow(pheno_remove))){
               l <-  append(l,
                            which(grepl(unlist(strsplit(unlist(str_split(pheno_remove$pheno[i],"_"))[length(unlist(str_split(pheno_remove$pheno[i],"_")))-1],"h"))[2], 
                                        panuk_dat$phenocode_annotate_lst)))
}

h = rep(0,nrow(panuk_dat))
for (i in a){
                h[i] <- 1
                
}

panuk_dat$to_remove <- h












######################################################################
library(fuzzywuzzyR)

aaaa <- panuk_manifest[which((panuk_manifest$trait_type == "phecode")& 
                                             (panuk_manifest$ICD10_code=="NO_ICD10_code")),]


a5 <- a4[which(a4$PheCode %in% aaaa$phenocode ),]

a5$ICD10.String + a5$Phenotype



library("stringdist")
d <- expand.grid(aaaa$description,a5$ICD10.String) # Distance matrix in long form
names(d) <- c("a_name","b_name")
d$dist <- stringdist(d$a_name,d$b_name, method="jw") # String edit distance (use your favorite function here)

# Greedy assignment heuristic (Your favorite heuristic here)
greedyAssign <- function(a,b,d){
                x <- numeric(length(a)) # assgn variable: 0 for unassigned but assignable, 
                # 1 for already assigned, -1 for unassigned and unassignable
                while(any(x==0)){
                                min_d <- min(d[x==0]) # identify closest pair, arbitrarily selecting 1st if multiple pairs
                                a_sel <- a[d==min_d & x==0][1] 
                                b_sel <- b[d==min_d & a == a_sel & x==0][1] 
                                x[a==a_sel & b == b_sel] <- 1
                                x[x==0 & (a==a_sel|b==b_sel)] <- -1
                }
                cbind(a=a[x==1],b=b[x==1],d=d[x==1])
}
aa22 <- data.frame(greedyAssign(as.character(d$a_name),as.character(d$b_name),d$dist))
aa33$ICD10 <- a5$ICD10[which(a5$ICD10.String %in% aa33$b)]


panuk_manifest$ICD10_code[which(panuk_manifest$description %in% aa33$a)] <- aa33$ICD10[which(panuk_manifest$description %in% aa33$a)]

for (i in seq(1,nrow(panuk_manifest[which(panuk_manifest$description %in% aa33$a),]))){
                print(aa33$ICD10[which(panuk_manifest[which(panuk_manifest$description %in% aa33$a),][i,"description"] == aa33$a)])
                
                panuk_manifest[which(panuk_manifest$description %in% aa33$a),][i,"ICD10_code"] <- unique(aa33$ICD10[which(panuk_manifest[which(panuk_manifest$description %in% aa33$a),][i,"description"] == aa33$a)])
}

write.csv(panuk_manifest,"/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/code/Pan_UKBB/ukbb_manifest_EUR_h2_05_both_sex_selected_pheno_new_added_annotated11111111345431.csv")
   




##################################

panuk_dat <- read.csv("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/code/Pan_UKBB/ukbb_manifest_EUR_h2_05_both_sex_selected_pheno_new_added_annotated2.csv")
a <- panuk_dat[which((panuk_dat$ICD10_code != "NO_ICD10_code")&
                                (panuk_dat$trait_type != "continuous")&
                                (panuk_dat$trait_type != "biomarkers")),]$new_pheno_annot
b <- sub("_[^_]+$", "", a)


panuk_dat[which((panuk_dat$ICD10_code != "NO_ICD10_code")&
                                (panuk_dat$trait_type != "continuous")&
                                (panuk_dat$trait_type != "biomarkers")),]$new_pheno_annot <- b

grepl("_h0.",b[which(!grepl("h0.0",b))])
panuk_dat$new_pheno_annot <- gsub("__","_",panuk_dat$new_pheno_annot)

write.csv(panuk_dat,"/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/code/Pan_UKBB/ukbb_manifest_EUR_h2_05_both_sex_selected_pheno_new_added_annotated3.csv")



##############################################################
library("readxl")
library(stringr)


panuk1 <- read_excel("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/code/Pan_UKBB/InputData.xlsx", sheet = "All")
panuk_phecode <- read_excel("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/code/Pan_UKBB/InputData.xlsx", sheet = "phecode")
panuk_categorical <- read_excel("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/code/Pan_UKBB/InputData.xlsx", sheet = "categorical")
panuk_prescriptions <- read_excel("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/code/Pan_UKBB/InputData.xlsx", sheet = "prescriptions")


panuk_dat <- rbind(panuk1,panuk_categorical,panuk_phecode,panuk_prescriptions)



for ( i in seq(1,nrow(panuk_dat))){
                if(!grepl("h0.",unlist(str_split(panuk_dat$new_annot[i],"_"))[length(unlist(str_split(panuk_dat$new_annot[i],"_")))])){
                               panuk_dat$new_annot[i] <- paste(unlist(str_split(panuk_dat$new_annot[i],"_"))[1:(length(unlist(str_split(panuk_dat$new_annot[i],"_"))))-1],collapse="_")
                }
}

panuk_dat$ICD10_Automation
str_extract(panuk_dat$ICD10_Automation[1128], '(?<=Code\\s)\\w+')
for (i in seq(1,nrow(panuk_dat[which((panuk_dat$trait_type == "categorical")|
                                     (panuk_dat$trait_type == "phecode")|
                                     (panuk_dat$trait_type == "prescriptions")),]))){
                print(str_extract(panuk_dat$ICD10_Automation[i], '(?<=Code\\s)\\w+.\\w+'))
                
                }

panuk_dat$ICD10_Automation <- str_extract(panuk_dat$ICD10_Automation, '(?<=Code\\s)\\w+.\\w+')


for (i in seq(1,nrow(panuk_dat))){
                panuk_dat$ICD10_Automation[i] <- unlist(str_split(panuk_dat$ICD10_Automation[i]," "))[1]
}




panuk_dat$ICD10_Automation[which(panuk_dat$trait_type == "icd10")] <- panuk_dat[which(panuk_dat$trait_type == "icd10"),]$phenocode

panuk_dat <- read.csv(panuk_dat,"/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/code/Pan_UKBB/ukbb_manifest_EUR_h2_05_both_sex_selected_pheno_new_added_annotated4.csv")
