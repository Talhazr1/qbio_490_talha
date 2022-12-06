getwd()
knitr::opts_knit$set(root.dir = normalizePath("../final_project_group#4/outputs")) 
if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
library("ggplot2")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library("BiocManager")
if (!require("TCGAbiolinks", quietly = TRUE))
  install.packages("TCGAbiolinks")
library("TCGAbiolinks")
if (!require("maftools", quietly = TRUE))
  install.packages("maftools")
library("maftools")
if (!require("SummarizedExperiment", quietly = TRUE))
  install.packages("SummarizedExperiment")
library("SummarizedExperiment") #use for future project
if (!require("survival", quietly = TRUE))
  install.packages("survival")
library("survival")
if (!require("survminer", quietly = TRUE))
  install.packages("survminer")
library("survminer")
if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
library("dplyr")
clinical_query <- GDCquery(project="TCGA-GBM", data.category = "Clinical", file.type = "xml")
clinical_data <- GDCdownload(clinical_query)
clinical <- GDCprepare_clinic(query = clinical_query, clinical.info = "patient")
clinical.drug <- GDCprepare_clinic(query = clinical_query, clinical.info = "drug")
maf_query <- GDCquery(
  project = "TCGA-GBM", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", 
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(maf_query)
maf <- GDCprepare(maf_query) 
colnames(clinical)[ colnames(clinical) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode" #rename barcode column for read.maf()
maf_object <- read.maf(maf = maf, 
                       clinicalData = clinical,
                       isTCGA = TRUE)
rna_query <- GDCquery(project = "TCGA-GBM",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")
GDCdownload(rna_query)
rna_se <- GDCprepare(rna_query)
maf_object@clinical.data$age_class <- ifelse(maf_object@clinical.data$age_at_initial_pathologic_diagnosis > 50, "Old", "Young")
clinical$age_class <- ifelse(clinical$age_at_initial_pathologic_diagnosis > 50, "Old", "Young")
young_mask = ifelse(clinical$age_class == "Young", TRUE, FALSE)
old_mask = ifelse(clinical$age_class == "Old", TRUE, FALSE)
young_clinical <- clinical[young_mask,]
old_clinical <- clinical[old_mask, ]
gender_clinical <- clinical[ifelse(clinical$gender == " ", F, T), ]
jpeg(filename = "../final_project_group#4/outputs/karn_dist_histogram.jpeg")
p <- ggplot(gender_clinical, aes(x=gender, y=age_at_initial_pathologic_diagnosis)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=2)
p + xlab("Sex") + ylab("Age at Initial Pathologic Diagnosis (years)") 
while (!is.null(dev.list()))  dev.off()
jpeg(filename = "../final_project_group#4/outputs/karn_age.jpeg")
p <- ggplot(clinical, aes(x=icd_10, y=age_at_initial_pathologic_diagnosis)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=2)
p + xlab("Malignant Neoplasm Site (ICD-10)") + ylab("Age at Initial Pathologic Diagnosis (years)")  + scale_x_discrete(breaks=c("C71.1","C71.2","C71.4", "C71.9"),
                                                                                                                       labels=c("Cerebrum", "Frontal Lobe", "Occipital Lobe", "Unspecified"))
while (!is.null(dev.list()))  dev.off()
clinical$karnofsky_class <- iflelse(clinical$karnofsky_performance_score < 50, "Poor", ifelse(clinical$karnofsky_performance_score < 80, "Moderate", "Excellent"))
clinical$survival_time <- ifelse(is.na(clinical$days_to_death),clinical$days_to_last_followup,clinical$days_to_death)
unique(clinical$survival_time) #contains some NA values 
surv_na_mask <- ifelse(is.na(clinical$survival_time), F, T)
na_text_mask <- ifelse(clinical$survival_time == "NA", F, T)
surv_inf_mask <- ifelse(is.infinite(clinical$survival_time), F, T)
surv_inf_text_mask <- ifelse(clinical$survival_time == "-Inf", F, T)
clinical <- clinical[surv_na_mask,]
clinical <- clinical[na_text_mask,]
clinical <- clinical[surv_inf_mask,]
clinical <- clinical[surv_inf_text_mask,]
clinical$death_event <- ifelse(is.na(clinical$days_to_death), F, T)
na_text_mask2 <- ifelse(is.na(clinical$death_event), F, T)
clinical <- clinical[na_text_mask2,]
# Initialize a 'survival' object, which contains the data we need.
surv_object_sex <- Surv(time = clinical$survival_time,
                        event = clinical$death_event)
# Create a fit object
# When writing formulas (x ~ y), x is what's being plotted and y is what is grouping x into categories
sex_fit <- surv_fit(surv_object_age ~ clinical$gender,
                    data = clinical)

# Feel free to play around with the margins and legend placement
survplot_sex = ggsurvplot(sex_fit,
                          pval=TRUE,
                          ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")),
                          xlab = "Time (days)",
                          legend = "right")

# when you create plots on your own, be sure to name them descriptively
KM_plot_sex = survplot_sex$plot +
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=10),
        legend.text = element_text(size=8))

jpeg(filename = "../final_project_group#4/outputs/siteKM.jpeg")
KM_plot_sex
while (!is.null(dev.list()))  dev.off()
# Initialize a 'survival' object, which contains the data we need.
surv_object_site <- Surv(time = clinical$survival_time,
                         event = clinical$death_event)
# Create a fit object
# When writing formulas (x ~ y), x is what's being plotted and y is what is grouping x into categories
site_fit <- surv_fit(surv_object_age ~ clinical$icd_10,
                     data = clinical)

# Feel free to play around with the margins and legend placement
survplot_site = ggsurvplot(site_fit,
                           pval=TRUE,
                           ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")),
                           xlab = "Time (days)",
                           legend = "right")

# when you create plots on your own, be sure to name them descriptively
KM_plot_site = survplot_site$plot +
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=10),
        legend.text = element_text(size=8))

jpeg(filename = "../final_project_group#4/outputs/siteKM.jpeg")
KM_plot_site
while (!is.null(dev.list()))  dev.off()
na_race_mask <- ifelse(clinical$race_list != "ASIAN" && clinical$race_list != "WHITE" && clinical$race_list != "BLACK OR AFRICAN AMERICAN", F, T)
race_clinical <- clinical[na_race_mask, ]

# Initialize a 'survival' object, which contains the data we need.
surv_object_race <- Surv(time = clinical$survival_time,
                         event = clinical$death_event)
# Create a fit object
# When writing formulas (x ~ y), x is what's being plotted and y is what is grouping x into categories
race_fit <- surv_fit(surv_object_race ~ clinical$race_list,
                     data = clinical)

# Feel free to play around with the margins and legend placement
survplot_race = ggsurvplot(race_fit,
                           pval=TRUE,
                           ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")),
                           xlab = "Time (days)",
                           legend = "right")

# when you create plots on your own, be sure to name them descriptively
KM_plot_race = survplot_race$plot +
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=10),
        legend.text = element_text(size=8))

jpeg(filename = "../final_project_group#4/outputs/raceKM.jpeg")
KM_plot_race
while (!is.null(dev.list()))  dev.off()
# Initialize a 'survival' object, which contains the data we need.
surv_object_age <- Surv(time = clinical$survival_time,
                        event = clinical$death_event)
# Create a fit object
# When writing formulas (x ~ y), x is what's being plotted and y is what is grouping x into categories
age_fit <- surv_fit(surv_object_age ~ clinical$age_class,
                    data = clinical )

# Feel free to play around with the margins and legend placement
survplot_age = ggsurvplot(age_fit,
                          pval=TRUE,
                          ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")),
                          xlab = "Time (days)",
                          legend = "right")

# when you create plots on your own, be sure to name them descriptively
KM_plot_age = survplot_age$plot +
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=10),
        legend.text = element_text(size=8))


jpeg(filename = "../final_project_group#4/outputs/ageKM.jpeg")
KM_plot_age
while (!is.null(dev.list()))  dev.off()
jpeg(filename ="../final_project_group#4/outputs/oncoplot.jpeg")
oncoplot(maf = maf_object,
         top= 10,
         clinicalFeatures = "age_class"
) 
while (!is.null(dev.list()))  dev.off()
#subsetting our MAF object into two MAF objects composed of old and young patients
young_mask = ifelse(maf_object@clinical.data$age_class == "Young", TRUE, FALSE)
young_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[young_mask]

old_mask = ifelse(maf_object@clinical.data$age_class == "Old", TRUE, FALSE)
old_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[old_mask]

young_maf <- subsetMaf(maf = maf_object,
                       tsb = young_patient_barcodes)

old_maf <- subsetMaf(maf = maf_object,
                     tsb = old_patient_barcodes)
jpeg(filename = "../final_project_group#4/outputs/co_oncoplot.jpeg")
coOncoplot(m1 = young_maf,
           m2 = old_maf,
           m1Name = "Young",
           m2Name = "Old"
)
while (!is.null(dev.list()))  dev.off()
jpeg(filename = "../final_project_group#4/outputs/atrx_lolli.jpeg")
lollipopPlot2(m1 = young_maf,
              m2 = old_maf,
              m1_name = "Young",
              m2_name = "Old",
              gene = "ATRX") ## pick any gene of your choosing to fill in here
while (!is.null(dev.list()))  dev.off()
jpeg(filename = "../final_project_group#4/outputs/1dh1_lolli.jpeg")
lollipopPlot2(m1 = young_maf,
              m2 = old_maf,
              m1_name = "Young",
              m2_name = "Old",
              gene = "IDH1") ## pick any gene of your choosing to fill in here
while (!is.null(dev.list()))  dev.off()
jpeg(filename = "../final_project_group#4/outputs/tp53_lolli.jpeg")
lollipopPlot2(m1 = young_maf,
              m2 = old_maf,
              m1_name = "Young",
              m2_name = "Old",
              gene = "TP53") ## pick any gene of your choosing to fill in here
while (!is.null(dev.list()))  dev.off()
jpeg(filename = "../final_project_group#4/outputs/maf_somatic.jpeg")
somaticInteractions(maf = maf_object,
                    top = 25,
                    pvalue = c(0.05, 0.1))
while (!is.null(dev.list()))  dev.off()
jpeg(filename = "../final_project_group#4/outputs/young_somatic.jpeg")
somaticInteractions(maf = young_maf,
                    top = 25,
                    pvalue = c(0.05, 0.1))
while (!is.null(dev.list()))  dev.off()
jpeg(filename = "../final_project_group#4/outputs/old_somatic.jpeg")
somaticInteractions(maf = old_maf,
                    top = 25,
                    pvalue = c(0.05, 0.1))
while (!is.null(dev.list()))  dev.off()
unique(clinical.drug$drug_name) #there are 54 drugs, we will only look at the top 4
clinical.drug$drug_name[clinical.drug$drug_name == "Temodar"] <- "Temozolomide"
clinical.drug$drug_name[clinical.drug$drug_name == "Temodor"] <- "Temozolomide"
clinical.drug$drug_name[clinical.drug$drug_name == "Temador"] <- "Temozolomide"
clinical.drug$drug_name[clinical.drug$drug_name == "Temozolamide"] <- "Temozolomide"
clinical.drug$drug_name[clinical.drug$drug_name == "CCNU"] <- "Lomustine"
clinical.drug$drug_name[clinical.drug$drug_name == "BCNU"] <- "Carmustine"
clinical.drug$drug_name[clinical.drug$drug_name == "Carmustine (BCNU)"] <- "Carmustine"
clinical.drug$drug_name[clinical.drug$drug_name == "Avastin"] <- "Bevacizumab"
counts <- count(clinical.drug,vars=  drug_name) #count the number of times each drug is used
counts
ordered_counts <- counts[order(counts$n),]
ordered_counts #order counts to obtain top 3 drugs administered, GO TO BACK OF TABLE
top_drugs = c("Temozolomide", "Bevacizumab", "Lomustine", "Carmustine")
top_frequency = ordered_counts[ordered_counts$vars %in% top_drugs,]
top_frequency
#merge clinical and clinical.drug df's using patient ID as common column; needed to visualize using ggplot
total_clinical <- merge (clinical, clinical.drug, by.x = 'Tumor_Sample_Barcode', by.y ='bcr_patient_barcode')
modified_clinical <- droplevels(total_clinical[total_clinical$drug_name %in% top_drugs,]) #filters df to include only T3 drugs
modified_clinical$drug_name <- tolower(modified_clinical$drug_name)
modified_clinical <- modified_clinical[ifelse(is.na(modified_clinical$drug_name) == TRUE | modified_clinical$drug_name == "", F, T),]
modified_clinical$survival_time <- ifelse(is.na(modified_clinical$days_to_death),modified_clinical$days_to_last_followup,modified_clinical$days_to_death)
unique(modified_clinical$survival_time) #contains some NA values 
surv_na_mask <- ifelse(is.na(modified_clinical$survival_time), F, T)
na_text_mask <- ifelse(modified_clinical$survival_time == "NA", F, T)
surv_inf_mask <- ifelse(is.infinite(modified_clinical$survival_time), F, T)
modified_clinical <- modified_clinical[surv_na_mask,] #removes NA values
modified_clinical <- modified_clinical[na_text_mask,] #removes "NA" text values
modified_clinical <- modified_clinical[surv_inf_mask,] #removes infinity values
modified_clinical$death_event <- ifelse(is.na(modified_clinical$days_to_death), F, T)
unique(modified_clinical$death_event)
na_text_mask2 <- ifelse(is.na(modified_clinical$death_event), F, T)
modified_clinical <- modified_clinical[na_text_mask2,]
surv_object_age <- Surv(time = modified_clinical$survival_time,
                        event = modified_clinical$death_event)
drug_fit <- surv_fit(surv_object_age ~ modified_clinical$drug_name,
                     data = modified_clinical)
survplot_drug = ggsurvplot(drug_fit,
                           pval=TRUE, 
                           ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                           legend = "right")
KM_plot_drug = survplot_drug$plot + 
  theme_bw() +  
  theme(axis.title = element_text(size=20), 
        axis.text = element_text(size=16),
        legend.title = element_text(size=8),
        legend.text = element_text(size=6))
jpeg(filename = "../final_project_group#4/outputs/overall_drug_KM.jpeg")
KM_plot_drug
while (!is.null(dev.list()))  dev.off()
#merge clinical and clinical.drug df's using patient ID as common column; needed to visualize using ggplot
young_total_clinical <- merge (clinical.drug, young_clinical, by.x = 'bcr_patient_barcode', by.y ='bcr_patient_barcode')
young_modified_clinical <- droplevels(young_total_clinical[young_total_clinical$drug_name %in% top_drugs,]) #filters df to include only T3 drugs
young_modified_clinical$drug_name <- tolower(young_modified_clinical$drug_name)
young_modified_clinical$survival_time <- ifelse(is.na(young_modified_clinical$days_to_death), young_modified_clinical$days_to_last_followup,young_modified_clinical$days_to_death)
surv_na_mask <- ifelse(is.na(young_modified_clinical$survival_time), F, T)
na_text_mask <- ifelse(young_modified_clinical$survival_time == "NA", F, T)
surv_inf_mask <- ifelse(is.infinite(young_modified_clinical$survival_time), F, T)
young_modified_clinical <- young_modified_clinical[surv_na_mask,] #removes NA values
young_modified_clinical <- young_modified_clinical[na_text_mask,] #removes "NA" text values
young_modified_clinical <- young_modified_clinical[surv_inf_mask,] #removes infinity values
young_modified_clinical$death_event <- ifelse(is.na(young_modified_clinical$days_to_death), F, T)
na_text_mask2 <- ifelse(is.na(young_modified_clinical$death_event), F, T)
young_modified_clinical <- young_modified_clinical[na_text_mask2,]
surv_object_age <- Surv(time = young_modified_clinical$survival_time,
                        event = young_modified_clinical$death_event)
drug_fit <- surv_fit(surv_object_age ~ young_modified_clinical$drug_name,
                     data = young_modified_clinical)
survplot_drug = ggsurvplot(drug_fit,
                           pval=TRUE, 
                           ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                           legend = "right")
KM_plot_drug_young = survplot_drug$plot + 
  theme_bw() +  
  theme(axis.title = element_text(size=20), 
        axis.text = element_text(size=16),
        legend.title = element_text(size=8),
        legend.text = element_text(size=6))
jpeg(filename = "../final_project_group#4/outputs/young_drug_KM.jpeg")
KM_plot_drug_young
while (!is.null(dev.list()))  dev.off()
old_total_clinical <- merge (old_clinical, clinical.drug, by.x = 'bcr_patient_barcode', by.y ='bcr_patient_barcode')
old_modified_clinical <- droplevels(old_total_clinical[old_total_clinical$drug_name %in% top_drugs,]) #filters df to include only T3 drugs
old_modified_clinical$drug_name <- tolower(old_modified_clinical$drug_name)
old_modified_clinical$survival_time <- ifelse(is.na(old_modified_clinical$days_to_death),old_modified_clinical$days_to_last_followup,old_modified_clinical$days_to_death)
surv_na_mask <- ifelse(is.na(old_modified_clinical$survival_time), F, T)
na_text_mask <- ifelse(old_modified_clinical$survival_time == "NA", F, T)
surv_inf_mask <- ifelse(is.infinite(old_modified_clinical$survival_time), F, T)
old_modified_clinical <- old_modified_clinical[surv_na_mask,] #removes NA values
old_modified_clinical <- old_modified_clinical[na_text_mask,] #removes "NA" text values
old_modified_clinical <- old_modified_clinical[surv_inf_mask,] #removes infinity values
old_modified_clinical$death_event <- ifelse(is.na(old_modified_clinical$days_to_death), F, T)
na_text_mask2 <- ifelse(is.na(old_modified_clinical$death_event), F, T)
old_modified_clinical <- old_modified_clinical[na_text_mask2,]
surv_object_age <- Surv(time = old_modified_clinical$survival_time,
                        event = old_modified_clinical$death_event)
drug_fit <- surv_fit(surv_object_age ~ old_modified_clinical$drug_name,
                     data = old_modified_clinical)
survplot_drug = ggsurvplot(drug_fit,
                           pval=TRUE, 
                           ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                           legend = "right")
KM_plot_drug_old = survplot_drug$plot + 
  theme_bw() +  
  theme(axis.title = element_text(size=20), 
        axis.text = element_text(size=16),
        legend.title = element_text(size=8),
        legend.text = element_text(size=6))
jpeg(filename = "../final_project_group#4/outputs/old_drug_KM.jpeg")
KM_plot_drug_old
while (!is.null(dev.list()))  dev.off()
write.csv(clinical, "../final_project_group#4/outputs/gbm_clinical.csv")
write.csv(maf, "../final_project_group#4/outputs/gbm_genes.csv")
write.csv(rna_query, "../final_project_group#4/outputs/gbm_rna_query.csv")
