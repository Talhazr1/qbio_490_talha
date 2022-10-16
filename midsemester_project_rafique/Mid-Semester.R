"---
project: Mid-Semester Project
author: Talha Rafique 
---"
#----------------------------------------------
#Output Folder Setup

#creates outputs folder
dir.create("/Users/talhaaccount/Desktop/Qbio_490_R/qbio_490_talha/midsemester_project_rafique/outputs")
#sets working directory
setwd("/Users/talhaaccount/Desktop/Qbio_490_R/qbio_490_talha/midsemester_project_rafique/outputs")

#----------------------------------------------
#Loading Packages

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15")
library(BiocManager)

if (!require("TCGAbiolinks", quietly = TRUE)) 
  BiocManager::install("TCGAbiolinks") 

library(TCGAbiolinks)
library(ggplot2)
if (!require(survival)){ 
  install.packages("survival")
}
library(survival)
if (!require(survminer)){ 
  install.packages("survminer")
}
library(survminer)
library(maftools)

#----------------------------------------------
#Query all my Data 

#Clinical stuff
clinical_query <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical", file.type = "xml")
GDCdownload(clinical_query) 
clinical <- GDCprepare_clinic(clinical_query, clinical.info = "patient") #placed into data frame
clinical.drug <- GDCprepare_clinic(query = clinical_query, clinical.info = "drug")



# Changed the name of "bcr_patient_barcode" to "Tumor_Sample_Barcode" since present in MAF too
colnames(clinical)[ colnames(clinical) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

#MAF stuff
maf_query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", 
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(maf_query)
maf <- GDCprepare(maf_query)
maf_object <- read.maf(maf = maf, 
                       clinicalData = clinical,  
                       isTCGA = TRUE)

#----------------------------------------------
#My first analysis of the Clinical data will look at  histological type (0 NA's present), 
#The goal of this analysis is to create a Kaplan-Meier plot to display the survival times
#depending on which histological type is present 

#Since the histological type have no NA's present, and are already stratified we are good on that

#To create a KM plot we need the survival time of each patient
clinical$survival_time <- ifelse(is.na(clinical$days_to_death), clinical$days_to_last_followup, clinical$days_to_death)
#We also need to know when each patient passed away
clinical$death_occurance <- ifelse(!is.na(clinical$vital_status),TRUE,FALSE)

#Using the death event and surival time we can now use our grouping variable, which are the histological
#types to analyze how the histological type affected patients lifespan

# Initialize a 'survival' object, which contains the data we need.
surv_object_age <- Surv(time = clinical$survival_time,
                        event = clinical$death_occurance)

# Create a fit object
# When writing formulas (x ~ y), x is what's being plotted and y is what is grouping x into categories
age_fit <- surv_fit( surv_object_age ~ clinical$histological_type,
                     data = clinical )

# the ggtheme and legend arguments are for formatting. 
# Feel free to play around with the margins and legend placement
survplot_age = ggsurvplot(age_fit, 
                          pval=TRUE, 
                          ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                          legend = "right")

# when you create plots on your own, be sure to name them descriptively
KM_plot_histological = survplot_age$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))

pdf("/Users/talhaaccount/Desktop/Qbio_490_R/qbio_490_talha/midsemester_project_rafique/outputs/KM_plot_histological.pdf")
KM_plot_histological
dev.off()
#----------------------------------------------
#This second section is a box plot graph based on the average survival time depending on 
#the type of histological cancer type

box_plot_histological = plot(clinical$histological_type, clinical$survival_time)
box_plot_histological

#relative path does NOT work so used absolute path
pdf("/Users/talhaaccount/Desktop/Qbio_490_R/qbio_490_talha/midsemester_project_rafique/outputs/box_plot_histological.pdf")
box_plot_histological
dev.off()

#----------------------------------------------
#Here I used a oncoplot as my MAF plot 
onc = oncoplot(maf = maf_object,
         top = 10,
         clinicalFeatures = "days_to_last_known_alive")

#relative path does NOT work so used absolute path
pdf("/Users/talhaaccount/Desktop/Qbio_490_R/qbio_490_talha/midsemester_project_rafique/outputs/onc.pdf")
onc
dev.off()


#----------------------------------------------
#usage of Clinical drugs

box_plot_drug_name = plot(clinical$histological_type, clinical$survival_time)
box_plot_drug_name
#relative path does NOT work so used absolute path

pdf("/Users/talhaaccount/Desktop/Qbio_490_R/qbio_490_talha/midsemester_project_rafique/outputs/box_plot_drug_name.pdf")
box_plot_drug_name
dev.off()

#----------------------------------------------




