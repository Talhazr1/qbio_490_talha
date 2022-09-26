"---
title: Group Project
author: Talha Rafique 
---"
#ALL THE SETUP 
setwd("/Users/talhaaccount/Desktop/Qbio_490_R/qbio_490_talha/analysis_data/")
clinical <- read.csv("/Users/talhaaccount/Desktop/Qbio_490_R/qbio_490_talha/analysis_data/brca_clinical_data.csv")
library(BiocManager)
library(TCGAbiolinks)
clinical_query <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical", file.type = "xml")
clinical.drug <- GDCprepare_clinic(query = clinical_query, clinical.info = "drug")
clinical.rad <- GDCprepare_clinic(query = clinical_query, clinical.info = "radiation")
#-------------------------------------
#check for NA's in breast cancer 
naInLymph <- is.na(clinical$lymph_node_examined_count)
sum(naInLymph) # = 139 NA's

#checking for NA's in therapy ongoing
naInStage <- is.na(clinical$stage_event_pathologic_stage)
sum(naInStage) # = 0 NA's
#-------------------------------------

#1 plotting both variables against each other
lymphNode <- (clinical$lymph_node_examined_count) #place all lymph node counts into variable 
stage <- (clinical$stage_event_pathologic_stage)  #place all stage tyypes into variable 
boxplot(formula = lymphNode ~ stage,              #ploting the lymphnode count against stage of cancer
        data = clinical,
        xlab = "Stage",                #x-axis is the categorical variable of what stage the breast cancer is at
        ylab = "Lymph Node Count")     #y-axis is the discrete variables of how many lymphnodes a person has

#---------------------------------------------------

# Importing necessary survival packages
library(ggplot2)
if (!require(survival)){ 
  install.packages("survival")
}
library(survival)
if (!require(survminer)){ 
  install.packages("survminer")
}
library(survminer)

#---------------------------------------------------

#stratification for lymph nodes, so we get a more readable KM Plot 
# If lymph node count < 10 then "Minimal"
# If lymph node count < 20 then "Low"
# If lymph node count < 30 then "Medium"
# If lymph node count < 40 then "Elevated"
# If lymph node count > 40 then "High"
clinical$lymphConcentration <- ifelse(lymphNode <= 10, "Minimal", ifelse(lymphNode <= 20, "Low",ifelse(lymphNode <= 30, "Medium", ifelse(lymphNode <= 40, "Elavated", "High"))))

#---------------------------------------------------

#2 - survival analysis on first variable, stage of cancer
clinical$death_event <- ifelse(!is.na(clinical$vital_status),TRUE,FALSE)     #death event
clinical$survival_time <- ifelse(is.na(clinical$days_to_death), clinical$days_to_last_followup, clinical$days_to_death)  #total survival time

?Surv
surv_object_age <- Surv(time = clinical$survival_time,
                        event = clinical$death_event)

# Create a fit object
# When writing formulas (x ~ y), x is what's being plotted and y is what is grouping x into categories
?surv_fit
age_fit <- surv_fit( surv_object_age ~ stage,          #comparing survival to stage of breast cancer
                     data = clinical )

# the ggtheme and legend arguments are for formatting. 
# Feel free to play around with the margins and legend placement
?ggsurvplot
survplot_age = ggsurvplot(age_fit, 
                          pval=TRUE, 
                          ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                          legend = "right")

# when you create plots on your own, be sure to name them descriptively
KM_plot_stage = survplot_age$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))

KM_plot_stage
#---------------------------------------------------

#3 - survival analysis on second variable, lymphnode count
?Surv
surv_object_age <- Surv(time = clinical$survival_time,
                        event = clinical$death_event)

# Create a fit object
# When writing formulas (x ~ y), x is what's being plotted and y is what is grouping x into categories
?surv_fit
age_fit <- surv_fit( surv_object_age ~ clinical$lymphConcentration,    #comparing survival to lympnode concentration
                     data = clinical )

# the ggtheme and legend arguments are for formatting. 
# Feel free to play around with the margins and legend placement
?ggsurvplot
survplot_age = ggsurvplot(age_fit, 
                          pval=TRUE, 
                          ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                          legend = "right")

# when you create plots on your own, be sure to name them descriptively
KM_plot_lymph = survplot_age$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))

KM_plot_lymph
#---------------------------------------------------

#Exporting plots and saving data frames
write.csv(clinical, "/Users/talhaaccount/Desktop/Qbio_490_R/qbio_490_talha/week5/brca_clinical_data.csv", row.names = FALSE)

#PDF saving of lymph Node count compared against stage of breast Cancer
pdf("/Users/talhaaccount/Desktop/Qbio_490_R/qbio_490_talha/week5/lymphToStageCompare.pdf")
boxplot(formula = lymphNode ~ stage,
        data = clinical,
        xlab = "anatomic neoplasm subdivisions",
        ylab = "Age of Patient at diagnosis")
dev.off()

#PDF saving of Breast Cancer stage survival Analysis
pdf("/Users/talhaaccount/Desktop/Qbio_490_R/qbio_490_talha/week5/stageSurvival.pdf")
KM_plot_stage  
dev.off()

#PDF saving of Lymph Node count survival Analysis
pdf("/Users/talhaaccount/Desktop/Qbio_490_R/qbio_490_talha/week5/lymphSurvival.pdf")
KM_plot_lymph
dev.off()
#---------------------------------------------------










