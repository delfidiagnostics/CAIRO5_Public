################################################################################
#################### Setup #####################################################

# loc
setwd("~/git/CAIRO5_Public")

# libraries
library(readxl)
library(data.table)
library(gt)
library(survival)
library(survminer)

# Data
BloodSampleCharacteristics <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=2, skip=1))
ClinicalDemographics <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=1, skip=1))
SumOfLongestDiameters <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=7, skip=1))

################################################################################
################# pair down data to required cols ##############################

# Blood Samples
keep <- c('Sequence ID', 'Subject ID', 'Sample Time Point', 'DELFI-TF', 'CEA')
BloodSampleCharacteristics <- BloodSampleCharacteristics[,..keep]

# clinical data
keep <- c('Subject ID', 'OS Registration', 'OS Status (0=alive; 1=death)', 'Primary Tumor Sidedness', 'Age Bin')
ClinicalDemographics <- ClinicalDemographics[,..keep]

# Sum of Longest Diameter, pair down to be within 60 days
SumOfLongestDiameters <- SumOfLongestDiameters[abs(`Time Between SLD and Liquid Biopsy (Days)`) <= 60]
SumOfLongestDiameters[,`Time Between SLD and Liquid Biopsy (Days)` := NULL]

# merge data together
finalData <- merge(BloodSampleCharacteristics, ClinicalDemographics, by='Subject ID', all.x=T)
finalData <- merge(finalData, SumOfLongestDiameters, by=c('Subject ID', 'Sample Time Point'), all.x=T)

################################################################################
############### format data for analysis #######################################

# restrict to baseline
finalData <- finalData[`Sample Time Point` == 'T0']

# factor cols so we know the order
finalData[,`Age Bin` := factor(`Age Bin`, levels=c("< 65", "≥ 65"))]
finalData[,`Primary Tumor Sidedness` := factor(`Primary Tumor Sidedness`, levels=c("Left", "Right"))]

################################################################################
########### perform survival analysis ##########################################

# calculate hazard ratio
finalData <- na.omit(finalData)
coxMultivariate <- summary(coxph(Surv(`OS Registration`, `OS Status (0=alive; 1=death)`) ~ `Age Bin` + `Sum of the Longest Diameters` + CEA + `DELFI-TF` + `Primary Tumor Sidedness`, data = finalData))
coxMultivariate

# make a DT
coxphDT <- data.table('HR'=c(coxMultivariate$conf.int[1,1],
                             coxMultivariate$conf.int[2,1],
                             coxMultivariate$conf.int[3,1],
                             coxMultivariate$conf.int[4,1],
                             coxMultivariate$conf.int[5,1]),
                      'low'=c(coxMultivariate$conf.int[1,3],
                              coxMultivariate$conf.int[2,3],
                              coxMultivariate$conf.int[3,3],
                              coxMultivariate$conf.int[4,3],
                              coxMultivariate$conf.int[5,3]),
                      'high'=c(coxMultivariate$conf.int[1,4],
                               coxMultivariate$conf.int[2,4],
                               coxMultivariate$conf.int[3,4],
                               coxMultivariate$conf.int[4,4],
                               coxMultivariate$conf.int[5,4]),
                      'p'=c(signif(coxMultivariate$coefficients[1,5], 3),
                                  signif(coxMultivariate$coefficients[2,5], 3),
                                  signif(coxMultivariate$coefficients[3,5], 3),
                                  signif(coxMultivariate$coefficients[4,5], 3),
                                  signif(coxMultivariate$coefficients[5,5], 3)),
                      'Description'=c('< 65 v ≥ 65',
                                      'SLD',
                                      'CEA Level',
                                      'DELFI-TF',
                                      'Left v Right'))

# reformat table to 2/3 sigfigs
coxphDT[,HR := signif(HR, 3)]
coxphDT[,low := signif(low, 3)]
coxphDT[,high := signif(high, 3)]

# sort table by HR and re-arrange cols
coxphDT <- coxphDT[order(-HR)]
coxphDT <- coxphDT[,.(Description, HR, low, high, `p`)]

################################################################################
######### Make table ###########################################################

# make the table
coxph_gt <- gt(coxphDT) %>%
  tab_header(title = "Overall Survival HR", subtitle = paste0("Multivariate, (n=", nrow(finalData), ")")) %>%
  tab_options(heading.background.color = "#1B465A", heading.title.font.size = 25, heading.subtitle.font.size = 20) %>%
  tab_spanner(label="Confidence Interval 95%", columns=c("low", "high")) %>%
  tab_footnote(footnote = paste0("Log-Rank Test: ", ifelse(coxMultivariate$sctest[3] < .0001, 'p<0.0001'))) %>%
  cols_align(align = c("center"), columns = c("HR", "low", "high", "p"))

# save the table
gtsave(coxph_gt, filename = "Fig3E.png", path="Figures/")
gtsave(coxph_gt, filename = "Fig3E.pdf", path="Figures/")

