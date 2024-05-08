################################################################################
#################### Setup #####################################################

# loc
setwd("~/git/CAIRO5_Public")

# libraries
library(readxl)
library(data.table)
library(ggplot2)
library(scales)

# Data
BloodSampleCharacteristics <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=2, skip=1))
NonCancerControl <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=4, skip=1))
ClinicalDemographics <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=1, skip=1))
recistCategory <- readRDS('Data/Fig2b_AuxData.rds')

################################################################################
################# Combine Data #################################################

# Blood Samples
keep <- c('Sequence ID', 'Subject ID', 'Sample Time Point', 'DELFI-TF')
BloodSampleCharacteristics <- BloodSampleCharacteristics[,..keep]
BloodSampleCharacteristics[,SampleType := 'Cancer']

# Non Cancers
keep <- c('Sequence ID', 'DELFI-TF')
NonCancerControl <- NonCancerControl[,..keep]
NonCancerControl[,SampleType := 'Control']

# clinical data
keep <- c('Subject ID', 'Study Arm')
ClinicalDemographics <- ClinicalDemographics[,..keep]

# merge data together
finalData <- merge(BloodSampleCharacteristics, ClinicalDemographics, by='Subject ID', all.x=T)
finalData <- rbindlist(list(finalData, NonCancerControl), fill=TRUE)

################################################################################
############### final data formatting ##########################################

# restrict to relevant time points
finalData <- finalData[SampleType == 'Control' | (SampleType == 'Cancer' & `Sample Time Point` %chin% c('T0', 'T1'))]


# add labels with counts
SampleCounts <- finalData[,.N,by=.(`Sample Time Point`, SampleType)]
setkeyv(SampleCounts, c("Sample Time Point", "SampleType"))
setkeyv(finalData, c("Sample Time Point", "SampleType"))
finalData <- merge(finalData, SampleCounts)
finalData[SampleType == 'Cancer' & `Sample Time Point` == 'T0', Label := paste0("Colorectal Cancers\nPre-Treatment (T0)\n(n=", N, ")")]
finalData[SampleType == 'Cancer' & `Sample Time Point` == 'T1', Label := paste0("Colorectal Cancers\nPost-Treatment (T1)\n(n=", N, ")")]
finalData[SampleType == 'Control', Label := paste0("Non-Cancer\nReference\n(n=", N, ")")]

# Annotate Mutation Status from Arm, ARM1 is MT
finalData[,MutationStatus := ifelse(`Study Arm` == 1, 'MT', 'WT')]

################################################################################
############## Also create boxplot data for RECIST Categories ##################

# merge RECIST categories within 60 day windows
BloodSampleCharacteristics <- merge(BloodSampleCharacteristics, ClinicalDemographics, by='Subject ID', all.x=T)
recistFinalData <- merge(BloodSampleCharacteristics, recistCategory, by.x='Sequence ID', by.y='delfiID')

# restrict down to relevant rows and columns
recistFinalData <- recistFinalData[TimePoint == 'Progression',.(`Sequence ID`, `DELFI-TF`, `Study Arm`)]

# Annotate Mutation Status from Arm, ARM1 is MT
recistFinalData[,MutationStatus := ifelse(`Study Arm` == 1, 'MT', 'WT')]

# make label
recistFinalData[, Label := paste0("RECIST\nProgression\n(n=", nrow(recistFinalData), ")")]

################################################################################
############## plot formatting #################################################

# bind RECIST boxplot data on
finalData <- rbindlist(list(finalData, recistFinalData), use.names=T, fill=T)

# set factor order (note if anything changes with the counts change this)
finalData[,Label := factor(Label, levels=c("Non-Cancer\nReference\n(n=153)",
                                           "Colorectal Cancers\nPre-Treatment (T0)\n(n=128)",
                                           "Colorectal Cancers\nPost-Treatment (T1)\n(n=151)",
                                           "RECIST\nProgression\n(n=73)"))]

################################################################################
############### Run Stats ######################################################

# run stats
statsT0_T1 <- wilcox.test(finalData[`Sample Time Point` == 'T1' & SampleType == 'Cancer']$`DELFI-TF`, finalData[`Sample Time Point` == 'T0' & SampleType == 'Cancer']$`DELFI-TF`)

statsT0_Non <- wilcox.test(finalData[SampleType == 'Control']$`DELFI-TF`, finalData[`Sample Time Point` == 'T0' & SampleType == 'Cancer']$`DELFI-TF`)

statsT1_Non <- wilcox.test(finalData[SampleType == 'Control']$`DELFI-TF`, finalData[`Sample Time Point` == 'T1' & SampleType == 'Cancer']$`DELFI-TF`)

statsT0_MT <- wilcox.test(finalData[`Sample Time Point` == 'T0' & SampleType == 'Cancer' & MutationStatus == 'WT']$`DELFI-TF`, finalData[`Sample Time Point` == 'T0' & SampleType == 'Cancer' & MutationStatus == 'MT']$`DELFI-TF`)

statsT1_MT <- wilcox.test(finalData[`Sample Time Point` == 'T1' & SampleType == 'Cancer' & MutationStatus == 'WT']$`DELFI-TF`, finalData[`Sample Time Point` == 'T1' & SampleType == 'Cancer' & MutationStatus == 'MT']$`DELFI-TF`)

statsT1_Tprogression <- wilcox.test(finalData[`Sample Time Point` == 'T1' & SampleType == 'Cancer']$`DELFI-TF`, finalData[Label == 'RECIST\nProgression\n(n=73)']$`DELFI-TF`)

statsTprogression_MT <- wilcox.test(finalData[Label == 'RECIST\nProgression\n(n=73)' & MutationStatus == 'MT']$`DELFI-TF`, finalData[Label == 'RECIST\nProgression\n(n=73)' & MutationStatus == 'WT']$`DELFI-TF`)

# switch statement to re-define p-values
relabel_p <- function(x){
  if(x < .0001) return('****')
  if(x < .001) return('***')
  if(x < .01) return('**')
  if(x < .05) return('*')
  return('ns')
}

# contruct data.table to plot p-vals
#!!!!!!!! if factor order changes this needs to change !!!!!!!!
statsT0_T1 <- data.table(x=c('Colorectal Cancers\nPre-Treatment (T0)\n(n=128)'),
                         xend=c('Colorectal Cancers\nPost-Treatment (T1)\n(n=151)'),
                         y=c(.62),
                         label=relabel_p(statsT0_T1$p.value),
                         label_x=c(2.5),
                         label_y=c(.64))
statsT0_Non <- data.table(x=c('Colorectal Cancers\nPre-Treatment (T0)\n(n=128)'),
                         xend=c('Non-Cancer\nReference\n(n=153)'),
                         y=c(.76),
                         label=relabel_p(statsT0_Non$p.value),
                         label_x=c(1.5),
                         label_y=c(.78))
statsT1_Non <- data.table(x=c('Colorectal Cancers\nPost-Treatment (T1)\n(n=151)'),
                          xend=c('Non-Cancer\nReference\n(n=153)'),
                          y=c(.69),
                          label=relabel_p(statsT1_Non$p.value),
                          label_x=c(2),
                          label_y=c(.71))
statsT1_Tprogression <- data.table(x=c('Colorectal Cancers\nPost-Treatment (T1)\n(n=151)'),
                          xend=c('RECIST\nProgression\n(n=73)'),
                          y=c(.55),
                          label=relabel_p(statsT1_Tprogression$p.value),
                          label_x=c(3.5),
                          label_y=c(.57))
statData <- rbindlist(list(statsT0_T1, statsT0_Non, statsT1_Non, statsT1_Tprogression))

################################################################################
############## Make Plot #######################################################

# make the plot
Fig2B <- ggplot() +
  geom_boxplot(data=finalData, mapping=aes(x=Label, y=`DELFI-TF`, fill=MutationStatus), outlier.shape=NA) +
  geom_point(data=finalData, mapping=aes(x=Label, y=`DELFI-TF`, fill=MutationStatus), position=position_jitterdodge(jitter.width = .2), size=3, alpha=.75, color='grey30', show.legend = FALSE) +
  geom_segment(data=statData, mapping=aes(x=x, xend=xend, y=y, yend=y), linewidth=1.75, linetype='solid', color='gray65', lineend = "round") +
  geom_text(data=statData, mapping=aes(x=label_x, y=label_y, label=label), size=12, color='grey40') +
  ylab("DELFI-TF") +
  scale_fill_manual("RAS/BRAF Mutation", values=c("WT"="#E6A0C4", "MT"="#CB4335"), na.translate = F) +
  scale_y_continuous(labels=scales::percent) +
  theme_bw() +
  theme(legend.position = "top") +
  theme(legend.text=element_text(size=36), legend.title=element_text(size=42)) +
  theme(legend.key.size = unit(2, 'cm')) +
  theme(axis.text.x=element_text(size=34), axis.text.y=element_text(size=34)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_text(size=42)) +
  theme(panel.grid=element_blank())

# save plot
png(filename="Figures/Fig2B.png", height=14, width=20, units="in", res=150)
Fig2B
dev.off()

cairo_pdf(file="Figures/Fig2B.pdf", height=14, width=20)
Fig2B
dev.off()

