################################################################################
#################### Setup #####################################################

# loc
setwd("~/git/CAIRO5_Public")

# libraries
library(readxl)
library(data.table)
library(ggplot2)
library(scales)
library(gtools)

# Data
BloodSampleCharacteristics <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=2, skip=1))
ClinicalDemographics <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=1, skip=1))

################################################################################
################# pair down data to required cols ##############################

# Blood Samples
keep <- c('Sequence ID', 'Subject ID', 'Sample Time Point', 'DELFI-TF')
BloodSampleCharacteristics <- BloodSampleCharacteristics[,..keep]

# clinical data
keep <- c('Subject ID', 'Study Arm', 'Metastasis Type')
ClinicalDemographics <- ClinicalDemographics[,..keep]

# merge data together
finalData <- merge(BloodSampleCharacteristics, ClinicalDemographics, by='Subject ID', all.x=T)

################################################################################
############## format data #####################################################

# restrict to Baseline and make sure we have a DELFI-TF value
finalData <- finalData[`Sample Time Point` == 'T0']
finalData <- finalData[!(is.na(`DELFI-TF`))]

################################################################################
############## run stats #######################################################

# switch statement to re-define p-values
relabel_p <- function(x){
  if(x < .0001) return('****')
  if(x < .001) return('***')
  if(x < .01) return('**')
  if(x < .05) return('*')
  return('ns')
}

# wilcox test
delfiLabel <- wilcox.test(finalData[`Metastasis Type` == 'Synchronous']$`DELFI-TF`, finalData[`Metastasis Type` == 'Metachronous']$`DELFI-TF`)
delfiLabel <- relabel_p(delfiLabel$p.value)

# labelDT
statDT <- data.table(variable=c("DELFI-TF"),
                     y=c(.65),
                     label=c(delfiLabel))

################################################################################
############## Make Plot #######################################################

# melt data for plot
finalData <- melt(finalData, measure.vars=c('DELFI-TF'))

# make the plot
Fig3C <- ggplot() +
  geom_boxplot(data=finalData, mapping=aes(x=variable, y=value, fill=`Metastasis Type`)) +
  geom_jitter(data=finalData, mapping=aes(x=variable, y=value, fill=`Metastasis Type`), position = position_jitterdodge(jitter.width=.2), size=4, color='grey30', show.legend = FALSE) +
  geom_text(data=statDT, mapping=aes(x=variable, y=y, label=label), size=20, color='grey30') +
  annotate("segment", x=.75, xend=1.25, y=.62, yend=.62, linewidth=1.75, linetype='solid', color='gray65', lineend = "round") +
  ylab("Tumor Fraction") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual("Metastasis\nTiming", values=c("Metachronous"="#899DA4", "Synchronous"="#74A089")) +
  guides(fill=guide_legend(ncol=1)) +
  theme_bw() +
  theme(legend.position = "top") +
  theme(legend.title=element_text(size=40), legend.text=element_text(size=36)) +
  theme(legend.key.size = unit(2, 'cm')) +
  theme(axis.text.x=element_text(size=36), axis.text.y=element_text(size=36)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_text(size=40)) +
  theme(strip.background=element_rect(fill="black")) +
  theme(strip.text=element_text(color="white", size=36, face="bold")) +
  theme(panel.grid=element_blank())

# save plot
png(filename="Figures/Fig3C.png", height=16, width=13, units="in", res=150)
Fig3C
dev.off()

cairo_pdf(file="Figures/Fig3C.pdf", height=16, width=13)
Fig3C
dev.off()



