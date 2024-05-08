################################################################################
#################### Setup #####################################################

# loc
setwd("~/git/CAIRO5_Public")

# libraries
library(readxl)
library(data.table)
library(ggplot2)
library(gtools)
library(scales)

# Data
BloodSampleCharacteristics <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=2, skip=1))
ClinicalDemographics <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=1, skip=1))

################################################################################
################# pair down data to required cols ##############################

# Blood Samples
keep <- c('Sequence ID', 'Subject ID', 'Sample Time Point', 'DELFI-TF')
BloodSampleCharacteristics <- BloodSampleCharacteristics[,..keep]

# clinical data
keep <- c('Subject ID', 'Study Arm', 'PFS Status (0=alive; 1=progression or death)')
ClinicalDemographics <- ClinicalDemographics[,..keep]

# merge data together
finalData <- merge(BloodSampleCharacteristics, ClinicalDemographics, by='Subject ID', all.x=T)

################################################################################
################### final formatting ###########################################

# name progression status for plot
finalData[,progressionStatus := ifelse(`PFS Status (0=alive; 1=progression or death)` == 1, 'Ever Progressor', 'Never Progressor')]

# factor timepoint
finalData[,`Sample Time Point` := factor(`Sample Time Point`, levels=mixedsort(unique(`Sample Time Point`)))]

# keep only samples with T0 and another timepoint
baselineSamples <- finalData[`Sample Time Point` == 'T0']$`Subject ID`
multiTimePointSamples <- finalData[,.N,by=.(`Subject ID`)][N > 1]$`Subject ID`

keep <- baselineSamples[baselineSamples %in% multiTimePointSamples]
finalData <- finalData[`Subject ID` %chin% keep]

################################################################################
################## run stats ###################################################

# switch statement to re-define p-values
relabel_p <- function(x){
  if(x < .0001) return('p < 0.0001')
  if(x < .001) return('p < 0.001')
  if(x < .01) return('p < 0.01')
  if(x < .05) return('p < 0.05')
  return(paste0('p = ', signif(x, 3)))
}

# stats
statTest <- wilcox.test(finalData[`Sample Time Point` != 'T0' & progressionStatus == 'Never Progressor', median(`DELFI-TF`, na.rm=T), by=.(`Subject ID`)]$V1,
                        finalData[`Sample Time Point` != 'T0' & progressionStatus != 'Never Progressor', median(`DELFI-TF`, na.rm=T), by=.(`Subject ID`)]$V1)

statDT <- data.table(x='T11',
                     y=.38,
                     label=paste0(relabel_p(statTest$p.value)),
                     progressionStatus="Ever Progressor")

################################################################################
################### Make Plot ##################################################

# make plot
FigS8A <- ggplot() +
  geom_line(data=finalData, mapping=aes(x=`Sample Time Point`, y=`DELFI-TF`, color=progressionStatus, group=`Subject ID`), alpha=.75, linewidth=1) +
  geom_text(data=statDT, aes(x=x, y=y, label=label), size=10) +
  facet_wrap(~progressionStatus, ncol=1) +
  scale_y_continuous(trans="sqrt", labels=scales::percent) +
  geom_hline(yintercept=.006, color="grey70", linetype=5, linewidth=1) +
  scale_color_manual("", values=c("#C93312", "#81A88D")) +
  xlab("Blood Collection Time Point") +
  facet_wrap(~progressionStatus, ncol=1) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(axis.text.x=element_text(size=28, angle=0), axis.text.y=element_text(size=28)) +
  theme(axis.title=element_text(size=34)) +
  theme(legend.position = "top") +
  theme(legend.title=element_text(size=28), legend.text=element_text(size=28))+
  theme(strip.text=element_text(size=28, color='white')) +
  theme(strip.background = element_rect(fill="#1B465A"))

# save plots
png(filename="Figures/FigS8A.png", height=10, width=14, units="in", res=150)
FigS8A
dev.off()

cairo_pdf(file="Figures/FigS8A.pdf", height=10, width=18)
FigS8A
dev.off()
