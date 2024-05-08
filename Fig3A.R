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
auxRecistData <- readRDS('Data/Fig3A_AuxData.rds')

################################################################################
################# pair down data to required cols ##############################

# Blood Samples
keep <- c('Sequence ID', 'Subject ID', 'Sample Time Point', 'DELFI-TF', 'MAF')
BloodSampleCharacteristics <- BloodSampleCharacteristics[,..keep]

# clinical data
keep <- c('Subject ID', 'Study Arm')
ClinicalDemographics <- ClinicalDemographics[,..keep]

# merge data together
finalData <- merge(BloodSampleCharacteristics, ClinicalDemographics, by='Subject ID', all.x=T)

################################################################################
############## format data #####################################################

# restrict to arm 1 and Baseline
finalData <- finalData[`Study Arm` == 1 & `Sample Time Point` == 'T0']

# add on the RECIST categories
finalData <- merge(finalData, auxRecistData, by='Sequence ID', all.x=TRUE)

# anything with an NA RECIST had no consecutive RECIST and should be removed
finalData <- finalData[!is.na(`Best`)]

# melt data
finalData <- melt(finalData, measure.vars=c("DELFI-TF", "MAF"))

# factor for consistent order
finalData[,variable := factor(variable, levels=c('DELFI-TF', 'MAF'))]

################################################################################
############# Run Stats ########################################################

# switch statement to re-define p-values
relabel_p <- function(x){
  if(x < .0001) return('****')
  if(x < .001) return('***')
  if(x < .01) return('**')
  if(x < .05) return('*')
  return('ns')
}

# run wilcox
delfiStat <- wilcox.test(finalData[variable == "DELFI-TF" & Category == "SD/PD"]$value, finalData[variable == "DELFI-TF" & Category == "CR/PR"]$value)
delfiStat <- relabel_p(delfiStat$p.value)
mafStat <- wilcox.test(finalData[variable == "MAF" & Category == "SD/PD"]$value, finalData[variable == "MAF" & Category == "CR/PR"]$value)
mafStat <- relabel_p(mafStat$p.value)

# labelDT
statDT <- data.table(variable=c("DELFI-TF", "MAF"),
                     y=c(.93, .93),
                     label=c(delfiStat, mafStat))

################################################################################
########### Plot ###############################################################

Fig3A <- ggplot() +
  geom_boxplot(data=finalData, mapping=aes(x=variable, y=value, fill=Category), outlier.shape = NA) +
  geom_point(data=finalData, mapping=aes(x=variable, y=value, fill=Category, color=Best), alpha=.85, position=position_jitterdodge(jitter.width=.1, dodge.width=.75), size=7) +
  geom_text(data=statDT, mapping=aes(x=variable, y=y, label=label), size=16, color='grey30') +
  annotate("segment", x=.75, xend=1.25, y=.9, yend=.9, linewidth=1.75, linetype='solid', color='gray65', lineend = "round") +
  annotate("segment", x=1.75, xend=2.25, y=.9, yend=.9, linewidth=1.75, linetype='solid', color='gray65', lineend = "round") +
  scale_y_continuous(labels = scales::percent) +
  ylab("Tumor Fraction") +
  scale_fill_manual(values=c("CR/PR"="#F1BB7B", "SD/PD"="#e6e8fa")) +
  scale_color_manual("Best\nResponse",values=c("CR"="#8fbc8f", "PR"="#004953", "PD"="#ab4a37", 'SD'='#704214')) +
  guides(color=guide_legend(ncol=2), fill="none") +
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
png(filename="Figures/Fig3A.png", height=16, width=13, units="in", res=150)
Fig3A
dev.off()

cairo_pdf(file="Figures/Fig3A.pdf", height=16, width=13)
Fig3A
dev.off()




