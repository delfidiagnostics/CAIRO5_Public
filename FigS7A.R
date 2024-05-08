################################################################################
#################### Setup #####################################################

# loc
setwd("~/git/CAIRO5_Public")

# libraries
library(readxl)
library(data.table)
library(ggplot2)

# Data
BloodSampleCharacteristics <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=2, skip=1))
ClinicalDemographics <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=1, skip=1))
SumOfLongestDiameters <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=7, skip=1))

################################################################################
################# pair down data to required cols ##############################

# Blood Samples
keep <- c('Sequence ID', 'Subject ID', 'Sample Time Point', 'DELFI-TF', 'MAF')
BloodSampleCharacteristics <- BloodSampleCharacteristics[,..keep]

# clinical data
keep <- c('Subject ID', 'Study Arm', 'Progression Status (0=no progression; 1=progression)')
ClinicalDemographics <- ClinicalDemographics[,..keep]

# Sum of Longest Diameter, pair down to be within 60 days
SumOfLongestDiameters <- SumOfLongestDiameters[abs(`Time Between SLD and Liquid Biopsy (Days)`) <= 60]
SumOfLongestDiameters[,`Time Between SLD and Liquid Biopsy (Days)` := NULL]

# merge data together
finalData <- merge(BloodSampleCharacteristics, ClinicalDemographics, by='Subject ID', all.x=T)
finalData <- merge(finalData, SumOfLongestDiameters, by=c('Subject ID', 'Sample Time Point'), all.x=T)

################################################################################
################### final formatting ###########################################

# restrict to ARM 1 samples at baseline
finalData <- finalData[`Study Arm` == 1 & `Sample Time Point` == 'T0']

# name progression status for plot
finalData[,progressionStatus := ifelse(`Progression Status (0=no progression; 1=progression)` == 1, 'Ever Progressor', 'Never Progressor')]

# change col name to fit on plot later
setnames(finalData, "Sum of the Largest Diameters", "Sum of the\nLongest Diameters")

# melt the data into a long format for plotting
finalData <- melt(finalData, measure.vars=c("DELFI-TF", "MAF", "Sum of the\nLongest Diameters"))

################################################################################
################## run the stats ###############################################

# switch statement to re-define p-values
relabel_p <- function(x){
  if(x < .0001) return('****')
  if(x < .001) return('***')
  if(x < .01) return('**')
  if(x < .05) return('*')
  return('ns')
}

# run stats
delfiStat <- wilcox.test(finalData[variable == 'DELFI-TF' & `Progression Status (0=no progression; 1=progression)` == 1]$value,
                         finalData[variable == 'DELFI-TF' & `Progression Status (0=no progression; 1=progression)` == 0]$value)
mafStat <- wilcox.test(finalData[variable == 'MAF' & `Progression Status (0=no progression; 1=progression)` == 1]$value,
                         finalData[variable == 'MAF' & `Progression Status (0=no progression; 1=progression)` == 0]$value)
sldStat <- wilcox.test(finalData[variable == 'Sum of the\nLongest Diameters' & `Progression Status (0=no progression; 1=progression)` == 1]$value,
                       finalData[variable == 'Sum of the\nLongest Diameters' & `Progression Status (0=no progression; 1=progression)` == 0]$value)


# make dt to plot
statData <- data.table("variable"=c('DELFI-TF', 'MAF', 'Sum of the\nLongest Diameters'),
           "x"=c(1.5),
           "y"=c(.7, .95, 280),
           "p_value"=c(relabel_p(delfiStat$p.value), relabel_p(mafStat$p.value), relabel_p(sldStat$p.value)))
segmentData <- data.table("variable"=c('DELFI-TF', 'MAF', 'Sum of the\nLongest Diameters'),
                          "x"=c(1),
                          "xend"=c(2),
                          "y"=c(.65, .9, 260),
                          "yend"=c(.65, .9, 260))

################################################################################
################### Construct the Violin plot ##################################

# violin
FigS7A <- ggplot() +
  geom_boxplot(data=finalData, mapping=aes(x=progressionStatus, y=value, fill=progressionStatus), outlier.shape=NA) +
  geom_point(data=finalData, mapping=aes(x=progressionStatus, y=value, fill=progressionStatus), position=position_jitter(width=.2), size=4, show.legend = FALSE, color='grey30', alpha=.75) +
  geom_segment(data=segmentData, mapping=aes(x=x, xend=xend, y=y, yend=yend), linewidth=1.75, linetype='solid', color='gray65', lineend = "round") +
  geom_text(data=statData, mapping=aes(x=x, y=y, label=p_value), size=14) +
  scale_fill_manual("Progression Status", values=c("#81A88D", "#C93312")) +
  facet_wrap(~variable, scales="free_y", ncol=1, strip.position = "right") +
  xlab("Progression Status") + ylab("") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(strip.text=element_text(color="white", face="bold", size=42)) +
  theme(axis.text.x=element_text(size=36), axis.text.y=element_text(size=36)) +
  theme(axis.title=element_text(size=42)) +
  theme(strip.background=element_rect(fill="#1B465A")) +
  theme(panel.grid=element_blank())

# save plots
png(filename="Figures/FigS7A.png", height=17, width=17, units="in", res=150)
FigS7A
dev.off()

cairo_pdf(file="Figures/FigS7A.pdf", height=17, width=17)
FigS7A
dev.off()


