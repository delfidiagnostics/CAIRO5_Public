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
ClinicalDemographics <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=1, skip=1))

################################################################################
################# pair down data to required cols ##############################

# Blood Samples
keep <- c('Sequence ID', 'Subject ID', 'Sample Time Point', 'DELFI-TF', 'MAF')
BloodSampleCharacteristics <- BloodSampleCharacteristics[,..keep]

# clinical data
keep <- c('Subject ID', 'Study Arm', 'Surgery (1=Yes; 0=No)', 'Complete Resection (1=Yes; 0=No)')
ClinicalDemographics <- ClinicalDemographics[,..keep]

# merge data together
finalData <- merge(BloodSampleCharacteristics, ClinicalDemographics, by='Subject ID', all.x=T)

################################################################################
################## Format Data #################################################

# restrict to arm 1 and Baseline
finalData <- finalData[`Study Arm` == 1 & `Sample Time Point` == 'T0']

# label the resection status
finalData[,`Complete Resection (1=Yes; 0=No)` := as.numeric(`Complete Resection (1=Yes; 0=No)`)]
finalData[,`Complete Resection Status` := ifelse(`Complete Resection (1=Yes; 0=No)` == 1, "Complete\nResection", "Incomplete\nResection")]
finalData[is.na(`Complete Resection Status`), `Complete Resection Status` := "No Resection"]
finalData[,'Complete Resection Status' := factor(`Complete Resection Status`, levels=c("Complete\nResection", "Incomplete\nResection", "No Resection"))]

################################################################################
################# run Stats ####################################################

# switch statement to re-define p-values
relabel_p <- function(x){
  if(x < .0001) return('****')
  if(x < .001) return('***')
  if(x < .01) return('**')
  if(x < .05) return('*')
  return('ns')
}

# run kruskal-wallis test
delfiLabel <- kruskal.test(finalData$`DELFI-TF`, finalData$`Complete Resection Status`)
delfiLabel <- relabel_p(delfiLabel$p.value)
mafLabel <- kruskal.test(finalData$MAF, finalData$`Complete Resection Status`)
mafLabel <- relabel_p(mafLabel$p.value)

# labelDT
statDT <- data.table(variable=c("DELFI-TF", "MAF"),
                      y=c(.93, .93),
                      label=c(delfiLabel, mafLabel))

################################################################################
############## Make Plot #######################################################

# melt data for plot
finalData <- melt(finalData, measure.vars = c('DELFI-TF', 'MAF'))

# make the plot
Fig3B <- ggplot() +
  geom_boxplot(data=finalData, mapping=aes(x=variable, fill=`Complete Resection Status`, y=value), position=position_dodge(), outlier.shape = NA) +
  geom_point(data=finalData, mapping=aes(x=variable, fill=`Complete Resection Status`, y=value), position=position_jitterdodge(jitter.width = .2), size=5, alpha=.75, color='grey30', show.legend = FALSE) +
  geom_text(data=statDT, aes(x=variable, y=y, label=label), size=16) +
  annotate("segment", x=.75, xend=1.25, y=.9, yend=.9, linewidth=1.75, linetype='solid', color='gray65', lineend = "round") +
  annotate("segment", x=1.75, xend=2.25, y=.9, yend=.9, linewidth=1.75, linetype='solid', color='gray65', lineend = "round") +
  scale_y_continuous(labels = scales::percent) +
  ylab("Tumor Fraction") +
  scale_fill_manual("Resection\nStatus",values=c("Complete\nResection"="#DD8D29", "Incomplete\nResection"="#E2D200", "No Resection"="#46ACC8")) +
  guides(fill=guide_legend(ncol=2)) +
  theme_bw() +
  theme(legend.position = "top") +
  theme(legend.title=element_text(size=40), legend.text=element_text(size=36)) +
  theme(legend.key.size = unit(2, 'cm')) +
  theme(legend.spacing.y = unit(2, 'cm')) +
  theme(axis.text.x=element_text(size=36), axis.text.y=element_text(size=36)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_text(size=40)) +
  theme(strip.background=element_rect(fill="black")) +
  theme(strip.text=element_text(color="white", size=36, face="bold")) +
  theme(panel.grid=element_blank())

# save plot
png(filename="Figures/Fig3B.png", height=16, width=13, units="in", res=150)
Fig3B
dev.off()

cairo_pdf(file="Figures/Fig3B.pdf", height=16, width=13)
Fig3B
dev.off()
