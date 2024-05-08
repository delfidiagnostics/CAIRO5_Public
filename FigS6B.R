################################################################################
#################### Setup #####################################################

# loc
setwd("~/git/CAIRO5_Public")

# libraries
library(readxl)
library(data.table)
library(ggplot2)
library(scales)
library(ggallin)

# Data
BloodSampleCharacteristics <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=2, skip=1))
ClinicalDemographics <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=1, skip=1))

################################################################################
################# pair down data to required cols ##############################

# Blood Samples
keep <- c('Sequence ID', 'Subject ID', 'Sample Time Point', 'DELFI-TF', 'MAF', 'CEA')
BloodSampleCharacteristics <- BloodSampleCharacteristics[,..keep]

# clinical data
keep <- c('Subject ID', 'Study Arm')
ClinicalDemographics <- ClinicalDemographics[,..keep]

# merge data together
finalData <- merge(BloodSampleCharacteristics, ClinicalDemographics, by='Subject ID', all.x=T)

################################################################################
############### format for plot ################################################

# restrict to arm1 baseline samples, and make sure everything has a value
finalData <- finalData[`Study Arm` == 1 & `Sample Time Point` == 'T0']
finalData <- finalData[!(is.na(`DELFI-TF`) | is.na(MAF) | is.na(CEA))]

# melt the data into the expected format
finalData <- melt(finalData, measure.vars=c("DELFI-TF", "MAF"))

################################################################################
################ run stats #####################################################

# switch statement to re-define p-values
relabel_p <- function(x){
  if(x < .0001) return('p < 0.0001')
  if(x < .001) return('p < 0.001')
  if(x < .01) return('p < 0.01')
  if(x < .05) return('p < 0.05')
  return(paste0('p = ', signif(x, 3)))
}

# correlation
delfiStat <- cor.test(finalData[variable == 'DELFI-TF']$value, finalData[variable == 'DELFI-TF']$CEA)
mafStat <- cor.test(finalData[variable == 'MAF']$value, finalData[variable == 'MAF']$CEA)

# labels
delfiStat <- paste0("r=", signif(delfiStat$estimate, 2), "\n", relabel_p(delfiStat$p.value))
mafStat <- paste0("r=", signif(mafStat$estimate, 2), "\n", relabel_p(mafStat$p.value))

# make dt to plot
statDT <- data.table(variable=c('DELFI-TF', 'MAF'),
                     x=c(400),
                     y=c(.67, .83),
                     label=c(delfiStat, mafStat))

################################################################################
############ Construct Plot ####################################################

# plot
FigS6B <- ggplot(finalData, aes(x=CEA, y=value, color=variable)) +
  geom_point(size=5, alpha=.65) +
  geom_smooth(method="lm", alpha=.1) +
  geom_text(data=statDT, mapping=aes(x=x, y=y, label=label, color=variable), size=14, key_glyph="point") +
  scale_color_manual(values=c("DELFI-TF"="#046C9A", "MAF"="#D69C4E")) +
  scale_y_continuous("Mutation Allele Frequency", sec.axis = sec_axis(trans=~., name="DELFI-TF", labels=scales::percent), labels=scales::percent) +
  scale_x_continuous(oob=squish, limits=c(0, 500)) +
  xlab("Carcinoembryonic Antigen (μg/L)") +
  facet_wrap(~'Baseline') +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid=element_blank()) +
  theme(axis.text.x=element_text(size=36), axis.text.y=element_text(size=36)) +
  theme(axis.title.x=element_text(size=42)) +
  theme(axis.title.y=element_text(size=42, color="#D69C4E")) +
  theme(axis.title.y.right=element_text(size=42, color="#046C9A")) +
  theme(strip.background=element_rect(fill="#1B465A")) +
  theme(strip.text=element_text(color="white", size=42, face="bold"))

# save plot
png(filename="Figures/FigS6B.png", height=14, width=16, units="in", res=150)
FigS6B
dev.off()

cairo_pdf(file="Figures/FigS6B.pdf", height=14, width=16)
FigS6B
dev.off()



