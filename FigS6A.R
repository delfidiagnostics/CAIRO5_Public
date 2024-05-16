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
keep <- c('Subject ID', 'Study Arm')
ClinicalDemographics <- ClinicalDemographics[,..keep]

# Sum of Longest Diameter, pair down to be within 60 days
SumOfLongestDiameters <- SumOfLongestDiameters[abs(`Time Between SLD and Liquid Biopsy (Days)`) <= 60]
SumOfLongestDiameters[,`Time Between SLD and Liquid Biopsy (Days)` := NULL]

# merge data together
finalData <- merge(BloodSampleCharacteristics, ClinicalDemographics, by='Subject ID', all.x=T)
finalData <- merge(finalData, SumOfLongestDiameters, by=c('Subject ID', 'Sample Time Point'), all.x=T)

################################################################################
############### format for plot ################################################

# restrict to arm1 baseline samples, and make sure everything has a value
finalData <- finalData[`Study Arm` == 1 & `Sample Time Point` == 'T0']
finalData <- finalData[!(is.na(`DELFI-TF`) | is.na(MAF) | is.na(`Sum of the Longest Diameters`))]

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
delfiStat <- cor.test(finalData[variable == 'DELFI-TF']$value, finalData[variable == 'DELFI-TF']$`Sum of the Longest Diameters`)
mafStat <- cor.test(finalData[variable == 'MAF']$value, finalData[variable == 'MAF']$`Sum of the Longest Diameters`)

# labels
delfiStat <- paste0("r=", signif(delfiStat$estimate, 2), "\n", relabel_p(delfiStat$p.value))
mafStat <- paste0("r=", signif(mafStat$estimate, 2), "\n", relabel_p(mafStat$p.value))

# make dt to plot
statDT <- data.table(variable=c('DELFI-TF', 'MAF'),
                     x=c(185),
                     y=c(.65, .8),
                     label=c(delfiStat, mafStat))

################################################################################
############ Construct Plot ####################################################

# plot
FigS6A <- ggplot(finalData, aes(x=`Sum of the Longest Diameters`, y=value, color=variable)) +
  geom_point(size=5, alpha=.65) +
  geom_smooth(method="lm", alpha=.1) +
  geom_text(data=statDT, mapping=aes(x=x, y=y, label=label, color=variable), size=14, key_glyph="point") +
  scale_color_manual(values=c("DELFI-TF"="#046C9A", "MAF"="#D69C4E")) +
  scale_y_continuous("Mutation Allele Frequency", labels=scales::percent, sec.axis = sec_axis(trans=~., name="DELFI-TF", labels=scales::percent)) +
  facet_wrap(~'Baseline') +
  xlab("Sum of The Longest Diameters (mm)") +
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
png(filename="Figures/FigS6A.png", height=14, width=16, units="in", res=150)
FigS6A
dev.off()

cairo_pdf(file="Figures/FigS6A.pdf", height=14, width=16)
FigS6A
dev.off()



