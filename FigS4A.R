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
keep <- c('Sequence ID', 'Subject ID', 'Sample Time Point', 'DELFI-TF', 'MAF', 'ichorCNA')
BloodSampleCharacteristics <- BloodSampleCharacteristics[,..keep]

# clinical data
keep <- c('Subject ID', 'Study Arm')
ClinicalDemographics <- ClinicalDemographics[,..keep]

# merge data together
finalData <- merge(BloodSampleCharacteristics, ClinicalDemographics, by='Subject ID', all.x=T)

################################################################################
################### final formatting ###########################################

# restrict to only values with MAF, DELFI-TF, and ichorCNA
finalData <- finalData[(!is.na(ichorCNA) & !is.na(`DELFI-TF`) & !is.na(MAF))]

# calculate cut points
mafQuant <- quantile(finalData[MAF != 0]$MAF)

# select only those values where MAF is 0 "undetectable" and also above Q2 MAF from the boxplots in FigS3A, close to the apparent ichorCNA floor
finalData <- finalData[MAF == 0 & ichorCNA >= mafQuant[3]]

################################################################################
############## run correlation #################################################

# switch statement to re-define p-values
relabel_p <- function(x){
  if(x < .0001) return('p < 0.0001')
  if(x < .001) return('p < 0.001')
  if(x < .01) return('p < 0.01')
  if(x < .05) return('p < 0.05')
  return(paste0('p = ', signif(x, 3)))
}

# run cor
corStat <- cor.test(finalData$`DELFI-TF`, finalData$ichorCNA)

# make label
corStat <- paste0("r=", signif(corStat$estimate, 2), "\n", relabel_p(corStat$p.value))

################################################################################
############# make final plot ##################################################

# scatter plot
FigS4A <- ggplot(finalData, aes(x=ichorCNA, y=`DELFI-TF`)) +
  geom_point(color="#d48e20", alpha=.75, size=4) +
  geom_smooth(method="lm", color="black") +
  scale_y_sqrt(label=scales::percent, limits=c(0, .5), breaks=c(0, .05, .1, .2, .3, .4, .5)) +
  scale_x_sqrt(label=scales::percent, breaks=c(.05, .1, .15, .2, .25, .3)) +
  annotate("text", x=.1, y=.4, label=corStat, size=10) +
  facet_wrap(~'Undetectable by MAF') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(strip.text=element_text(color="white", size=28), strip.background = element_rect(fill="#1B465A")) +
  theme(axis.text = element_text(size=24), axis.title=element_text(size=28))

# save plot
png(filename="Figures/FigS4A.png", height=10, width=12, units="in", res=150)
FigS4A
dev.off()

cairo_pdf(file="Figures/FigS4A.pdf", height=10, width=12)
FigS4A
dev.off()

