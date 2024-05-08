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

################################################################################
################# pair down data to required cols ##############################

# Blood Samples
keep <- c('Sequence ID', 'Subject ID', 'Sample Time Point', 'Days From Randomization', 'DELFI-TF', 'MAF')
BloodSampleCharacteristics <- BloodSampleCharacteristics[,..keep]

# clinical data
keep <- c('Subject ID', 'Study Arm', 'Progression Status (0=no progression; 1=progression)', 'PFS Registration')
ClinicalDemographics <- ClinicalDemographics[,..keep]

# merge data together
finalData <- merge(BloodSampleCharacteristics, ClinicalDemographics, by='Subject ID', all.x=T)

################################################################################
################# format data for plotting #####################################

# restrict to arm 1
finalData <- finalData[`Study Arm` == 1]

# remove samples without a T0 baseline
keep <- finalData[`Sample Time Point` == 'T0']$`Subject ID`
finalData <- finalData[`Subject ID` %chin% keep]

# remove samples with only 2 observations
remove <- finalData[,.N,by=.(`Subject ID`)][N <= 2]$`Subject ID`
finalData <- finalData[!`Subject ID` %chin% remove]

# convert days to months for plot
daysInMonth <- 30.4167
finalData[,`Months From Registration` := `Days From Randomization`/daysInMonth]
finalData[,`Months From Registration PFS` := `PFS Registration`/daysInMonth]

# split out separate DT for progression lines
progessionDT <- finalData[`Progression Status (0=no progression; 1=progression)` == 1,.(`Subject ID`, `Months From Registration PFS`)]
progessionDT <- unique(progessionDT)

# melt data
finalData <- melt(finalData, measure.vars = c('DELFI-TF', 'MAF'))

################################################################################
############# make plots #######################################################

# make plot
FigS9A <- ggplot() +
  geom_vline(data=progessionDT, mapping=aes(xintercept=`Months From Registration PFS`), linetype=2, color='gray60', linewidth=1) +
  geom_line(data=finalData, mapping=aes(x=`Months From Registration`, y=value, color=variable, group=variable)) +
  geom_point(data=finalData, mapping=aes(x=`Months From Registration`, y=value, color=variable)) +
  scale_color_manual("",values=c("DELFI-TF"="#046C9A", "MAF"="#D69C4E")) +
  facet_wrap(~`Subject ID`, scales = 'free_x', ncol=7) +
  xlab("Months From Registration") +
  ylab("Score") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(strip.text=element_text(color='white', face='bold', size=20), strip.background = element_rect(fill="#1B465A")) +
  theme(legend.position = 'top') +
  theme(legend.key.size = unit(1, 'cm')) +
  theme(legend.text=element_text(size=18)) +
  theme(axis.text=element_text(size=18)) +
  theme(axis.title=element_text(size=20))

# save plot
png(filename="Figures/FigS9A.png", height=18, width=14, units="in", res=150)
FigS9A
dev.off()

cairo_pdf(file="Figures/FigS9A.pdf", height=18, width=14)
FigS9A
dev.off()
