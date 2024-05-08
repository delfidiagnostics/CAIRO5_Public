################################################################################
#################### Setup #####################################################

# loc
setwd("~/git/CAIRO5_Public")

# libraries
library(readxl)
library(data.table)
library(ggplot2)
library(gridExtra)
library(grid)
library(gtable)
library(ggallin)

# Data
BloodSampleCharacteristics <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=2, skip=1))
ClinicalDemographics <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=1, skip=1))
Slopes <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=8, skip=1))

################################################################################
################# pair down data to required cols ##############################

# Blood Samples
keep <- c('Sequence ID', 'Subject ID', 'Sample Time Point', 'Days From Randomization')
BloodSampleCharacteristics <- BloodSampleCharacteristics[,..keep]

# clinical data
keep <- c('Subject ID', 'OS Registration')
ClinicalDemographics <- ClinicalDemographics[,..keep]

# merge data together
finalData <- merge(ClinicalDemographics, Slopes, by='Subject ID', all.x=T)

################################################################################
################ format for plot ###############################################

# restrict down to eligible subjects, subjects 11 and 23 are ineligible
finalData <- finalData[!is.na(Velocity)]
finalData <- finalData[!`Subject ID` %in% c(11, 123)]

# annotate slope direction
medianVelocity <- quantile(finalData$Velocity)[3]
finalData[,slopeThreshold := ifelse(Velocity > medianVelocity, "Above Median", "Below Median")]

# determine sample order based on OS
sampleOrder <- finalData[order(slopeThreshold, `OS Registration`)]$`Subject ID`
finalData <- finalData[,`Subject ID` := factor(`Subject ID`, levels=sampleOrder)]

# transform days to weeks for OS
finalData[,`OS Registration` := `OS Registration`/7]

################################################################################
############# format additional data for biopsies times ########################

# convert to weeks
BloodSampleCharacteristics[,`Weeks From Registration` := `Days From Randomization`/7]

# restrict relevant sampels and factor them
BloodSampleCharacteristics <- BloodSampleCharacteristics[`Subject ID` %in% sampleOrder]
BloodSampleCharacteristics <- BloodSampleCharacteristics[,`Subject ID` := factor(`Subject ID`, levels=sampleOrder)]

################################################################################
################# plots ########################################################

# make the swimmer plot
p1 <- ggplot() +
  geom_segment(data=finalData, mapping=aes(y=`Subject ID`, yend=`Subject ID`, x=0, xend=`OS Registration`), linewidth=4.5, alpha=.7, lineend="round", color="#CB2314") +
  geom_point(data=BloodSampleCharacteristics, mapping=aes(x=`Weeks From Registration`, y=`Subject ID`), shape=1, size=4.5, stroke=1.5) +
  xlab("Weeks On Study") + 
  scale_x_continuous(limits=c(0, 310), expand=c(0, 0), breaks=seq(0, 300, by=50)) +
  theme_bw() +
  theme(axis.title.y=element_blank()) +
  theme(panel.grid = element_blank()) +
  theme(axis.text.x=element_text(size=24)) +
  theme(axis.text.y=element_text(size=24)) +
  theme(axis.title = element_text(size=32)) +
  theme(legend.text=element_text(size=24)) +
  theme(legend.title = element_text(size=32)) +
  theme(legend.position = "top") +
  theme(legend.box="vertical") +
  theme(plot.margin = margin(t = 0, r = 15, b = 0, l = 10, unit = "pt"))

# make the slope plot
p2 <- ggplot() +
  geom_vline(xintercept = medianVelocity, linetype=2, linewidth=2, color='grey70') +
  geom_point(data=finalData, mapping=aes(x=Velocity, y=`Subject ID`, color=slopeThreshold), size=4.5) +
  scale_color_manual(values=c("Above Median"="#da7b3d", "Below Median"="#8da8bf")) +
  scale_x_continuous(trans=ssqrt_trans, breaks=c(-0.001, 0, 0.001)) +
  xlab("DELFI-TF\nSlope") +
  theme_bw() +
  theme(axis.title.y=element_blank()) +
  theme(legend.position = "none") +
  theme(panel.grid = element_blank()) +
  theme(axis.text.x=element_text(size=24)) +
  theme(axis.text.y=element_blank()) + 
  theme(axis.ticks.y=element_blank()) +
  theme(axis.title = element_text(size=32)) +
  theme(plot.margin = margin(t = 0, r = 9, b = 0, l = 15, unit = "pt"))

# convert to grob
p1 <- ggplotGrob(p1)
p2 <- ggplotGrob(p2)

# align plot heights
maxHeight <- unit.pmax(p1$heights, p2$heights)
p1$heights <- maxHeight
p2$heights <- maxHeight

# arrange the plot
finalPlot <- arrangeGrob(p2, p1, ncol=2,widths=c(.15, .85))

# save the results
png(filename="Figures/Fig4A.png", height=22, width=32, units="in", res=150)
grid.draw(finalPlot)
dev.off()

cairo_pdf(file="Figures/Fig4A.pdf", height=22, width=32)
grid.draw(finalPlot)
dev.off()



