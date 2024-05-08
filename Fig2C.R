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
################# Combine Data #################################################

# Blood Samples
keep <- c('Sequence ID', 'Subject ID', 'Sample Time Point', 'DELFI-TF', 'MAF')
BloodSampleCharacteristics <- BloodSampleCharacteristics[,..keep]

# clinical data
keep <- c('Subject ID', 'Study Arm')
ClinicalDemographics <- ClinicalDemographics[,..keep]

# merge data together
finalData <- merge(BloodSampleCharacteristics, ClinicalDemographics, by='Subject ID', all.x=T)

################################################################################
################## format final data ############################################

# restrict to arm 1 and make sure we have MAF and DELFI-TF values
finalData <- finalData[`Study Arm` == 1]
finalData <- finalData[!(is.na(`DELFI-TF`) | is.na(MAF))]

# create label for undetectable by MAF
finalData[,StatusMAF := ifelse(MAF == 0, 'Undetectable', 'Detectable')]

################################################################################
################# Run correlations for labels ##################################

# switch statement to re-define p-values
relabel_p <- function(x){
  if(x < .0001) return('p < 0.0001')
  if(x < .001) return('p < 0.001')
  if(x < .01) return('p < 0.01')
  if(x < .05) return('p < 0.05')
  return(paste0('p = ', signif(x, 3)))
}

# correlation of MAF and DELFI-TF excluding undetectable samples
corStat <- cor.test(finalData[StatusMAF == 'Detectable']$MAF, finalData[StatusMAF == 'Detectable']$`DELFI-TF`, method='pearson')

# label p-value
p_value <- signif(corStat$p.value, 3)
label <- paste0("r = ", signif(corStat$estimate, 2), "\n", relabel_p(corStat$p.value))
corLabel <- data.table("x"=.6, "y"=.09, "label"=label)

################################################################################
################ Simple plot ###################################################

# make plot without zoom
p2 <- ggplot() +
  geom_jitter(data=finalData[StatusMAF == 'Undetectable'], size=3, aes(color=StatusMAF, x=MAF, y=`DELFI-TF`), width=.001) +
  geom_point(data=finalData[StatusMAF == 'Detectable'], size=3, aes(color=StatusMAF, x=MAF, y=`DELFI-TF`)) +
  geom_smooth(data=finalData[StatusMAF == 'Detectable'], aes(x=MAF, y=`DELFI-TF`), method="glm", color="gray20") +
  geom_text(data=corLabel, mapping=aes(x=x, y=y, label=label), size=14) +
  scale_y_continuous(trans=ssqrt_trans, limits=c(0, .9), labels = scales::percent, breaks=c(.025, .1, .25, .4, .6, .8)) +
  scale_x_continuous(trans=ssqrt_trans, labels = scales::percent, breaks=c(0, .025, .1, .25, .4, .6, .8)) +
  xlab("Mutant Allele Frequency") + ylab("DELFI-TF") +
  scale_color_manual(values=c("#81A88D", "#972D15")) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(axis.text.x=element_text(size=36), axis.text.y=element_text(size=36)) +
  theme(axis.title=element_text(size=42)) +
  theme(legend.position = "none")

# save the plot
png(filename="Figures/Fig2C.png", height=14, width=16, units="in", res=150)
p2
dev.off()

cairo_pdf(file="Figures/Fig2C.pdf", height=14, width=16)
p2
dev.off()
