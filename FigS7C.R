################################################################################
#################### Setup #####################################################

# loc
setwd("~/git/CAIRO5_Public")

# libraries
library(readxl)
library(data.table)
library(ggplot2)
library(survival)
library(survminer)

# Data
BloodSampleCharacteristics <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=2, skip=1))
ClinicalDemographics <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=1, skip=1))

################################################################################
################# pair down data to required cols ##############################

# Blood Samples
keep <- c('Sequence ID', 'Subject ID', 'Sample Time Point', 'CEA')
BloodSampleCharacteristics <- BloodSampleCharacteristics[,..keep]

# clinical data
keep <- c('Subject ID', 'Study Arm', 'PFS Registration', 'PFS Status (0=alive; 1=progression or death)')
ClinicalDemographics <- ClinicalDemographics[,..keep]

# merge data together
finalData <- merge(BloodSampleCharacteristics, ClinicalDemographics, by='Subject ID', all.x=T)

################################################################################
############### format data for analysis #######################################

# restrict to ARM 1 samples at baseline
finalData <- finalData[`Sample Time Point` == 'T0']

# convert OS from supp. to days
daysInMonth <- 30.4167
finalData[,`PFS Registration` := `PFS Registration`/daysInMonth]

# classify high/low delfi-tf
ceaQuantile <- quantile(finalData$CEA, na.rm=T)
finalData[,CEA := ifelse(CEA > ceaQuantile[3], "High CEA", "Low CEA")]

# factor
factorOrder <- c("Low CEA", "High CEA")
finalData[,CEA := factor(CEA, levels=factorOrder)]

################################################################################
########### perform survival analysis ##########################################

# switch statement to re-define p-values
relabel_p <- function(x){
  if(x < .0001) return('p < 0.0001')
  if(x < .001) return('p < 0.001')
  if(x < .01) return('p < 0.01')
  if(x < .05) return('p < 0.05')
  return(paste0('p = ', signif(x, 3)))
}

# do the survival fit
fit <- survfit(Surv(`PFS Registration`, `PFS Status (0=alive; 1=progression or death)`) ~ CEA, data = finalData)

# calculate hazard ratio
cox <-coxph(Surv(`PFS Registration`, `PFS Status (0=alive; 1=progression or death)`) ~ CEA, data = finalData)
cox <- summary(cox)


# label p-value
p_value <- relabel_p(cox$sctest[3])
p_value <- paste0("Log-Rank\n", p_value)

################################################################################
######### Make plots ###########################################################

# make themes for the plot
tableTheme <- theme_void() +
  theme(axis.text.y=element_text(size=28), plot.title = element_blank())

mainTheme <- theme_classic() +
  theme(axis.text=element_text(size=28), axis.title=element_text(size=34)) +
  theme(legend.text=element_text(size=34), legend.title = element_blank())

# make the plot
FigS7C <- ggsurvplot(
  fit,
  data = finalData,
  xlab = "Time in Months",
  xlim=c(0,36),
  legend=c(.75, .75),
  ylab = "Progression-Free Survival",
  break.time.by = 6,
  size = 3,                 # change line size
  palette = c("#8cd29c", "#c16f84"), # custom color palettes
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs = factorOrder,
  risk.table.height = 0.2, # Useful to change when you have multiple groups
  ggtheme = mainTheme, # Change ggplot2 theme
  tables.theme = tableTheme,
  surv.median.line = "hv",
  fontsize=15,
  censor.shape="|",
  censor.size=14,
  risk.table.y.text=FALSE
)

# annotate Hazzard ratio
FigS7C$plot <- FigS7C$plot+
  ggplot2::annotate("text",
                    x = 4, y = 0.25, # x and y coordinates of the text
                    label = p_value, size = 15) +
  theme(legend.text=element_text(size=50)) +
  theme(axis.title = element_text(size=50)) +
  theme(axis.text = element_text(size=40))
FigS7C$table <- FigS7C$table + theme(axis.text.y.left = element_text(size=40))

# save plot
png(filename="Figures/FigS7C.png", height=16, width=18, units="in", res=150)
FigS7C
dev.off()

cairo_pdf(file="Figures/FigS7C.pdf", height=16, width=18)
FigS7C
dev.off()
