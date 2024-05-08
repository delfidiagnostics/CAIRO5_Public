################################################################################
#################### Setup #####################################################

# loc
setwd("~/git/CAIRO5_Public")

# libraries
library(readxl)
library(data.table)
library(survival)
library(survminer)

# Data
ClinicalDemographics <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=1, skip=1))
Slopes <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=8, skip=1))

################################################################################
################# pair down data to required cols ##############################

# clinical data
keep <- c('Subject ID', 'OS Registration', 'OS Status (0=alive; 1=death)')
ClinicalDemographics <- ClinicalDemographics[,..keep]

# merge data together
finalData <- merge(ClinicalDemographics, Slopes, by='Subject ID', all.x=T)

################################################################################
################ final formatting for survival #################################

# restrict down to eligible subjects, subjects 11 and 23 are ineligible
finalData <- finalData[!is.na(Velocity)]
finalData <- finalData[!`Subject ID` %in% c(11, 123)]

# annotate slope direction
medianVelocity <- quantile(finalData$Velocity)[3]
finalData[,slopeThreshold := ifelse(Velocity > medianVelocity, "Above\nDELFI-TF Slope", "Below\nDELFI-TF Slope")]

# factor
factorOrder <- c("Below\nDELFI-TF Slope", "Above\nDELFI-TF Slope")
finalData[,slopeThreshold := factor(slopeThreshold, levels=factorOrder)]

# convert PFS from days to months
daysInMonth <- 30.4167
finalData[,`OS Registration` := `OS Registration`/daysInMonth]

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
fit <- survfit(Surv(`OS Registration`, `OS Status (0=alive; 1=death)`) ~ slopeThreshold, data = finalData)

# calculate hazard ratio
cox <-coxph(Surv(`OS Registration`, `OS Status (0=alive; 1=death)`) ~ slopeThreshold, data = finalData)
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
Fig4C <- ggsurvplot(
  fit,
  data = finalData,
  xlab = "Time in Months",
  xlim=c(0,60),
  legend=c(.75,.8),
  ylab = "Overall Survival",
  break.time.by = 6,
  size = 3,                 # change line size
  palette = c("#8da8bf", "#da7b3d"),# custom color palettes
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

# annotate Hazzard ratio and adjust plot
Fig4C$plot <- Fig4C$plot+
  ggplot2::annotate("text",
                    x = 6, y = 0.25, # x and y coordinates of the text
                    label = p_value, size = 15) +
  theme(legend.text=element_text(size=50)) +
  theme(axis.title = element_text(size=50)) +
  theme(axis.text = element_text(size=40))
Fig4C$table <- Fig4C$table + theme(axis.text.y.left = element_text(size=40))


png(filename="Figures/Fig4C.png", height=16, width=18, units="in", res=150)
Fig4C
dev.off()

cairo_pdf(file="Figures/Fig4C.pdf", height=16, width=18)
Fig4C
dev.off()
