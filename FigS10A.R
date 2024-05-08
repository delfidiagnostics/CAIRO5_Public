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
ClinicalDemographics <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=1, skip=1))
Slopes <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=8, skip=1))

################################################################################
################# pair down data to required cols ##############################

# clinical data
keep <- c('Subject ID', 'Study Arm', 'PFS Registration', 'PFS Status (0=alive; 1=progression or death)', 'Durable Clinical Benefit')
ClinicalDemographics <- ClinicalDemographics[,..keep]

# merge data together
finalData <- merge(Slopes, ClinicalDemographics, by='Subject ID', all.x=T)

################################################################################
############### format data for analysis #######################################

# restrict down to eligible subjects, subjects 11 and 23 are ineligible
finalData <- finalData[!is.na(Velocity)]
finalData <- finalData[!`Subject ID` %in% c(11, 123)]

# calculate the cutpoint for groups
slopeQuant <- quantile(finalData$Velocity)
finalData[,slopeThreshold := ifelse(Velocity > slopeQuant[3], "Above\nDELFI-TF Slope", "Below\nDELFI-TF Slope")]

# restrict to those samples with a DCB
finalData <- finalData[`Durable Clinical Benefit` == 'DCB']

# convert OS from supp. to days
daysInMonth <- 30.4167
finalData[,`PFS Registration` := `PFS Registration`/daysInMonth]

# factor
factorOrder <- c("Below\nDELFI-TF Slope", "Above\nDELFI-TF Slope")
finalData[,slopeThreshold := factor(slopeThreshold, levels=factorOrder)]

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
fit <- survfit(Surv(`PFS Registration`, `PFS Status (0=alive; 1=progression or death)`) ~ slopeThreshold, data = finalData)

# calculate hazard ratio
cox <-coxph(Surv(`PFS Registration`, `PFS Status (0=alive; 1=progression or death)`) ~ slopeThreshold, data = finalData)
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
  theme(panel.grid=element_blank()) +
  theme(axis.text=element_text(size=28), axis.title=element_text(size=34)) +
  theme(legend.text=element_text(size=34), legend.title = element_blank())

# make the plot
FigS10A <- ggsurvplot(
  fit,
  data = finalData,
  xlab = "Time in Months",
  ylab = "Progression-Free Survival",
  legend=c(.73, .85),
  break.time.by = 6,
  xlim=c(0, 36),
  size = 3,                 # change line size
  palette = c("#8da8bf", "#da7b3d"),# custom color palettes
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs = c("Below\nDELFI-TF Slope:DCB", "Above\nDELFI-TF Slope:DCB"),
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
FigS10A$plot <- FigS10A$plot+
  ggplot2::annotate("text",
                    x = 4, y = 0.25, # x and y coordinates of the text
                    label = p_value, size = 15) +
  theme(legend.text=element_text(size=50)) +
  theme(axis.title = element_text(size=50)) +
  theme(axis.text = element_text(size=40))
FigS10A$table <- FigS10A$table + theme(axis.text.y.left = element_text(size=40))

# save plot
png(filename="Figures/FigS10A.png", height=16, width=18, units="in", res=150)
FigS10A
dev.off()

cairo_pdf(file="Figures/FigS10A.pdf", height=16, width=18)
FigS10A
dev.off()
