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
keep <- c('Subject ID', 'Study Arm', 'OS Registration', 'OS Status (0=alive; 1=death)', 'Complete Resection (1=Yes; 0=No)')
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
finalData[,slopeThreshold := ifelse(Velocity > slopeQuant[3], "Above Median", "Below Median")]

# rename CR to make more sense
finalData[,crStatus := ifelse(as.numeric(`Complete Resection (1=Yes; 0=No)`) == 1, "CR", "NCR")]

# if CR is NA there was no surgery, classify this as NCR
finalData[is.na(crStatus), crStatus := 'NCR']

# convert OS from supp. to days
daysInMonth <- 30.4167
finalData[,`OS Registration` := `OS Registration`/daysInMonth]

# factor
factorOrder <- c("Below Median", "Above Median")
finalData[,slopeThreshold := factor(slopeThreshold, levels=factorOrder)]

factorOrder <- c("CR", "NCR")
finalData[,crStatus := factor(crStatus, levels=factorOrder)]

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
fit <- survfit(Surv(`OS Registration`, `OS Status (0=alive; 1=death)`) ~ slopeThreshold + crStatus, data = finalData)

# calculate hazard ratio
cox <-coxph(Surv(`OS Registration`, `OS Status (0=alive; 1=death)`) ~ slopeThreshold + crStatus, data = finalData)
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
  theme(legend.text=element_text(size=30), legend.title = element_blank())

# make the plot
FigS10C <- ggsurvplot(
  fit,
  data = finalData,
  xlab = "Time in Months",
  ylab = "Overall Survival",
  legend=c(.18, .25),
  break.time.by = 6,
  xlim=c(0, 60),
  size = 3,                 # change line size
  palette = c("#8b4513", "#556b2f", "#ff8c00", "#cd5c5c"),# custom color palettes
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs = c("Below\nDELFI-TF Slope:CR", "Below\nDELFI-TF Slope:NCR", "Above\nDELFI-TF Slope:CR", "Above\nDELFI-TF Slope:NCR"),
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
FigS10C$plot <- FigS10C$plot+
  ggplot2::annotate("text",
                    x = 54, y = 0.85, # x and y coordinates of the text
                    label = p_value, size = 15) +
  theme(legend.text=element_text(size=30)) +
  theme(axis.title = element_text(size=50)) +
  theme(axis.text = element_text(size=40))
FigS10C$table <- FigS10C$table + theme(axis.text.y.left = element_text(size=30))


# save plot
png(filename="Figures/FigS10C.png", height=16, width=18, units="in", res=150)
FigS10C
dev.off()

cairo_pdf(file="Figures/FigS10C.pdf", height=16, width=18)
FigS10C
dev.off()




