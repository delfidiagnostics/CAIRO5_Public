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
keep <- c('Sequence ID', 'Subject ID', 'Sample Time Point', 'DELFI-TF', 'MAF')
BloodSampleCharacteristics <- BloodSampleCharacteristics[,..keep]

# clinical data
keep <- c('Subject ID', 'Study Arm', 'OS Registration', 'OS Status (0=alive; 1=death)')
ClinicalDemographics <- ClinicalDemographics[,..keep]

# merge data together
finalData <- merge(BloodSampleCharacteristics, ClinicalDemographics, by='Subject ID', all.x=T)

################################################################################
############### format data for analysis #######################################

# restrict to ARM 1 samples at baseline
finalData <- finalData[`Study Arm` == 1 & `Sample Time Point` == 'T0']

# convert OS from supp. to days
daysInMonth <- 30.4167
finalData[,`OS Registration` := `OS Registration`/daysInMonth]

# classify high/low MAF
finalData[,mafStatus := as.character(NA)]
finalData[MAF <= quantile(finalData$MAF)[2], "mafStatus"] <- "Low MAF"
finalData[MAF > quantile(finalData$MAF)[2], "mafStatus"] <- "High MAF"

# classify high/low delfi-tf
finalData[,dmsStatus := as.character(NA)]
finalData[`DELFI-TF` <= quantile(finalData$`DELFI-TF`)[2], "dmsStatus"] <- "Low DELFI-TF"
finalData[`DELFI-TF` > quantile(finalData$`DELFI-TF`)[2], "dmsStatus"] <- "High DELFI-TF"

# melt
finalData <- melt(finalData, id.vars=c('Subject ID', 'OS Registration', 'OS Status (0=alive; 1=death)'), measure.vars = c('mafStatus', 'dmsStatus'))
factorOrder <- c("Low MAF", "Low DELFI-TF", "High MAF", "High DELFI-TF")
finalData[,value := factor(value, levels=factorOrder)]

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
fit <- survfit(Surv(`OS Registration`, `OS Status (0=alive; 1=death)`) ~ value, data = finalData)

# calculate hazard ratio
cox <-coxph(Surv(`OS Registration`, `OS Status (0=alive; 1=death)`) ~ value, data = finalData)
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
# appears to be harmless warning with new ggplot: https://github.com/kassambara/survminer/issues/643
Fig3D <- ggsurvplot(
  fit,
  data = finalData,
  xlab = "Time in Months",
  legend=c(.75, .75),
  xlim=c(0,60),
  ylab = "Overall Survival",
  break.time.by = 6,
  size = 3,                 # change line size
  palette = c("#efa926", '#6e7e40', "#25b1c0",  '#7b0f23'),# custom color palettes
  legend.labs=factorOrder,
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
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
Fig3D$plot <- Fig3D$plot+
  ggplot2::annotate("text",
                    x = 4, y = 0.25, # x and y coordinates of the text
                    label =p_value, size = 15) +
  theme(legend.text=element_text(size=50)) +
  theme(axis.title = element_text(size=50)) +
  theme(axis.text = element_text(size=40))
Fig3D$table <- Fig3D$table + theme(axis.text.y.left = element_text(size=40))

# save plot
png(filename="Figures/Fig3D.png", height=16, width=18, units="in", res=150)
Fig3D
dev.off()

cairo_pdf(file="Figures/Fig3D.pdf", height=16, width=18)
Fig3D
dev.off()
