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

################################################################################
################# pair down data to required cols ##############################

# clinical data
keep <- c('Subject ID', 'Study Arm', 'PFS Registration', 'PFS Status (0=alive; 1=progression or death)', 'First RECIST Assessment')
ClinicalDemographics <- ClinicalDemographics[,..keep]
setnames(ClinicalDemographics, 'First RECIST Assessment', 'First_RECIST_Assessment')

################################################################################
############### format data for analysis #######################################

# convert OS from supp. to days
daysInMonth <- 30.4167
ClinicalDemographics[,`PFS Registration` := `PFS Registration`/daysInMonth]

# factor
factorOrder <- c("PR", "SD", "PD")
ClinicalDemographics <- ClinicalDemographics[First_RECIST_Assessment %in% factorOrder]
ClinicalDemographics[,First_RECIST_Assessment := factor(First_RECIST_Assessment, levels=factorOrder)]

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
fit <- survfit(Surv(`PFS Registration`, `PFS Status (0=alive; 1=progression or death)`) ~ First_RECIST_Assessment, data = ClinicalDemographics)

# calculate hazard ratio
cox <- coxph(Surv(`PFS Registration`, `PFS Status (0=alive; 1=progression or death)`) ~ First_RECIST_Assessment, data = ClinicalDemographics)
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
FigS10B <- ggsurvplot(
  fit,
  data = ClinicalDemographics,
  xlab = "Time in Months",
  legend=c(.75, .75),
  xlim=c(0,36),
  ylab = "Progression-Free Survival",
  break.time.by = 6,
  size = 3,                 # change line size
  palette = c("#87c6c3", "#e0af46", "#5e292f"),# custom color palettes
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
  risk.table.y.text=F
)

# annotate Hazzard ratio
FigS10B$plot <- FigS10B$plot+
  ggplot2::annotate("text",
                    x = 7, y = 0.3, # x and y coordinates of the text
                    label = p_value, size = 15) +
  theme(legend.text=element_text(size=50)) +
  theme(axis.title = element_text(size=50)) +
  theme(axis.text = element_text(size=40))
FigS10B$table <- FigS10B$table + theme(axis.text.y.left = element_text(size=40))

# save plot
png(filename="Figures/FigS10B.png", height=16, width=18, units="in", res=150)
FigS10B
dev.off()

cairo_pdf(file="Figures/FigS10B.pdf", height=16, width=18)
FigS10B
dev.off()
