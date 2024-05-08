################################################################################
#################### Setup #####################################################

# loc
setwd("~/git/CAIRO5_Public")

# libraries
library(readxl)
library(data.table)
library(ggplot2)
library(ggallin)

# Data
LungCancerValidation <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=6, skip=1))

################################################################################
###########  run correlation ###################################################

# switch statement to re-define p-values
relabel_p <- function(x){
  if(x < .0001) return('p < 0.0001')
  if(x < .001) return('p < 0.001')
  if(x < .01) return('p < 0.01')
  if(x < .05) return('p < 0.05')
  return(paste0('p = ', signif(x, 3)))
}

# run correlation
corStat <- cor.test(LungCancerValidation$`DELFI-TF`, LungCancerValidation$`Maximum MAF`)
corStat <- paste0("r=", signif(corStat$estimate, 2), "\n", relabel_p(corStat$p.value))

################################################################################
############ make plot #########################################################

FigS5C <- ggplot(LungCancerValidation, aes(x=`Maximum MAF`, y=`DELFI-TF`)) +
  geom_point(alpha=.5, color='firebrick', size=4) +
  geom_smooth(method="glm", color='grey20') +
  annotate("text", x=.015, y=.4, label=corStat, size=14) +
  scale_y_continuous(trans=ssqrt_trans, labels=scales::percent, breaks=c(0, .025, .1, .2, .4, .6), limits=c(0, .6)) +
  scale_x_continuous(trans=ssqrt_trans, labels=scales::percent, breaks=c(0, .025, .1, .2, .4, .6), limits=c(0, .6)) +
  ylab("DELFI-TF") +
  xlab("Maximum Observed MAF") +
  facet_wrap(~"Validation Cohort") +
  theme_bw() +
  theme(axis.text.x=element_text(size=36), axis.text.y=element_text(size=36)) +
  theme(axis.title=element_text(size=42)) +
  theme(legend.position = "none") +
  theme(panel.grid=element_blank()) +
  theme(strip.background=element_rect(fill="#1B465A")) +
  theme(strip.text=element_text(color="white", size=42, face="bold"))

# save output
png(filename="Figures/FigS5C.png", height=14, width=16, units="in", res=150)
FigS5C
dev.off()

cairo_pdf(file="Figures/FigS5C.pdf", height=14, width=16)
FigS5C
dev.off()