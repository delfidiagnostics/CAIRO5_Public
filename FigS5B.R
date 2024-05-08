################################################################################
#################### Setup #####################################################

# loc
setwd("~/git/CAIRO5_Public")

# libraries
library(readxl)
library(data.table)
library(ggplot2)
library(ggallin)
library(scales)

# Data
BloodSampleCharacteristics <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=2, skip=1))
SummaryCopyNumber <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=5, skip=1))

################################################################################
################# pair down data to required cols ##############################

# Blood Samples
keep <- c('Subject ID', 'Sample Time Point', 'DELFI-TF', 'MAF', 'ichorCNA')
BloodSampleCharacteristics <- BloodSampleCharacteristics[,..keep]
setkeyv(BloodSampleCharacteristics, c('Subject ID', 'Sample Time Point'))

# CNA data
setkeyv(SummaryCopyNumber, c('Subject ID', 'Sample Time Point'))

# merge data together
finalData <- merge(SummaryCopyNumber, BloodSampleCharacteristics)

################################################################################
################ format data ###################################################

# melt the data first by technology and then by CNA feature
finalData <- melt(finalData, measure.vars=c('DELFI-TF', 'MAF'), variable.name="Technology", value.name="Value")
finalData <- melt(finalData, measure.vars=c('PLCG1 Copy Number', 'MBD1 Copy Number'), variable.name="Gene", value.name="CN")

# rename copy number value names
finalData[Gene == 'PLCG1 Copy Number', Gene := 'PLCG1 (Chr20q)']
finalData[Gene == 'MBD1 Copy Number', Gene := 'MBD1 (Chr18q)']

# rename MAF for plot
finalData[Technology == 'MAF', Technology := 'RAS/BRAF MAF']

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

# run correlations
MAF_PLCG1 <- cor.test(finalData[Technology == 'RAS/BRAF MAF' & Gene == 'PLCG1 (Chr20q)']$Value,
                      finalData[Technology == 'RAS/BRAF MAF' & Gene == 'PLCG1 (Chr20q)']$CN)
MAF_PLCG1 <- paste0("r=", signif(MAF_PLCG1$estimate, 2), "\n", relabel_p(MAF_PLCG1$p.value))

delfi_PLCG1 <- cor.test(finalData[Technology == 'DELFI-TF' & Gene == 'PLCG1 (Chr20q)']$Value,
                        finalData[Technology == 'DELFI-TF' & Gene == 'PLCG1 (Chr20q)']$CN)
delfi_PLCG1 <- paste0("r=", signif(delfi_PLCG1$estimate, 2), "\n", relabel_p(delfi_PLCG1$p.value))

MAF_mbd1 <- cor.test(finalData[Technology == 'RAS/BRAF MAF' & Gene == 'MBD1 (Chr18q)']$Value,
                     finalData[Technology == 'RAS/BRAF MAF' & Gene == 'MBD1 (Chr18q)']$CN)
MAF_mbd1 <- paste0("r=", signif(MAF_mbd1$estimate, 2), "\n", relabel_p(MAF_mbd1$p.value))

delfi_mbd1 <- cor.test(finalData[Technology == 'DELFI-TF' & Gene == 'MBD1 (Chr18q)']$Value,
                       finalData[Technology == 'DELFI-TF' & Gene == 'MBD1 (Chr18q)']$CN)
delfi_mbd1 <- paste0("r=", signif(delfi_mbd1$estimate, 2), "\n", relabel_p(delfi_mbd1$p.value))


# make statistics DT
statDT_mbd1 <- data.table(Technology=c('DELFI-TF', 'RAS/BRAF MAF'),
                          Gene=c('MBD1 (Chr18q)'),
                          x=c(.12, .15),
                          y=c(-1),
                          label=c(delfi_mbd1, MAF_mbd1))

statDT_PLCG1 <- data.table(Technology=c('DELFI-TF', 'RAS/BRAF MAF'),
                           Gene=c('PLCG1 (Chr20q)'),
                           x=c(.12, .15),
                           y=c(.85),
                           label=c(delfi_PLCG1, MAF_PLCG1))


statDT <- rbindlist(list(statDT_mbd1, statDT_PLCG1))

################################################################################
############## Make final plots ################################################

# construct the plot
FigS5B <- ggplot() +
  geom_point(data=finalData, mapping=aes(x=Value, y=CN, color=Technology), alpha=.5, size=3) +
  geom_smooth(data=finalData, mapping=aes(x=Value, y=CN, color=Technology), method='glm') +
  scale_color_manual("",values=c("DELFI-TF"="#046C9A", "RAS/BRAF MAF"="#D69C4E")) +
  scale_x_continuous(labels=scales::percent) +
  geom_text(data=statDT, mapping=aes(x=x, y=y, label=label), size=12) +
  facet_grid(Gene~Technology, scales="free") +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  ylab("cfDNA Log2 Copy Number Ratio") +
  xlab("Tumor Fraction") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(strip.text=element_text(color="white", size=42)) +
  theme(strip.background = element_rect(fill="#1B465A")) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=36), axis.title=element_text(size=42)) +
  theme(legend.text=element_text(size=36))

# save the results
png(filename="Figures/FigS5B.png", height=14, width=16, units="in", res=150)
FigS5B
dev.off()

cairo_pdf(file="Figures/FigS5B.pdf", height=14, width=16)
FigS5B
dev.off()
