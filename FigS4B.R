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
library(scales)
library(ggallin)
library(rstatix)

# Data
BloodSampleCharacteristics <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=2, skip=1))
ClinicalDemographics <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=1, skip=1))

################################################################################
################# pair down data to required cols ##############################

# Blood Samples
keep <- c('Sequence ID', 'Subject ID', 'Sample Time Point', 'DELFI-TF', 'MAF', 'ichorCNA')
BloodSampleCharacteristics <- BloodSampleCharacteristics[,..keep]

# clinical data
keep <- c('Subject ID', 'Study Arm')
ClinicalDemographics <- ClinicalDemographics[,..keep]

# merge data together
finalData <- merge(BloodSampleCharacteristics, ClinicalDemographics, by='Subject ID', all.x=T)

################################################################################
################### final formatting ###########################################

# split out MT and WT arms for formatting and clean up WT where we have NA for MAF
finalBoxMT <- finalData[`Study Arm` == 1]

finalBoxWT <- finalData[`Study Arm` == 3]
finalBoxWT[,MAF := NULL]

# restrict to only values with MAF, DELFI-TF, and ichorCNA and remove undetectable MAF
finalBoxMT <- finalBoxMT[(!is.na(ichorCNA) & !is.na(`DELFI-TF`) & !is.na(MAF))]
finalBoxMT <- finalBoxMT[MAF != 0]

# calculate cut points for x-axis categories
mafQuant <- quantile(finalBoxMT$MAF)

# create MAF categories
finalBoxMT[MAF > 0 & MAF <= mafQuant[2], `:=`(mafCategory = 'Group 1', Label = 'MT Arm')]
finalBoxMT[MAF > mafQuant[2] & MAF <= mafQuant[3], `:=`( mafCategory = 'Group 2', Label='MT Arm')] 
finalBoxMT[MAF > mafQuant[3] & MAF <= mafQuant[4], `:=`(mafCategory = 'Group 3', Label='MT Arm')] 
finalBoxMT[MAF > mafQuant[4], `:=`(mafCategory='Group 4', Label='MT Arm')]

# define DELFI-TF categories based on MAF categories from Arm 1
finalBoxWT <- finalBoxWT[!is.na(`DELFI-TF`) & !is.na(ichorCNA)]
finalBoxWT[`DELFI-TF` <= mafQuant[2], `:=`(mafCategory = 'Group 1', Label='WT Arm')]
finalBoxWT[`DELFI-TF` > mafQuant[2] & `DELFI-TF` <= mafQuant[3], `:=`(mafCategory = 'Group 2', Label='WT Arm')]
finalBoxWT[`DELFI-TF` > mafQuant[3] & `DELFI-TF` <= mafQuant[4], `:=`(mafCategory = 'Group 3', Label='WT Arm')]
finalBoxWT[`DELFI-TF` > mafQuant[4], `:=`(mafCategory = 'Group 4', Label='WT Arm')]

# melt the data
finalBoxMT <- melt(finalBoxMT, measure.vars = c('DELFI-TF', 'ichorCNA', 'MAF'))
finalBoxWT <- melt(finalBoxWT, measure.vars=c('DELFI-TF', 'ichorCNA'))

# bind everything together and remove unassigned values
finalBox <- rbindlist(list(finalBoxMT, finalBoxWT))

# factor everything for plotting
finalBox[,mafCategory := factor(mafCategory, levels=c('Group 1', 'Group 2', 'Group 3', 'Group 4'))]
finalBox[,Label := factor(Label, levels=c("MT Arm", "WT Arm"))]
finalBox[,variable := factor(variable, levels=c('MAF', 'DELFI-TF', 'ichorCNA'))]
finalBox[,`Sequence ID` := factor(`Sequence ID`)]

################################################################################
############# Run Statistics on MT Arm to annotate #############################

# we need separate DT for each MAF category so the paired test works correctly
test1 <- finalBox[mafCategory=='Group 1' & Label == 'MT Arm']
test1[,`Sequence ID` := factor(`Sequence ID`)]
test1[,variable := factor(variable)]
test1 <- test1[order(`Sequence ID`, variable)]
test1 <- wilcox_test(value ~ variable, data=test1, paired=T, p.adjust.method = "bonferroni")

test2 <- finalBox[mafCategory=='Group 2' & Label == 'MT Arm']
test2[,`Sequence ID` := factor(`Sequence ID`)]
test2[,variable := factor(variable)]
test2 <- test2[order(`Sequence ID`, variable)]
test2 <- wilcox_test(value ~ variable, data=test2, paired=T, p.adjust.method = "bonferroni")

test3 <- finalBox[mafCategory=='Group 3' & Label == 'MT Arm']
test3[,`Sequence ID` := factor(`Sequence ID`)]
test3[,variable := factor(variable)]
test3 <- test3[order(`Sequence ID`, variable)]
test3 <- wilcox_test(value ~ variable, data=test3, paired=T, p.adjust.method = "bonferroni")

test4 <- finalBox[mafCategory=='Group 4' & Label == 'MT Arm']
test4[,`Sequence ID` := factor(`Sequence ID`)]
test4[,variable := factor(variable)]
test4 <- test4[order(`Sequence ID`, variable)]
test4 <- wilcox_test(value ~ variable, data=test4, paired=T, p.adjust.method = "bonferroni")

# make statistic DT
statDTMT_test1 <- data.table(mafCategory=c('Group 1'),
                             pvalue=c(test1$p.adj[1], test1$p.adj[2], test1$p.adj[3]),
                             xpos=c(.875, 1, 1.125),
                             ypos=c(1, 1.1, 1.2))
statDTMT_test2 <- data.table(mafCategory=c('Group 2'),
                             pvalue=c(test2$p.adj[1], test2$p.adj[2], test2$p.adj[3]),
                             xpos=c(1.875, 2, 2.125),
                             ypos=c(1, 1.1, 1.2))
statDTMT_test3 <- data.table(mafCategory=c('Group 3'),
                             pvalue=c(test3$p.adj[1], test3$p.adj[2], test3$p.adj[3]),
                             xpos=c(2.875, 3, 3.125),
                             ypos=c(1, 1.1, 1.2))
statDTMT_test4 <- data.table(mafCategory=c('Group 4'),
                             pvalue=c(test4$p.adj[1], test4$p.adj[2], test4$p.adj[3]),
                             xpos=c(3.875, 4, 4.125),
                             ypos=c(1, 1.1, 1.2))
statDTMT <- rbindlist(list(statDTMT_test1, statDTMT_test2, statDTMT_test3, statDTMT_test4))
statDTMT[,mafCategory := factor(mafCategory, levels=c('Group 1', 'Group 2', 'Group 3', 'Group 4'))]

################################################################################
############# Run Statistics on WT Arm to annotate #############################

# we need separate DT for each MAF category so the paired test works correctly
test1 <- finalBox[mafCategory=='Group 1' & Label == 'WT Arm']
test1[,`Sequence ID` := factor(`Sequence ID`)]
test1[,variable := factor(variable)]
test1 <- test1[order(`Sequence ID`, variable)]
test1 <- wilcox_test(value ~ variable, data=test1, paired=T, p.adjust.method = "bonferroni")

test2 <- finalBox[mafCategory=='Group 2' & Label == 'WT Arm']
test2[,`Sequence ID` := factor(`Sequence ID`)]
test2[,variable := factor(variable)]
test2 <- test2[order(`Sequence ID`, variable)]
test2 <- wilcox_test(value ~ variable, data=test2, paired=T, p.adjust.method = "bonferroni")

test3 <- finalBox[mafCategory=='Group 3' & Label == 'WT Arm']
test3[,`Sequence ID` := factor(`Sequence ID`)]
test3[,variable := factor(variable)]
test3 <- test3[order(`Sequence ID`, variable)]
test3 <- wilcox_test(value ~ variable, data=test3, paired=T, p.adjust.method = "bonferroni")

test4 <- finalBox[mafCategory=='Group 4' & Label == 'WT Arm']
test4[,`Sequence ID` := factor(`Sequence ID`)]
test4[,variable := factor(variable)]
test4 <- test4[order(`Sequence ID`, variable)]
test4 <- wilcox_test(value ~ variable, data=test4, paired=T, p.adjust.method = "bonferroni")

# make statistic DT
statDTWT <- data.table(mafCategory=c('Group 1', 'Group 2', 'Group 3', 'Group 4'),
                       pvalue=c(test1$p, test2$p, test3$p, test4$p),
                       ypos=c(1.1))
statDTWT[,mafCategory := factor(mafCategory, levels=c('Group 1', 'Group 2', 'Group 3', 'Group 4'))]

################################################################################
############### final formatting ###############################################

# switch statement to re-define p-values
relabel_p <- function(x){
  if(x < .0001) return('****')
  if(x < .001) return('***')
  if(x < .01) return('**')
  if(x < .05) return('*')
  return('ns')
}
statDTMT[,pvalue := sapply(statDTMT$pvalue, relabel_p)]
statDTWT[,pvalue := sapply(statDTWT$pvalue, relabel_p)]

################################################################################
############## construct plot ##################################################

# mt arm plot
p1 <- ggplot() +
  geom_hline(yintercept = c(mafQuant[2], mafQuant[3], mafQuant[4]), linetype=4, color='grey30', alpha=.5, linewidth=1.5) +
  geom_boxplot(data=finalBox[Label == 'MT Arm'], mapping=aes(x=mafCategory, y=value, fill=variable), outlier.shape=NA, position=position_dodge2(preserve = "single")) +
  geom_point(data=finalBox[Label == 'MT Arm'], mapping=aes(x=mafCategory, y=value, fill=variable), position=position_jitterdodge(jitter.width=.2), alpha=.65, show.legend = FALSE, size=2.5) +
  geom_text(data=statDTMT, mapping=aes(x=xpos, y=ypos, label=pvalue), size=10) +
  annotate("segment", x=0.75, xend=1, y=.95, yend=.95, linewidth=1.2, lineend="round", color="grey20") +
  annotate("segment", x=0.75, xend=1.25, y=1.05, yend=1.05, linewidth=1.2, lineend="round", color="grey20") +
  annotate("segment", x=1, xend=1.25, y=1.15, yend=1.15, linewidth=1.2, lineend="round", color="grey20") +
  annotate("segment", x=1.75, xend=2, y=.95, yend=.95, linewidth=1.2, lineend="round", color="grey20") +
  annotate("segment", x=1.75, xend=2.25, y=1.05, yend=1.05, linewidth=1.2, lineend="round", color="grey20") +
  annotate("segment", x=2, xend=2.25, y=1.15, yend=1.15, linewidth=1.2, lineend="round", color="grey20") +
  annotate("segment", x=2.75, xend=3, y=.95, yend=.95, linewidth=1.2, lineend="round", color="grey20") +
  annotate("segment", x=2.75, xend=3.25, y=1.05, yend=1.05, linewidth=1.2, lineend="round", color="grey20") +
  annotate("segment", x=3, xend=3.25, y=1.15, yend=1.15, linewidth=1.2, lineend="round", color="grey20") +
  annotate("segment", x=3.75, xend=4, y=.95, yend=.95, linewidth=1.2, lineend="round", color="grey20") +
  annotate("segment", x=3.75, xend=4.25, y=1.05, yend=1.05, linewidth=1.2, lineend="round", color="grey20") +
  annotate("segment", x=4, xend=4.25, y=1.15, yend=1.15, linewidth=1.2, lineend="round", color="grey20") +
  scale_y_continuous(labels=scales::percent, trans=ssqrt_trans, limits=c(0, 1.21), breaks=c(0, .01, .04, .08, .15, .25, .5, .75)) +
  scale_fill_manual(values=c("DELFI-TF"="#046C9A", "MAF"="#D69C4E", "ichorCNA"="#7d988d")) +
  facet_wrap(~Label) +
  ylab('Tumor Fraction') +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(strip.text=element_text(color="white", size=42)) +
  theme(strip.background = element_rect(fill="#1B465A")) +
  theme(legend.position="top") +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x=element_text(size=32), axis.text.y=element_text(size=32), axis.title.y=element_text(size=38)) +
  theme(legend.text=element_text(size=32), legend.title=element_blank()) +
  theme(legend.key.size = unit(3,"line"))

# wt arm plot
p2 <- ggplot() +
  geom_hline(yintercept = c(mafQuant[2], mafQuant[3], mafQuant[4]), linetype=4, color='grey30', alpha=.5, linewidth=1.5) +
  geom_boxplot(data=finalBox[Label == 'WT Arm'], mapping=aes(x=mafCategory, y=value, fill=variable), outlier.shape=NA, position=position_dodge2(preserve = "single")) +
  geom_point(data=finalBox[Label == 'WT Arm'], mapping=aes(x=mafCategory, y=value, fill=variable), position=position_jitterdodge(jitter.width=.2), alpha=.65, show.legend = FALSE, size=2.5) +
  geom_text(data=statDTWT, mapping=aes(x=mafCategory, y=ypos, label=pvalue), size=10) +
  annotate("segment", x=.8, xend=1.2, y=1.05, yend=1.05, linewidth=1.2, lineend="round", color="grey20") +
  annotate("segment", x=1.8, xend=2.2, y=1.05, yend=1.05, linewidth=1.2, lineend="round", color="grey20") +
  annotate("segment", x=2.8, xend=3.2, y=1.05, yend=1.05, linewidth=1.2, lineend="round", color="grey20") +
  annotate("segment", x=3.8, xend=4.2, y=1.05, yend=1.05, linewidth=1.2, lineend="round", color="grey20") +
  scale_y_continuous(labels=scales::percent, trans=ssqrt_trans, limits=c(0, 1.21), sec.axis=dup_axis(name="MAF Quartiles", breaks=c(mafQuant[2], mafQuant[3], mafQuant[4]), labels=c(paste0('Q1:', round(mafQuant[2] * 100, 1), '%'),paste0('Q2:', round(mafQuant[3] * 100, 1), '%'),paste0('Q3:', round(mafQuant[4] * 100, 1), '%')))) +
  scale_fill_manual(values=c("DELFI-TF"="#046C9A", "MAF"="#D69C4E", "ichorCNA"="#7d988d")) +
  facet_wrap(~Label) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(strip.text=element_text(color="white", size=42)) +
  theme(strip.background = element_rect(fill="#1B465A")) +
  theme(legend.position="top") +
  theme(axis.text.x=element_text(size=32)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.y.right = element_text(size=32), axis.title.y.right=element_text(size=38)) +
  theme(axis.text.y.left = element_blank(), axis.title.y.left=element_blank(), axis.ticks.y.left = element_blank()) +
  theme(legend.text=element_text(size=32), legend.title=element_blank()) +
  theme(legend.key.size = unit(3,"line"))

# save plots
png(filename='Figures/FigS4B.png', height=20, width=32, units="in", res=150)
grid.draw(arrangeGrob(p1, p2, ncol=2))
dev.off()

cairo_pdf(file='Figures/FigS4B.pdf', height=20, width=32)
grid.draw(arrangeGrob(p1, p2, ncol=2))
dev.off()
