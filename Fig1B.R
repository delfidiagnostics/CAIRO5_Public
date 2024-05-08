################################################################################
#################### Setup #####################################################

# loc
setwd("~/git/CAIRO5_Public")

# libraries
library(readxl)
library(data.table)
library(ggplot2)
library(scales)
library(gtools)
library(ggallin)
library(viridis)
library(grid)
library(gridExtra)

# Data
BloodSampleCharacteristics <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=2, skip=1))
FragmentScores <- readRDS('Data/SimulatedFragmentRatios.rds')

################################################################################
################# pair down data to required cols ##############################

# Blood Samples
keep <- c('Sequence ID', 'Subject ID', 'Sample Time Point', 'DELFI-TF')
BloodSampleCharacteristics <- BloodSampleCharacteristics[,..keep]

################################################################################
############## rescale and combine the model data ##############################

# format and rescale Fragment Scores
FragmentScores[,`Chromosome Arm` := NULL]
setnames(FragmentScores, 'Simulated Fragment Ratio', 'Value')
FragmentScores[,Variable := 'S/L Ratio']
FragmentScores[,Bin := paste0('f_', Bin)]

FragmentScores[,Value := rescale(Value, to=c(-1, 1))]
FragmentScores[,Value := Value - median(Value)]

# merge all data together
finalData <- merge(FragmentScores, BloodSampleCharacteristics, by.x='id', by.y='Sequence ID', all.x=T)

# set order for samples based on DELFI-TF
delfiScoreOrder <- rev(BloodSampleCharacteristics[order(`DELFI-TF`)]$`Sequence ID`)
finalData[,id := factor(id, levels=delfiScoreOrder)]

################################################################################
################## Cluster Bins ################################################

# convert to wide format
heatmapDataWide <- dcast(finalData, Bin~id, value.var = "Value")
rowNames <- heatmapDataWide$Bin
heatmapDataWide <- as.matrix(heatmapDataWide[,2:ncol(heatmapDataWide)])
rownames(heatmapDataWide) <- rowNames

##############################################################################
################ Perform Hiearchical clustering ##############################

# compute distance matrix
heatmapDistanceBins <- dist(heatmapDataWide)

# perform clustering
heatmapClusterBins <- hclust(heatmapDistanceBins, method="ward.D")

# bin dendogram
heatmapBinModel <- as.dendrogram(heatmapClusterBins)

# relevel heatmap data based on clusters and make plot, note we relevel samples by delfi-tf value above
finalData[,Bin := factor(Bin, levels=rev(labels(heatmapBinModel)))]

################################################################################
############ make plots ########################################################

# main heatmap plot
heatmap <- ggplot(finalData, aes(x=id, y=Bin, fill=Value)) +
  geom_raster() +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_discrete(expand=c(0, 0)) +
  scale_fill_gradientn('Feature\nValues', trans=ssqrt_trans, colors=c("blue", "white", "red")) +
  theme(axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank()) +
  theme(legend.title=element_text(size=96, face="bold")) +
  theme(legend.text=element_text(size=88, hjust=1)) +
  theme(legend.key.size = unit(1, 'cm')) +
  theme(legend.key.height = unit(18, 'cm')) +
  theme(legend.key.width = unit(3, 'cm'))

# delfi-tf sub-plot
clinicalContinuous <- ggplot(unique(finalData[,.(id, `DELFI-TF`)]), aes(x=id, y=1, fill=`DELFI-TF`)) +
  geom_raster() +
  scale_fill_viridis(option="G", trans=ssqrt_trans, direction=-1) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_discrete(expand=c(0, 0)) +
  guides(fill=guide_colourbar(reverse = T)) +
  theme(axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank()) +
  theme(plot.margin = margin(0, 0, 0, 0, "cm")) +
  theme(legend.position = "none")

# extract legend for sub-plot
clinicalDELFITFLegend <- ggplot(unique(finalData[,.(id, `DELFI-TF`)]), aes(x=id, y=1, fill=`DELFI-TF`)) +
  geom_raster() +
  scale_fill_viridis(option="G", "DELFI-TF", trans=ssqrt_trans, direction=-1, breaks=c(.6, .4, .2, .01)) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_discrete(expand=c(0, 0)) +
  guides(fill=guide_colourbar(reverse = T)) +
  theme(axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank()) +
  theme(plot.margin = margin(0, 0, 0, 0, "cm")) +
  theme(legend.position = "top") +
  theme(legend.title=element_text(size=96, face="bold")) +
  theme(legend.text=element_text(size=88, hjust=1)) +
  theme(legend.key.width = unit(16, 'cm')) +
  theme(legend.key.height = unit(2, 'cm')) +
  theme(legend.justification=c(0,0)) +
  theme(legend.margin = margin(0, 0, 0, 0, "cm"))
clinicalDELFITFLegend <- ggplotGrob(clinicalDELFITFLegend)$grobs[[18]]

################################################################################
############### combine all plots together #####################################

# convert everything to grobs
blankPanel <- grid.rect(gp=gpar(col="white"))
heatmapGrob <- ggplotGrob(heatmap)
clinicalContinuousGrob <- ggplotGrob(clinicalContinuous)

# make sure every width between all grobs is the same
maxWidth <- unit.pmax(heatmapGrob$widths, clinicalContinuousGrob$widths)
heatmapGrob$widths <- as.list(maxWidth)
clinicalContinuousGrob$widths <- as.list(maxWidth)

# arrange plot
Fig1B <- grid.arrange(clinicalContinuousGrob, heatmapGrob, blankPanel, clinicalDELFITFLegend, ncol=1, heights=c(.025, .95, .01, .025))

# save plot
png(filename="Figures/Fig1B.png", height=60, width=120, units="in", res=150)
grid::grid.draw(Fig1B)
dev.off()

cairo_pdf(file="Figures/Fig1B.pdf", height=60, width=120)
grid::grid.draw(Fig1B)
dev.off()

