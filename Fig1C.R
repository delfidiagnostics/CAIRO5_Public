################################################################################
#################### Setup #####################################################

# loc
setwd("~/git/CAIRO5_Public")

# libraries
library(readxl)
library(data.table)
library(ggplot2)
library(gtools)
library(gridExtra)
library(grid)

# Data
BloodSampleCharacteristics <- as.data.table(read_xlsx("Data/CAIRO5_Supplemental.xlsx", sheet=2, skip=1))
FragmentRatios <- readRDS('Data/SimulatedFragmentRatios.rds')
recistCategory <- readRDS('Data/Fig1c_AuxData.rds')
binCoords <- readRDS('Data/hg19_bin_coordinates.rds')

################################################################################
############# pair down to required cols for fragmentation profile #############

# bin coords
binCoords[,range := NULL]
binCoords[,chr := NULL]

# perform merges
FragmentRatios <- merge(FragmentRatios, binCoords, by.x='Bin', by.y='5Mb_bin', all.x=TRUE)
FragmentRatios <- merge(FragmentRatios, recistCategory, by.x='id', by.y='delfiID')

# clean up biopsy data
BloodSampleCharacteristics <- BloodSampleCharacteristics[,.(`Sequence ID`, `DELFI-TF`)]
BloodSampleCharacteristics <- merge(BloodSampleCharacteristics, recistCategory, by.x='Sequence ID', by.y='delfiID')

################################################################################
############# format for plotting ##############################################

# find the central coord between start and end
FragmentRatios[,Center := ((end - start)/2) + start]
FragmentRatios[,`:=`(end = NULL, start = NULL)]

# factor arms to be in a sensible order
chromArmOrder <- mixedsort(unique(FragmentRatios$`Chromosome Arm`))
FragmentRatios[,`Chromosome Arm` := factor(`Chromosome Arm`, levels=chromArmOrder)]

# factor timepoints to be in a sensible order
FragmentRatios[,TimePoint := factor(TimePoint, levels=c("Baseline", "Response", "Stable", "Progression"))]
BloodSampleCharacteristics[,TimePoint := factor(TimePoint, levels=c('Baseline', 'Response', 'Stable', 'Progression'))]

# calculate counts
sampleCounts <- unique(FragmentRatios[,.(id, TimePoint)])[,.N,by=.(TimePoint)]
FragmentRatios <- merge(FragmentRatios, sampleCounts, by='TimePoint')

# manually add some padding to a couple arms for plotting, this is just shifting the data a bit for facet labels
FragmentRatios[Bin == 447, "Center"] <- FragmentRatios[Bin == 447, "Center"] + 10000000
FragmentRatios[Bin == 460, "Center"] <- FragmentRatios[Bin == 460, "Center"] + 25000000
FragmentRatios[Bin == 475, "Center"] <- FragmentRatios[Bin == 475, "Center"] + 10000000

################################################################################
########### Construct Plots ####################################################

# make the actual plot
Fig1C <- ggplot(FragmentRatios, aes(x=Center, y=`Simulated Fragment Ratio`, group=id)) +
  geom_line(alpha=.55, color='#00e6c7') +
  facet_grid(TimePoint + N ~ `Chromosome Arm`, scales="free_x", space="free_x", switch="both") +
  ylab("cfDNA Fragmentation Ratio") +
  ylim(c(-.26, .26)) +
  theme_bw() +
  theme(legend.position = "top", legend.title=element_text(size=28), legend.text=element_text(size=24)) +
  theme(axis.text.x=element_blank(), axis.text.y=element_text(size=24)) +
  theme(axis.ticks.x=element_blank(), axis.title.x=element_blank()) +
  theme(axis.title.y=element_text(size=32)) +
  theme(legend.position = "bottom") +
  theme(legend.key.width = unit(4, 'cm')) +
  theme(panel.grid.minor = element_blank()) +
  theme(strip.text.y.left=element_text(color="white", face="bold", size=24)) +
  theme(strip.text.x.bottom=element_text(color="white", face="bold", size=20)) +
  theme(strip.text.x = element_text(angle=90)) +
  theme(strip.background = element_rect(fill="#1B465A", color="#1B465A")) +
  theme(panel.border=element_blank()) +
  theme(panel.grid=element_blank())

png(filename="Figures/Fig1C.png", height=15, width=32, units="in", res=150)
Fig1C
dev.off()

cairo_pdf(file="Figures/Fig1C.pdf", height=10, width=32)
Fig1C
dev.off()
