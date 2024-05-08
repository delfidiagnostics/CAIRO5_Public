# ################################################################################
# ############### Setup ##########################################################
# 
# # setwd
# setwd('~/git/CAIRO5_Public')
# 
# # libs
# library(data.table)
# library(scales)
# library(readxl)
# 
# # data
# healthy <- c('CGPLLU272P_NovaSeq_10x', 'CGPLLU273P_NovaSeq_10x', 'CGPLLU276P_NovaSeq_10x', 'CGPLLU277P_NovaSeq_10x',
#              'CGPLLU278P_NovaSeq_10x', 'CGPLLU282P_NovaSeq_10x', 'CGPLLU287P_NovaSeq_10x', 'CGPLLU288P_NovaSeq_10x',
#              'CGPLLU289P_NovaSeq_10x', 'CGPLLU291P_NovaSeq_10x')
# cancerSamps <- c("DL003194CRD0_AMSeq1", "DL003227CRD0_AMSeq1", "DL003189CRD0_AMSeq1",
#                  "DL003200CRD0_AMSeq1", "DL003215CRD0_AMSeq1", "DL003195CRD0_AMSeq1",
#                  "DL003239CRD0_AMSeq1", "DL003181CRD0_AMSeq1", "DL003185CRD0_AMSeq1",
#                  "DL003242CRD0Seq1")
# keep <- c(healthy, cancerSamps)
# 
# referenceBins <- fread('Data/kb100_genomic_coord_ref.txt')
# 
# # chromosome to keep
# keepChr <- 'chr22'
# 
# ################################################################################
# ################# define helper functions ######################################
# 
# # apply gc correction to coverages
# gc.correct <- function(coverage, bias) {
#   # GC content x-axis: from min to max GC-content by .001
#   i <- seq(min(bias, na.rm=TRUE), max(bias, na.rm=TRUE), by = 0.001)
#   coverage.trend <- loess(coverage ~ bias)
#   coverage.model <- loess(predict(coverage.trend, i) ~ i)
#   coverage.pred <- predict(coverage.model, bias)
#   coverage.corrected <- coverage - coverage.pred + median(coverage)
#   return(coverage.corrected)
# }
# 
# # smoothing function
# .meanSmoother <- function(x, k=1, iter=2, na.rm=TRUE){
#   meanSmoother.internal <- function(x, k=1, na.rm=TRUE){
#     n <- length(x)
#     y <- rep(NA,n)
#     
#     window.mean <- function(x, j, k, na.rm=na.rm){
#       if (k>=1){
#         return(mean(x[(j-(k+1)):(j+k)], na.rm=na.rm))
#       } else {
#         return(x[j])
#       }
#     }
#     
#     for (i in (k+1):(n-k)){
#       y[i] <- window.mean(x,i,k, na.rm)
#     }
#     for (i in 1:k){
#       y[i] <- window.mean(x,i,i-1, na.rm)
#     }
#     for (i in (n-k+1):n){
#       y[i] <- window.mean(x,i,n-i,na.rm)
#     }
#     y
#   }
#   
#   for (i in 1:iter){
#     x <- meanSmoother.internal(x, k=k, na.rm=na.rm)
#   }
#   x
# }
# 
# ################################################################################
# ################# calculate fragment profiles ##################################
# 
# # calculate short profiles
# buildProfile <- function(sampleID, method=c("short", "di", "mon")){
#   
#   # read in data
#   countData <- fread(paste0("Data/counts/", sampleID, ".counts.tsv.gz"))
#   gcData <- fread(paste0("Data/gc/", sampleID, ".gc.tsv.gz"))
#   
#   # select method, short, mono, di
#   if(method == "short") {
#     countData <- countData[width >= 100 & width <= 150]
#   } else if(method == 'di') {
#     countData <- countData[width >= 260 & width <= 420]
#   } else if(method == 'mon') {
#     countData <- countData[width >= 100 & width <= 220]
#   } else {
#     stop("ERROR")
#   }
#   
#   # aggregate count, combine with gc, and apply gc correction
#   countData <- countData[,.(raw_count = sum(count, na.rm=T)), by=.(bin_num)]
#   gcData <- merge(gcData, countData, by="bin_num", all.x=TRUE)
#   gcData[is.na(raw_count), raw_count := 0]
#   gcData <- gcData[,.(sample_id, bin_num, armlevel, gc, frag.gc, raw_count)]
#   gcData[,corrected := gc.correct(raw_count, frag.gc)]
#   
#   return(gcData)
# }
# short <- lapply(keep, buildProfile, method="short")
# short <- rbindlist(short)
# di <- lapply(keep, buildProfile, method="di")
# di <- rbindlist(di)
# mon <- lapply(keep, buildProfile, method="mon")
# mon <- rbindlist(mon)
# 
# ################################################################################
# ################ format the coad and hic compartments ##########################
# 
# data <- fread("https://raw.githubusercontent.com/Jfortin1/TCGA_AB_Compartments/master/data/coad_tumor_compartments_100kb.txt")
# data <- data[,.(chr, start, end, eigen)]
# setnames(data, 'eigen', 'coad')
# data <- merge(data, referenceBins, by=c("chr", "start", "end"))
# data <- data[,.(bin_num, arm, coad)]
# data <- data[order(bin_num)]
# coad.dat <- data[,.(coad=.meanSmoother(coad), bin_num = bin_num), by=.(arm)]
# 
# data <- fread("https://raw.githubusercontent.com/Jfortin1/HiC_AB_Compartments/master/data/hic_compartments_100kb_ebv_2014.txt")
# data <- data[,.(chr, start, end, eigen)]
# setnames(data, 'eigen', 'lymph')
# data <- merge(data, referenceBins, by=c("chr", "start", "end"))
# data <- data[,.(bin_num, arm, lymph)]
# data <- data[order(bin_num)]
# lymph.dat <- data[,.(lymph=.meanSmoother(lymph), bin_num = bin_num), by=.(arm)]
# 
# ############################
# # Construct fragmentation profiles for tumor and healthy
# ##############################
# 
# #calculate fragmentation profile values
# fps <- short %>% select(sample_id, bin_num, armlevel, short = corrected) %>% 
#   inner_join(mon %>% select(sample_id, bin_num, armlevel, mononucleosome = corrected)) %>%
#   inner_join(di %>% select(sample_id, bin_num, armlevel, dinucleosome = corrected)) %>% 
#   mutate(sl.ratio = short/(mononucleosome-short), 
#          mono.cov = mononucleosome, 
#          mono.di = mononucleosome/dinucleosome) %>% 
#   group_by(sample_id) %>% 
#   mutate(centered.sl.ratio = scale(sl.ratio, center = T, scale = F)[,1], 
#          scaled.mono.cov = scale(mono.cov)[,1], 
#          scaled.mono.di = scale(mono.di)[,1]) %>% 
#   ungroup() %>% 
#   select(sample_id, bin_num, armlevel, centered.sl.ratio, scaled.mono.cov, scaled.mono.di)
# 
# 
# #median healthy
# ref.healthy <- fps %>% filter(grepl("CG", sample_id)) %>% 
#   group_by(armlevel, bin_num) %>% 
#   summarize(ref.sl.ratio = median(centered.sl.ratio), 
#             ref.mono.cov = median(scaled.mono.cov), 
#             ref.mono.di = median(scaled.mono.di)) %>% ungroup() %>% 
#   group_by(armlevel) %>% 
#   mutate(smoothed.sl.ratio = .meanSmoother(ref.sl.ratio), 
#          smoothed.mono.cov = .meanSmoother(ref.mono.cov), 
#          smoothed.mono.di = .meanSmoother(ref.mono.di))
# 
# obs.tumor <- fps %>% filter(grepl("DL", sample_id)) %>% 
#   group_by(armlevel, bin_num) %>% 
#   summarize(med.sl.ratio = median(centered.sl.ratio), 
#             med.mono.cov = median(scaled.mono.cov), 
#             med.mono.di = median(scaled.mono.di)) %>% ungroup() %>% 
#   group_by(armlevel) %>% 
#   mutate(smoothed.sl.ratio = .meanSmoother(med.sl.ratio), 
#          smoothed.mono.cov = .meanSmoother(med.mono.cov), 
#          smoothed.mono.di = .meanSmoother(med.mono.di))
# 
# #calculate the tumor component per tumor sample, 
# #then take the median
# tumor.components <- fps%>% filter(grepl("DL", sample_id)) %>% 
#   inner_join(BloodSampleCharacteristics %>% select(sample_id = delfiID.y, dtf = `DELFI-TF`)) %>% 
#   inner_join(ref.healthy) %>% 
#   mutate(sl.ratio.tumor.component = (centered.sl.ratio - ((1-dtf) * ref.sl.ratio))/ dtf, 
#          mono.cov.tumor.component = (scaled.mono.cov - ((1-dtf) * ref.mono.cov))/ dtf, 
#          mono.di.tumor.component = (scaled.mono.di - ((1-dtf) * ref.mono.di))/ dtf) %>% 
#   group_by(bin_num, armlevel) %>%
#   summarize(med.sl.ratio = median(sl.ratio.tumor.component), 
#             med.mono.cov = median(mono.cov.tumor.component), 
#             med.mono.di = median(mono.di.tumor.component)) %>% ungroup() %>% 
#   group_by(armlevel) %>% 
#   mutate(smoothed.sl.ratio = .meanSmoother(med.sl.ratio), 
#          smoothed.mono.cov = .meanSmoother(med.mono.cov), 
#          smoothed.mono.di = .meanSmoother(med.mono.di))
# 
# 
# final.plot.data <- rbind(lymph.dat %>% mutate(type = "lymphoblastoid") %>% select(arm, bin_num, eigen = lymph, type), 
#                          coad.dat %>% mutate(type = "coad") %>% select(arm, bin_num, eigen = coad, type), 
#                          ref.healthy %>% mutate(type = "reference.noncancer") %>% select(arm = armlevel, bin_num, eigen = smoothed.mono.di, type), 
#                          obs.tumor %>% mutate(type = "observed.cancer") %>% select(arm = armlevel, bin_num, eigen = smoothed.mono.di, type), 
#                          tumor.components %>% mutate(type = "tumor.component") %>% select(arm = armlevel, bin_num, eigen = smoothed.mono.di, type)) %>% 
#   spread(key = type, value = eigen) %>% 
#   filter(!is.na(coad), !is.na(reference.noncancer)) %>% 
#   gather(key = type, value = eigen, -arm, -bin_num) %>% 
#   mutate(type = factor(type, levels = c("coad", "lymphoblastoid", "tumor.component", "observed.cancer", "reference.noncancer")))
# 
# highlight <- final.plot.data %>% 
#   mutate(class = ifelse(eigen > 0, "Closed", "Open")) %>% 
#   select(-eigen) %>%
#   spread(key = type, value = class) %>% 
#   filter(!is.na(tumor.component), !is.na(coad)) %>% 
#   mutate(highlight = ifelse(tumor.component == coad & reference.noncancer == lymphoblastoid & coad != lymphoblastoid, "Yes", "No"))
# 
# 
# final.plot.data <- final.plot.data %>% 
#   inner_join(highlight %>% select(arm, bin_num, highlight)) %>%
#   filter(arm == '20p') %>%
#   mutate(state = ifelse(eigen > 0, "Closed", "Open"))
# 
# final.plot.data <- as.data.table(final.plot.data)
# 
# # rename type
# final.plot.data[type == 'lymphoblastoid', type := 'Lymphoblastoid cell Hi-C reference A/B compartments']
# final.plot.data[type == 'observed.cancer', type := 'CAIRO5 CRC patient cfDNA']
# final.plot.data[type == 'reference.noncancer', type := 'Individuals without cancer cfDNA']
# final.plot.data[type == 'coad', type := 'TCGA CRC Tissue Reference A/B Compartments']
# final.plot.data[type == 'tumor.component', type := 'CAIRO5 Extracted CRC patient cfDNA component']
# final.plot.data[,type := factor(type, levels=c('TCGA CRC Tissue Reference A/B Compartments', 'CAIRO5 Extracted CRC patient cfDNA component',
#                                                'CAIRO5 CRC patient cfDNA', 'Individuals without cancer cfDNA', 'Lymphoblastoid cell Hi-C reference A/B compartments'))]
# 
# ################################################################################
# ################### make the final plot ########################################
# 
# # make plots
# FigS3A <- ggplot() +
#   geom_bar(data=final.plot.data, aes(x=as.character(bin_num), y=eigen, fill=state, alpha=highlight), stat="identity") +
#   facet_wrap(~type, ncol = 1) + 
#   scale_alpha_manual(values=c('Yes'=.95, 'No'=.1), guide="none") + 
#   scale_fill_manual('A/B Compartments',values=c('Open'='gray50', 'Closed'='red4')) + 
#   ylab("Fragmentation Profile") +
#   xlab("Chromosome 20p") +
#   theme_bw() +
#   theme(panel.grid=element_blank()) +
#   theme(strip.text=element_text(color='white', face='bold', size=20), strip.background = element_rect(fill="#1B465A")) +
#   theme(legend.position = 'top') +
#   theme(legend.key.size = unit(1, 'cm')) +
#   theme(legend.title = element_text(size=20), legend.text=element_text(size=18)) +
#   theme(axis.text.y=element_text(size=18), axis.title.y=element_text(size=20)) +
#   theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_text(size=20))
# 
# # save plots
# png(filename="Figures/FigS3A.png", height=14, width=20, units="in", res=150)
# FigS3A
# dev.off()
# 
# pdf(file='Figures/FigS3A.pdf', height=14, width=20)
# FigS3A
# dev.off()
