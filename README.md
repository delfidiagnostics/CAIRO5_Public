# Cancer treatment monitoring using cell-free DNA fragmentomes - Analysis

## Description
This repository holds the analysis scripts for the Nature Communication manuscript entitled "Cancer treatment monitoring using cell-free DNA fragmentomes". The analysis is structured such that the majority of scripts will run off of the supplemental data excel file available at ****** and provided herein for convenience. Occasionally it is necessary to read in data not provided in the supplemental due to the presence of identifiable or other sensitive information, in such cases data has been simulated, sanitized, or otherwise changed to allow for users to still run code, see the next section for more details.

## Auxillary Data Files
- SimulatedFragmentRatios.rds: Contains simulated fragment ratios using a normal distribution mimicing each sample, fragment bins have also been shuffled from the original version.
  - Fig1B and Fig1C
- Fig2b_AuxData.rds: Contains the categorization of samples to their closest RECIST evaluation within a 60 day window.
  - Fig2B
- Fig3a_AuxData.rds: Contains the "confirmed" consecutive RECIST response categories for samples.
  - Fig3A

## Main Figures

* Fig1A - Diagram of sample collection
  - Omitted, constructed manualy

* Fig1B - Heatmap of DELFI-TF model features

* Fig1C - Fragmentation profile across sample RECIST categories

* Fig2A - DELFI-TF model schematic
  - Omitted, constructed manualy

* Fig2B -  Boxplot comparing DELFI-TF amongst Cancer and Controls

* Fig2C - Scatter plot comparing DELFI-TF to MAF

* Fig3A - Boxplot comparing consecutive RECIST evaluations with DELFI-TF and MAF

* Fig3B - Boxplot comparing resection status with DELFI-TF and MAF

* Fig3C - Boxplot comparing metastasis type with DELFI-TF

* Fig3D - Overall survival comparing DELFI-TF and MAF performance synchronicity at baseline

* Fig3E - Multivariate analysis for overall survival

* Fig4A - Overall patient status timeline
  - Specific RECIST date calculations are omitted, overall time on study is shown in lieu

* Fig4B - Progression free survival comparing DELFI-TF performance longitudinally across eligible subjects.

* Fig4C - Overall survival comparing DELFI-TF performance longitudinally across eligible subjects.

## Supplemental Figures

* FigS1A - Study design schematic
  - Omitted, constructed manualy

* FigS2A - Participant count schematic
  - Omitted, constructed manualy

* FigS3A - Fragmentation profiles compared to open and closed A/B compartments
- Omitted, underlying data size limitations, code is provided to replicate analysis but is un-runnable

* FigS4A - Scatter plot comparing DELFI-TF and ichorCNA for samples undetectable by MAF and above presumed LOD thresholds

* FigS4B - Boxplots comparing DELFI-TF v MAF v ichorCNA performance

* FigS5A - Comparison of single patient copy number calls to fragment ratios
  - Omitted

* FigS5B - Comparison of DELFI-TF and MAF to copy number calls for presumed amplified and deleted regions.

* FigS5C - Comparison of DELFI-TF and MAF for a validation cohort of lung cancers.

* FigS6A - Correlation between DELFI-TF/MAF and SLD

* FigS6B - Correlation between DELFI-TF/MAF and CEA

* FigS7A - Boxplot between ever/never progressors for DELFI-TF, MAF, and SLD

* FigS7B - PFS comparing DELFI-TF and MAF performance synchronicity at baseline

* FigS7C - PFS comparing CEA low/high

* FigS8A - Longitudinal comparison of DELFI-TF between ever/never progressors

* FigS8B - Single Subject logitudinal example with RECIST images
  - Omitted, see FigS9A for MAF/DELFI-TF trends

* FigS9A - Longitudinal comparison of DELFI-TF and MAF across subjects (arm 1)
  - Months to surgery is ommited

* FigS9B - Longitudinal comparison of DELFI-TF across subject (arm 3)
  - Months to surgery is ommited

* FigS10A - PFS for durable clinical benefit + DELFI-TF slope trajectory

* FigS10B - PFS for first RECIST assesment

* FigS10C - OS for complete resection + DELFI-TF slope trajectory

## Bulk Re-Analysis
Users wishing to replicate this analysis in bulk can run all scripts using the following command: `ls Fig*.R | awk '{print "Rscript "$1}' | /bin/bash`

The following R libraries will need to be installed through CRAN or Bioconductor, see sessionInfo() for library versions:
```
library(data.table)
library(ggallin)
library(ggplot2)
library(grid)
library(gridExtra)
library(gt)
library(gtable)
library(gtools)
library(readxl)
library(rstatix)
library(scales)
library(survival)
library(survminer)
library(viridis)
```

Further users will need to modify the `setwd()` commands at the top of each script to point to the top level of this repository

## Score Descriptions

### ichorCNA
ichorCNA values were calculated for all 689 samples on study + 153 healthy controls using default parameters, with the exception that the tumor fraction parameter for the EM step was tuned to account for low ctDNA samples `--normal "c(0.5,0.6,0.7,0.8,0.9,0.95,.96,.97,0.98,0.99,0.995,0.999)"`. The panel of normals was derived from 30 healthy samples from [cristiano et al.](https://pubmed.ncbi.nlm.nih.gov/31142840/)

### ddPCR
ddPCR values were obtained from the MT arm of the CAIRO5 study for tissue confirmed RAS/BRAF mutations using BIO-RAD kits #1863506, #12001626, #10049550, #10049047, #12001627, #12001006, #12001037

### DELFI-TF
DELFI-TF is a measure of ctDNA using fragmentomic and aneuploidy based features in a random forest model.

## sessionInfo()
```
R version 4.3.3 (2024-02-29)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Ventura 13.6.6

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] viridis_0.6.5      viridisLite_0.4.2  survminer_0.4.9    ggpubr_0.6.0       survival_3.5-8     scales_1.3.0       rstatix_0.7.2     
 [8] readxl_1.4.3       gtools_3.9.5       gtable_0.3.4       gt_0.10.1          gridExtra_2.3      ggallin_0.1.1.1001 ggplot2_3.4.3     
[15] data.table_1.15.2 

loaded via a namespace (and not attached):
 [1] utf8_1.2.4        generics_0.1.3    tidyr_1.3.1       xml2_1.3.6        lattice_0.22-5    digest_0.6.35     magrittr_2.0.3   
 [8] fastmap_1.1.1     cellranger_1.1.0  Matrix_1.6-5      backports_1.4.1   purrr_1.0.2       fansi_1.0.6       abind_1.4-5      
[15] cli_3.6.2         KMsurv_0.1-5      rlang_1.1.3       munsell_0.5.0     splines_4.3.3     withr_3.0.0       tools_4.3.3      
[22] ggsignif_0.6.4    dplyr_1.1.4       colorspace_2.1-0  km.ci_0.5-6       broom_1.0.5       vctrs_0.6.5       R6_2.5.1         
[29] zoo_1.8-12        lifecycle_1.0.4   car_3.1-2         pkgconfig_2.0.3   pillar_1.9.0      glue_1.7.0        xfun_0.43        
[36] tibble_3.2.1      tidyselect_1.2.1  knitr_1.45        rstudioapi_0.16.0 xtable_1.8-4      survMisc_0.5.6    htmltools_0.5.8  
[43] carData_3.0-5     compiler_4.3.3   
```
