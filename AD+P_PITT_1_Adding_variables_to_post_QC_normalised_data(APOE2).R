## LUKE WEYMOUTH ## 03/02/2025 ## ADDING VARIABLES TO POST QC NORMALISED DATA (QCmetrics) 

## Load R
module load R/4.2.0-foss-2021b
R

## Change working directory
setwd("/lustre/projects/Research_Project-MRC164847/Luke_W/APOE_methylation/APOE_meta_analysis/New_QC/EPIC_IDATS/AD+P_PITT/AD+P_PITT_PFC/3_normalised/")

## Load existing objects from normalised.rdata
load("normalised.rdata")
ls()
#"normbeta"  "QCmetrics"
head(QCmetrics)

## Install and load required package
install.packages("readxl")
library(readxl)

## Load additional variables
ADD_VARS <- read_excel("AD+P_PITT_ADD_COLUMNS.xlsx")
head(ADD_VARS)

## Ensure column names match for merging
colnames(ADD_VARS)[1] <- colnames(QCmetrics)[3]

## Merge additional variables into QCmetrics
QCmetrics_2_BS <- merge(QCmetrics, ADD_VARS, 
                        by.x = colnames(QCmetrics)[3],  
                        by.y = colnames(ADD_VARS)[1],  
                        all.x = TRUE)

## Remove rows with missing values in 'Number_of_2_alleles'
QCmetrics_2_BS <- QCmetrics_2_BS[!is.na(QCmetrics_2_BS$Number_of_2_alleles), ]

## Create an additional copy for EWAS without Braak stage as a cofactor 
QCmetrics_2_null <- QCmetrics_2_BS

## Save QCmetrics_2_BS into normalised.rdata alongside existing objects
save(list = ls(), file = "normalised.rdata")