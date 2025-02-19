##---------------------------------------------------------------------#
##
## Title: Run EWAS
##
## Purpose of script: Run linear regression on normalised data
##                    
##                    This scripted is adapted from one co-authored
##                    by EJH and EMW for the MRC schizophrenia project
##
##                   
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

# to run: sbatch..... 

# within cell Line models

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(lme4)
library(lmerTest)
library(dplyr)
library(doParallel)

#----------------------------------------------------------------------#
# DEFINE ANALYSIS FUNCTION
#----------------------------------------------------------------------#

runEWAS <- function(row, QCmetrics_4_BS) {
  APOE_4 <- lm(row ~ QCmetrics_4_BS$Number_of_4_alleles + QCmetrics_4_BS$Sex + 
                 QCmetrics_4_BS$Age + QCmetrics_4_BS$PMI + QCmetrics_4_BS$Braak_Stage_AD)
  
  # Extract coefficients, std.error, and p-value
  results <- c(summary(APOE_4)$coefficients["QCmetrics_4_BS$Number_of_4_alleles1", c(1,2,4)],
               summary(APOE_4)$coefficients["QCmetrics_4_BS$Number_of_4_alleles2", c(1,2,4)])
  
  return(results)
}

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

projDir <- "/lustre/projects/Research_Project-MRC164847/Luke_W/APOE_methylation/APOE_meta_analysis/New_QC/EPIC_IDATS/AD+P_PITT/AD+P_PITT_PFC/"
normData <- file.path(projDir, "3_normalised/normalised.rdata")

#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

setwd(projDir)
load(normData)

print(paste0("running EWAS"))

# Change to factor and order levels
QCmetrics_4_BS$Number_of_4_alleles <- factor(QCmetrics_4_BS$Number_of_4_alleles, levels = c("0", "1", "2"))

# Ensure levels exist
print(levels(QCmetrics_4_BS$Number_of_4_alleles))

# Subset beta matrix
normbeta <- normbeta[, QCmetrics_4_BS$Basename, drop = FALSE]

#----------------------------------------------------------------------#
# RUN IN PARALLEL ENV
#----------------------------------------------------------------------#

nCores <- detectCores()
cl <- makeCluster(nCores - 1)
registerDoParallel(cl)
clusterExport(cl, list("runEWAS", "QCmetrics_4_BS"))

# Run regression in parallel
outtab <- matrix(data = parApply(cl, normbeta, 1, runEWAS, QCmetrics_4_BS), ncol = 6, byrow = TRUE)
rownames(outtab) <- rownames(normbeta)

# Stop cluster
stopCluster(cl)

# Colnames
modTerms <- c("QCmetrics_4_BS$Number_of_4_alleles1", "QCmetrics_4_BS$Number_of_4_alleles2")

generate_colnames <- function(terms) {
  paste0(terms, c("_coeff", "_SE", "_P"))
}

colnames(outtab) <- unlist(lapply(modTerms, generate_colnames))

#----------------------------------------------------------------------#
# SAVE OUTPUT
#----------------------------------------------------------------------#

filePath <- file.path(projDir, "4_analysis/results", paste0("APOE_4_with_Braak", "_EWASout.rdat"))
save(outtab, file = filePath)

#----------------------------------------------------------------------#
# Explore results
#----------------------------------------------------------------------#

bonfP <- 0.05 / nrow(outtab)

# Function to count significant results
countSig <- function(x) {
  sig <- sum(x < bonfP)
  return(sig)
}

print("Script completed successfully!")

