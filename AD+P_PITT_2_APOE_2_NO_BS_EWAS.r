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

runEWAS <- function(row, QCmetrics_2_null) {
  APOE_2 <- lm(row ~ QCmetrics_2_null$Number_of_2_alleles + QCmetrics_2_null$Sex + 
                 QCmetrics_2_null$Age + QCmetrics_2_null$PMI)
  
  # Extract coefficients, std.error, and p-value
  results <- c(summary(APOE_2)$coefficients["QCmetrics_2_null$Number_of_2_alleles1", c(1,2,4)])
  
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
QCmetrics_2_null$Number_of_2_alleles <- factor(QCmetrics_2_null$Number_of_2_alleles, levels = c("0", "1"))

# Ensure levels exist
print(levels(QCmetrics_2_null$Number_of_2_alleles))

# Subset beta matrix
normbeta <- normbeta[, QCmetrics_2_null$Basename, drop = FALSE]

#----------------------------------------------------------------------#
# RUN IN PARALLEL ENV
#----------------------------------------------------------------------#

nCores <- detectCores()
cl <- makeCluster(nCores - 1)
registerDoParallel(cl)
clusterExport(cl, list("runEWAS", "QCmetrics_2_null"))

# Run regression in parallel
outtab <- matrix(data = parApply(cl, normbeta, 1, runEWAS, QCmetrics_2_null), ncol = 3, byrow = TRUE)
rownames(outtab) <- rownames(normbeta)

# Stop cluster
stopCluster(cl)

# Colnames
modTerms <- c("QCmetrics_2_null$Number_of_2_alleles1")

generate_colnames <- function(terms) {
  paste0(terms, c("_coeff", "_SE", "_P"))
}

colnames(outtab) <- unlist(lapply(modTerms, generate_colnames))

#----------------------------------------------------------------------#
# SAVE OUTPUT
#----------------------------------------------------------------------#

filePath <- file.path(projDir, "4_analysis/results", paste0("APOE_2_without_Braak", "_EWASout.rdat"))
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

