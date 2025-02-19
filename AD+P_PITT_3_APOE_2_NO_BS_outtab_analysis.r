##---------------------------------------------------------------------#
##
## Luke Weymouth 10/02/25 - EWAS output analysis
##
##---------------------------------------------------------------------#


### SET UP ###


#Load R
module load R/4.2.0-foss-2021b
R

#Set working directory
setwd("/lustre/projects/Research_Project-MRC164847/Luke_W/APOE_methylation/APOE_meta_analysis/New_QC/EPIC_IDATS/AD+P_PITT/AD+P_PITT_PFC/4_analysis/results/")

#Load R data
load("/lustre/projects/Research_Project-MRC164847/Luke_W/APOE_methylation/APOE_meta_analysis/New_QC/EPIC_IDATS/AD+P_PITT/AD+P_PITT_PFC/4_analysis/results/APOE_2_without_Braak_EWASout.rdat")

#Visualize the results table
head(outtab)

#Cleans up the column names of the results table (removes everything before and including the $ sign)
colnames(outtab) <- sub(".*\\$", "", colnames(outtab))

# Load the manifest file
manifest <- read.csv("/lustre/projects/Research_Project-MRC164847/Luke_W/APOE_methylation/APOE_meta_analysis/New_QC/EPIC_IDATS/EPIC_V1_manifest_file.csv",
                     stringsAsFactors = FALSE)

# Select only the relevant columns from the manifest
manifest_subset <- manifest[, c("IlmnID", "CHR", "MAPINFO", "UCSC_RefGene_Name")]

# Convert outtab to a data frame and keep row names
outtab_df <- data.frame(IlmnID = rownames(outtab), outtab, stringsAsFactors = FALSE)

# Merge outtab with the manifest information
outtab_merged_df <- merge(outtab_df, manifest_subset, by = "IlmnID", all.x = TRUE)

# Print and check first few rows
head(outtab_merged_df)

# Save the updated version
write.csv(outtab_merged_df, "APOE_2_NO_BS_outtab_with_manifest_info.csv", row.names = TRUE)


### RETURN LIST OF NOMINALLY SIGNIFICANT SITES ###


# 1. "Number_of_2_alleles1_P" has a value less than 0.05
APOE2_NO_BS_1_copy_nominally_sig_sites <- outtab_merged_df[outtab_merged_df$Number_of_2_alleles1_P < 0.05, ]

# Print and save
print(APOE2_NO_BS_1_copy_nominally_sig_sites)
write.csv(APOE2_NO_BS_1_copy_nominally_sig_sites, "APOE2_NO_BS_1_copy_nominally_sig_sites.csv", row.names = FALSE)


### RETURN LIST OF BONFERRONI SIGNIFICANT SITES ###


# Calculate Bonferroni-corrected threshold
bonferroni_threshold <- 0.05 / nrow(outtab_merged_df)

# 1. "Number_of_2_alleles1_P" is below the Bonferroni threshold
APOE2_NO_BS_1_copy_bonferroni_sig_sites <- outtab_merged_df[outtab_merged_df$Number_of_2_alleles1_P < bonferroni_threshold, ]

# Print and save
print(APOE2_NO_BS_1_copy_bonferroni_sig_sites)
write.csv(APOE2_NO_BS_1_copy_bonferroni_sig_sites, "APOE2_NO_BS_1_copy_bonferroni_sig_sites.csv", row.names = FALSE)


### PLOTS ###



