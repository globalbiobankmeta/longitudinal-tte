library(data.table)
library(PheWAS)
source("/nfs/turbo/precision-health/DataDirect/HUM00126227_Genetics_associated_with_CVD/pheWASResources/helpers.R")
source("/nfs/turbo/precision-health/DataDirect/HUM00126227_Genetics_associated_with_CVD/pheWASResources/customMapPhecodesToExclusions.R")
# Use list provided on Google Sheets just to get the list of study phecodes
googleMappings <- fread("/nfs/turbo/precision-health/DataDirect/HUM00126227_Genetics_associated_with_CVD/Time-to-event_Phenotype_sample_sizes - _Corrected_ Phecode_ICD_map.csv",
                        colClasses = rep("character",5))

# Find all child phecodes that map to the study phecode (rollup = 1)
# Uses PheWAS::phecode_rollup_map
phewasMappings <- getChildPhecodes(phenotypes = googleMappings$Phecode)

# Find all control exclusion phecodes
# Uses a custom version of PheWAS::mapPhecodesToExclusions()
phewasMappings <- getExclusionPhecodes(phewasMappings)

# Get ICDs for all case inclusions
phewas12Map <- fread("/nfs/turbo/precision-health/DataDirect/HUM00126227_Genetics_associated_with_CVD/pheWASResources/Phecode_map_v1_2_icd9_icd10cm.csv",
                     colClasses = c("character", "integer", rep("character",4)))
phewasMappings <- getCaseInclusionIcds(phewasMappings, phewas12Map)

# Get ICDs for all control exclusions
phewasMappings <- getControlExclusionIcds(phewasMappings, phewas12Map)

# For the phenotypes with a lot (more than 1000) fewer cases than Google sheet 
# version, add on the 1-2 missing case inclusion ICDs
# Heart failure NOS
hfIdx <- which(googleMappings$Phecode == "428.2")
googleHfExclusiveIcds <- setdiff(unlist(strsplit(googleMappings$ICD_case_include[hfIdx], ",")), phewasMappings$caseInclusionIcds[[hfIdx]])
# # Which heart failure ICDs are exclusive to phewas 1.2? Answer: none
# setdiff(phewasMappings$caseInclusionIcds[[hfIdx]], unlist(strsplit(googleMappings$ICD_case_include[hfIdx], ",")))
# Concatenate the phewas mappings for HF with the ICD codes only Google sheet has
phewasMappings$caseInclusionIcds[[hfIdx]] <- c(phewasMappings$caseInclusionIcds[[hfIdx]], googleHfExclusiveIcds)
# Update count to avoid confusion 
phewasMappings$numInclusionIcds[[hfIdx]] <- length(phewasMappings$caseInclusionIcds[[hfIdx]])
# COPD
copdIdx <- which(googleMappings$Phecode == "496.21")
googleCopdExclusiveIcds <- setdiff(unlist(strsplit(googleMappings$ICD_case_include[copdIdx], ",")), phewasMappings$caseInclusionIcds[[copdIdx]])
# # Which COPD ICDs are exclusive to phewas 1.2? Answer: "491.20"
# setdiff(phewasMappings$caseInclusionIcds[[copdIdx]], unlist(strsplit(googleMappings$ICD_case_include[copdIdx], ",")))
# Concat the COPD Google mappings onto phewas mappings
phewasMappings$caseInclusionIcds[[copdIdx]] <- c(phewasMappings$caseInclusionIcds[[copdIdx]], googleCopdExclusiveIcds)
phewasMappings$numInclusionIcds[[copdIdx]] <- length(phewasMappings$caseInclusionIcds[[copdIdx]])
# Concatenate the case inclusion column contents per row, collapsing using a comma
phewasMappings$caseInclusionIcdsCollapsed <- sapply(phewasMappings$caseInclusionIcds, function(x) paste(x, collapse = ","))
# Concatenate the control exclusion column contents per row, collapsing using a comma
phewasMappings$controlExclusionIcdsCollapsed <- sapply(phewasMappings$controlExclusionIcds, function(x) paste(x, collapse = ","))

# Create a new table with all the new mappings concatenated into a string, with
# one column for case inclusions and one for control exclusions.
# Also include the count of inclusion and exclusion ICDs for the new mappings,
# and place these columns next to the counts for the Google Sheets mappings.
mappingSummary <- data.table("Phecode" = phewasMappings$parentPhecode,
                             "Phecode_Description" = googleMappings$Phecode_Description,
                             "sex" = googleMappings$Sex,
                             "ICD_case_include" = phewasMappings$caseInclusionIcdsCollapsed,
                             "ICD_control_exclude" = phewasMappings$controlExclusionIcdsCollapsed,
                             "google_case_include" = googleMappings$ICD_case_include,
                             "google_control_exclude" = googleMappings$ICD_control_exclude,
                             "numInclusionIcds" = phewasMappings$numInclusionIcds,
                             "numExclusionIcds" = phewasMappings$numExclusionIcds,
                             "googleNumInclusionIcds" = phewasMappings$googleNumInclusionIcds,
                             "googleNumExclusionIcds" = phewasMappings$googleNumExclusionIcds)


# Create a new mapping table with the results, only including the equivalent 
# columns to the Google Sheets mapping file and using the same column names
newMapping <- data.table("Phecode" = phewasMappings$parentPhecode,
                         "Phecode_Description" = googleMappings$Phecode_Description,
                         "sex" = googleMappings$Sex,
                         "ICD_case_include" = phewasMappings$caseInclusionIcdsCollapsed,
                         "ICD_control_exclude" = phewasMappings$controlExclusionIcdsCollapsed)

# Write the new mapping table to a csv file
# This is the version without the ICDs manually added to heart failure and COPD 
# case inclusion lists
# fwrite(newMapping, "/nfs/turbo/precision-health/DataDirect/HUM00126227_Genetics_associated_with_CVD/pheWASResources/phewas1_2_mapping.csv",)

# This version has more ICDs in heart failure and COPD case inclusion
fwrite(newMapping, "/nfs/turbo/precision-health/DataDirect/HUM00126227_Genetics_associated_with_CVD/pheWASResources/phewas1_2_mappingPlus.csv",)