library(data.table)
library(dplyr)

# phecode definitions and ICD mappings
# downloaded from https://phewascatalog.org/phewas/#phe12 on March 3, 2026
phecodeDefinitions <- fread("/nfs/turbo/precision-health/DataDirect/HUM00126227_Genetics_associated_with_CVD/pheWASResources/phecode_definitions1.2.csv",
                            colClasses = c(rep("character",4), rep("integer",3), "character"))
pheIcdMap <- fread("/nfs/turbo/precision-health/DataDirect/HUM00126227_Genetics_associated_with_CVD/pheWASResources/Phecode_map_v1_2_icd9_icd10cm.csv",
                   colClasses = c("character", "integer", rep("character",4)))

# Note: some ICDs map to more than one Phecode
uniqueN(pheIcdMap$ICD) - nrow(pheIcdMap)

# List of study phenotypes
correctedMappings <- fread("/nfs/turbo/precision-health/DataDirect/HUM00126227_Genetics_associated_with_CVD/Time-to-event_Phenotype_sample_sizes - _Corrected_ Phecode_ICD_map.csv",
                           colClasses = rep("character", 5))
studyPhecodes <- correctedMappings$Phecode

# Subset phecode definitions file to our study's phenotypes
phecodeDefinitions <- phecodeDefinitions %>%
  filter(phecode %in% studyPhecodes) %>%
  select(-rollup, -leaf, -category_number)

# Subset phecode ICD map file to phecodes in our study
pheIcdMap <- pheIcdMap %>%
  filter(Phecode %in% studyPhecodes) %>%
  select(-PhecodeString, -PhecodeCategory)

table(pheIcdMap$Flag)

# Why is it so much smaller than the descriptions tab in the Google sheet?
# is it because we've only listed case inclusion criteria? Answer: no
### phecode-ICD description list (1 ICD per row)
comprehensivePhecodeMappings <- fread("/nfs/turbo/precision-health/DataDirect/HUM00126227_Genetics_associated_with_CVD/Time-to-event_Phenotype_sample_sizes - Phecode_ICD_description.csv",
                                      colClasses = rep("character",7))
sum(comprehensivePhecodeMappings$ICD_category == "ICD_case_include") # 3011

# Compare one phecode
sum(comprehensivePhecodeMappings$Phecode == "165.1" & comprehensivePhecodeMappings$ICD_category == "ICD_case_include") # 123
sum(pheIcdMap$Phecode == "165.1") # 44
# The PheWAS 1.2 mapping seems to be much more compact?
# Let's compare all phecodes
pheIcdComparison <- data.table("phecode" = studyPhecodes,
                               "oldDescriptions" = numeric(length(studyPhecodes)),
                               "phewas1.2" = numeric(length(studyPhecodes)))
for (p in seq_along(studyPhecodes)) {
  pheIcdComparison$oldDescriptions[p] = sum(comprehensivePhecodeMappings$Phecode == studyPhecodes[p] &
                                              comprehensivePhecodeMappings$ICD_category == "ICD_case_include")
  pheIcdComparison$phewas1.2[p] = sum(pheIcdMap$Phecode == studyPhecodes[p])
}

# I manually verified that number of ICD mappings is correct for Phecodes 585.3 and 440
# Looks like there are duplicate ICD codes in the descriptions sheet, e.g. 440.2
# is listed twice with different descriptions
icdDescriptionCounts <- table(comprehensivePhecodeMappings$ICD)
sum(icdDescriptionCounts > 1)
icdDescriptionCounts[icdDescriptionCounts > 1]
# How many duplicates in total? Subtract one for the first instance of each
sum(icdDescriptionCounts[icdDescriptionCounts > 1] - 1)
# Duplicates still don't explain much of the difference in number of mappings
nrow(comprehensivePhecodeMappings) - sum(icdDescriptionCounts[icdDescriptionCounts > 1] - 1) # 6317

