# Read list of phecodes and find all child phecodes which map to it
getChildPhecodes <- function(phenotypes) {
  allPhecodes <- data.table("parentPhecode" = phenotypes,
                            "withChildren" = vector(mode = "list"))
  for (p in seq_along(phenotypes)) {
    unroll <- subset(PheWAS::phecode_rollup_map, 
                     phecode_unrolled %in% phenotypes[p])
    allPhecodes$withChildren[[p]] <- unroll$code
  }
  allPhecodes
}

# Read in inclusion phecodes and use the PheWAS package's built-in 
# exclusion mapping to create columns holding exclusion phecodes
getExclusionPhecodes <- function(phecodeInfo) {
  phecodeInfo$exclusionPhecodes <- vector(mode = "list", 
                                          length = nrow(phecodeInfo))
  phecodeInfo$numExclusionPhecodes <- numeric(nrow(phecodeInfo))
  for (p in 1:nrow(phecodeInfo)) {
    exclusions <- mapPhecodesToExclusions(phecodes = phecodeInfo$withChildren[[p]])
    phecodeInfo$exclusionPhecodes[[p]] <- unique(exclusions$exclusion)
    phecodeInfo$numExclusionPhecodes[p] <- uniqueN(exclusions$exclusion)
  }
  phecodeInfo
}

# Read in vectors of *case inclusion* phecodes and a mapping file, and create 
# columns holding the ICDs they map to, and a column holding the number of 
# unique ICDs they map to
getCaseInclusionIcds <- function(phecodeInfo, phecodeIcdMapping) {
  phecodeInfo$caseInclusionIcds <- vector(mode = "list", 
                                          length = nrow(phecodeInfo))
  phecodeInfo$numInclusionIcds <- numeric(nrow(phecodeInfo))
  for (p in 1:nrow(phecodeInfo)) {
    mappings <- subset(phecodeIcdMapping,
                       Phecode %in% phecodeInfo$withChildren[[p]])
    phecodeInfo$caseInclusionIcds[[p]] <- unique(mappings$ICD)
    phecodeInfo$numInclusionIcds[p] <- uniqueN(mappings$ICD)

    # Also record the number of inclusion ICDs mapped by the Google sheet
    googleInclusionIcds <- unlist(strsplit(googleMappings$ICD_case_include[p], split = ","))
    phecodeInfo$googleNumInclusionIcds[p] <- length(googleInclusionIcds)
  }
  phecodeInfo
}

# Read in vectors of *control exclusion* phecodes and a mapping file, 
# union all of the ICDs they map to, and create a column holding the result 
# and a column holding the number of unique ICDs they map to
getControlExclusionIcds <- function(phecodeInfo, phecodeIcdMapping) {
  phecodeInfo$controlExclusionIcds <- vector(mode = "list", 
                                             length = nrow(phecodeInfo))
  phecodeInfo$numExclusionIcds <- numeric(nrow(phecodeInfo))
  for (p in 1:nrow(phecodeInfo)) {
    mappings <- subset(phecodeIcdMapping,
                       Phecode %in% phecodeInfo$exclusionPhecodes[[p]])
    phecodeInfo$controlExclusionIcds[[p]] <- unique(mappings$ICD)
    phecodeInfo$numExclusionIcds[p] <- uniqueN(mappings$ICD)

    # Also record the number of exclusion ICDs mapped by the Google sheet
    googleExclusionIcds <- unlist(strsplit(googleMappings$ICD_control_exclude[p], split = ","))
    phecodeInfo$googleNumExclusionIcds[p] <- length(googleExclusionIcds)
  }
  phecodeInfo
}
