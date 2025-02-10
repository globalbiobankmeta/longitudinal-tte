#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(argparser)

p <- arg_parser("Create endpoint phenotype files for UKBB")
p <- add_argument(p, "--phecode", help = "")
args <- parse_args(p)

# Read the latest UKBB data
phecode_data <- fread("UKB_Phecode_v1.2b1_ICD_Mapping.txt", header=T, colClasses=rep("character", 6), data.table=F)
diagcode_data <- fread("ukb.hesin_diag.20231204_phecode_age_collapse.txt", header=T, colClasses=c("character", "character", "numeric", "numeric", "character", "character", "numeric"), data.table=F)
death_data <- fread("Death_pheno.20231204_allpops.txt", header=T, data.table=F, select=c("s", "age_death_or_lastvisit"))

# Add death or last visit age to diagnose data
death_unique <- death_data[!duplicated(death_data),]
death_data <- death_unique
colnames(death_data) <- c("eid", "age_death_or_lastvisit")
death_data$eid <- as.character(death_data$eid)
diagcode_data <- diagcode_data %>% left_join(death_data, by = "eid")
diagcode_data <- diagcode_data[which(!is.na(diagcode_data$phecode)), ]
sampleAll <- sort(unique(diagcode_data$eid))

# Extract case samples, i.e. samples with phecode of interest entries
casedata <- diagcode_data[which(diagcode_data$phecode == args$phecode), c("eid", "phecode_age", "age_death_or_lastvisit")] 
colnames(casedata) <- c("eid","age", "age_death_or_lastvisit")
casedata$status <- 1
caselist <- sort(casedata$eid)

# Starting control samples
ctrllist0 <- setdiff(sampleAll, caselist)
ctrldata0 <- diagcode_data[diagcode_data$eid %in% ctrllist0, c("eid", "age_death_or_lastvisit", "phecode")]
# Exclude controls that have similar case phcode
excludecontrol_phecode <- unlist(strsplit(phecode_data$exclude_phecodes[which(phecode_data$phecode == args$phecode)], split=","))
ctrlexcludeeid <- ctrldata0[ctrldata0$phecode %in% excludecontrol_phecode,"eid"] 
# Final control samples
ctrllist <- setdiff(ctrllist0, sort(ctrlexcludeeid))
ctrldata0 <- ctrldata0[,c("eid", "age_death_or_lastvisit")]
ctrldata1 <- ctrldata0[ctrldata0$eid %in% ctrllist, ]
ctrldata <- ctrldata1[!duplicated(ctrldata1),]
colnames(ctrldata) <- c("eid","age")
ctrldata$status <- 0
ctrldata$age_death_or_lastvisit <- ctrldata$age # for final controls, ignore other phecode and ages, keep 1 age_death_or_lastvisit for each sample

data_age <- rbind(casedata, ctrldata)
data_age$event <- data_age$status 
colnames(data_age)[which(colnames(data_age) == "status")] <- paste0("X", args$phecode)
# Got a data frame of all samples for a phecode, case samples have event age and death age; control samples have age_death_or_lastvisit age

# Read covariates, keep gPC 1-20, population, related or no, sex
covdata <- fread("zcat final_samples.txt.bgz", header=T, data.table=F, colClasses=c(rep("character",6), rep("numeric",20), rep("character",2), rep("numeric",5)))
covdata <- covdata[,c("s",paste0("PC",1:20), "pop", "related", "sex")]
phenodata <- covdata %>% left_join(data_age, by = c("s"="eid"))

if(phecode_data$sex[which(phecode_data$phecode == args$phecode)] == "females"){
    phenodata = phenodata[which(phenodata$sex == 0), ]
}
if(phecode_data$sex[which(phecode_data$phecode == args$phecode)] == "males"){
        phenodata = phenodata[which(phenodata$sex == 1), ]
}


# Read sample year of birth
yobdata <- fread("eid_birthdate.txt.onlyyear", header=T, sep="\t", colClasses=c("character", "numeric"), data.table=F)

colnames(yobdata) <- c("eid","birthyear")

phenodata2 <- phenodata %>% inner_join(yobdata, by = c("s" = "eid"))
phenodata2 <- phenodata2[which(!is.na(phenodata2$event)),]
phenodata2 <- phenodata2[which(!is.na(phenodata2$age)),]
phenodata2 <- phenodata2[!duplicated(phenodata2),]

# Write sample count by population and sex to a file
write.table(table(phenodata2$pop, phenodata2$event), paste0("pheno_", args$phecode, "_NbyPop.txt"), col.names = T, row.names = T, quote = F, sep = "\t")
write.table(table(phenodata2$pop, phenodata2$event, phenodata2$sex), paste0("pheno_", args$phecode, "_NbyPopSex.txt"), col.names = T, row.names = T, quote = F, sep = "\t")

# Write final phenotype list of endpoint phecode
for (pop in c("AFR", "AMR", "CSA", "EAS", "EUR", "MID")){
    phenodata2_pop = phenodata2 %>% filter(pop == pop)
    write.table(phenodata2_pop, paste0("pheno_", args$phecode, "_", pop, "_ALL.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
    write.table(phenodata2_pop %>% filter(sex == 0), paste0("pheno_", args$phecode, "_", pop, "_F.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
    write.table(phenodata2_pop %>% filter(sex == 1), paste0("pheno_", args$phecode, "_", pop, "_M.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
}
write.table(phenodata2, paste0("pheno_", args$phecode, ".txt"), quote=F, col.names=T, row.names=F, sep="\t")
