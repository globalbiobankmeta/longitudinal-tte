#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(argparser)

p <- arg_parser("Create T2E phenotype files for UKBB")
p <- add_argument(p, "--firstEventFile", help = "")
p <- add_argument(p, "--firstEvent", help = "")
p <- add_argument(p, "--secondEventFile", help = "")
p <- add_argument(p, "--secondEvent", help = "")
#p <- add_argument(p, "--excludeEventFile", help = "")


args <- parse_args(p)

firstdata <- fread(args$firstEventFile,  header=T, data.table=F)
seconddata <- fread(args$secondEventFile, header=T,  data.table=F)
# Extract samples with case status of t0, i.e. "censors"
phenodatanew = firstdata %>% filter(event == 1)

if(("event" %in% colnames(seconddata))){
  seconddata = seconddata[which(seconddata$event == 1), ]
  seconddata = seconddata[,c("s","age")]
}
phenodatanew2 <- phenodatanew %>% left_join(seconddata, by = "s")
phenodatanew2$secondEvent <- rep(0,nrow(phenodatanew2))
# Samples with case status of t1 event, i.e. "cases"
phenodatanew2$secondEvent[which(!is.na(phenodatanew2$age.y))] <- 1
phenodatanew2$diagAge <- phenodatanew2$age.x
#age.x comes from phenodatanew, which are "censors", t0 cases; age.y comes from seconddata, NA if controls, age for cases
phenodatanew2$secondTime <- phenodatanew2$age.y - phenodatanew2$diagAge # number of years to develop second event
# for those didn't develop second event, secondTime is number of years to death from censor event
phenodatanew2$secondTime[which(is.na(phenodatanew2$secondTime))] <- phenodatanew2$age_death_or_lastvisit[which(is.na(phenodatanew2$secondTime))] - phenodatanew2$diagAge[which(is.na(phenodatanew2$secondTime))]


phenodatanew2 <- phenodatanew2[which(!is.na(phenodatanew2$secondTime)), ]
phenodatanew2 <- phenodatanew2[which(phenodatanew2$secondTime >0), ]

# if(args$excludeEventFile != ""){
# 	phenodataexclude = fread(args$excludeEventFile,  header=T, data.table=F)
# 	excludesamplelist = phenodataexclude$s[which(phenodataexclude$event == 1)]
# 	phenodatanew2 = phenodatanew2[which(!(phenodatanew2$s %in% excludesamplelist)),]
# }

# Write sample count by population and T2E status to a file
write.table(table(phenodatanew2$pop, phenodatanew2$secondEvent), paste0("pheno_", args$firstEvent, "_to_", args$secondEvent, "_NbyPop.txt"), col.names = T, row.names = T, quote = F, sep = "\t")
write.table(table(phenodatanew2$pop, phenodatanew2$secondEvent, phenodatanew2$sex), paste0("pheno_", args$firstEvent, "_to_", args$secondEvent, "_NbyPopSex.txt"), col.names = T, row.names = T, quote = F, sep = "\t")

# Write final phenotype list of T2E
write.table(phenodatanew2, paste0("pheno_", args$firstEvent, "_to_", args$secondEvent, ".txt"), quote=F, col.names=T, row.names=F, sep = "\t")

for (POP in c("AFR", "AMR", "CSA", "EAS", "EUR", "MID")){
    phenodatanew2_pop = phenodatanew2 %>% filter(pop == POP)
    write.table(phenodatanew2_pop, paste0("pheno_", args$firstEvent, "_to_", args$secondEvent, "_", POP, "_ALL.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
    write.table(phenodatanew2_pop %>% filter(sex == 0), paste0("pheno_", args$firstEvent, "_to_", args$secondEvent, "_", POP, "_F.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
    write.table(phenodatanew2_pop %>% filter(sex == 1), paste0("pheno_", args$firstEvent, "_to_", args$secondEvent, "_", POP, "_M.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
}