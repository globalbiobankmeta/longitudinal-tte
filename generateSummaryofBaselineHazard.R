#!/usr/bin/env Rscript

options(stringsAsFactors=F)

## load R libraries
library(survival)
require(optparse) #install.packages("optparse")

print(sessionInfo())

## set list of cmd line arguments
option_list <- list(
  make_option("--step1outputFile", type="character",default="",
    help="Path to step1 output file. .rda"),
  make_option("--outputFile_for_baselineHazard", type="character",default="",
    help="Path to output the baseline Hazard to an R .rda file"))

## list of options
parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)



#load step1 output to R
#modglmm <- readRDS(opt$step1outputFile)
load(opt$step1outputFile)
# Sample survival data
# Survival time
time <- modglmm$eventTime
# Event status (1 = event, 0 = censored)
status <- modglmm$y

# Create a survival object
surv_obj <- Surv(time, status)

# Fit Kaplan-Meier model
km_fit <- survfit(surv_obj ~ 1)

# Extract summary statistics
km_summary <- summary(km_fit)

# Create a summary data frame
 km_data <- data.frame(
  time = km_summary$time,
  at_risk = km_summary$n.risk,
  events = km_summary$n.event,
  survival_prob = km_summary$surv  
)

# Print summary statistics
print(head(km_data))

#  time at_risk events survival_prob
#1    2       7      1      0.85714
#2    8       6      2      0.57143
#3   12       3      1      0.38095
#4   15       2      1      0.19048

modelsummary = list(km_data = km_data, 
minEventTime = modglmm$minEventTime, 
eventTimeBinSize = modglmm$eventTimeBinSize)

#Generate summary data for plotting a Kaplan-Meier curve
save(modelsummary, file=opt$outputFile_for_baselineHazard)
