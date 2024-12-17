# longitudinal-tte
This repo provides example scripts to curate Time-to-event phenotypes for GWAS on disease onset and progression based on phecodes. 
- To create endpoint phenotype file of UKBB, `Rscript ukbb_create_surv_phenotype.R --phecode $1`, where `$1` is phecode.
- To create T2E phenotype file of UKBB, `Rscript ukbb_create_surv_secondEvent.R --firstEventFile $1 --firstEvent $2 --secondEventFile $3 --secondEvent $4`, where `$1` and `$3` are paths to output files from the previous script, and `$2` and `$4` are phecodes.
