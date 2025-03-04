# GBMI Time-to-Event GWAS
This repository provides example scripts for curating **time-to-event phenotypes** and running **GWAS on disease onset and progression** using biobank data. The examples use UK Biobank (UKBB) data, but the approach can be adapted to other biobanks.

## 1. Disease Onset Phenotype Curation 
To generate **onset phenotype files**, use the following command:
```
Rscript ukbb_create_surv_phenotype.R --phecode ${phecode_of_interest}
```
### Example Output
See `example_ukbb_pheno_files/pheno_411.2_EUR_ALL.txt.gz` for an example onset phenotype file. The Rscript will write population- and sex-stratified phenotype files. A file listing case and control sample sizes for each population will also be generated.

## 2. Disease Progression Phenotype Curation
To generate **progression phenotype files**, use the following command:
```
Rscript ukbb_create_surv_secondEvent.R \
    --firstEventFile ${T0_phenotype_file} \
    --firstEvent ${T0_phecode} \
    --secondEventFile ${T1_phenotype_file} \
    --secondEvent ${T1_phecode}
```
### Example Output
See `example_ukbb_pheno_files/pheno_411.2_to_428.2_EUR_ALL.txt.gz` for an example progression phenotype file. The Rscript will write population- and sex-stratified phenotype files. A file listing case and control sample sizes for each population will also be generated.

## 3. Running Time-to-Event GWAS
To perform time-to-event (TTE) GWAS, refer to:
- `ukbb.smk` for example commands to run SAIGE GATE, including:
  - Step 1 and Step 2 for both **onset and progression GWAS**
  - Population and sex-stratified analyses `(F/M/ALL)`
- `WDLs/` for example commands to run SAIGE GATE on AllofUs platform.
- Please note only  GWAS with `N_event>=50` is to be performed. 
- Please note phecode "185" is male only, so no `F` or `ALL` analysis should be run.
### Key Parameters & Recommendations
- **Choosing** `--eventTimeBinSize`: 
  - If the unit of event time is **year**, you could choose `--eventTimeBinSize=1` 
  - If the unit of event time is **month**, `--eventTimeBinSize=1/12` could be a preferred choice.
- **Recommended covariates:**
  - **Onset GWAS:** `gPC + sex + birthyear` + any biobank-specific covariates
  - **Progression GWAS:** `gPC + sex + birthyear +`**`age_of_T0_disease_diagnosis`** + any biobank-specific covariates
  - Sex-stratified analysis: Please remove `sex` from covariates 
- **SAIGE GATE Input Columns:**
  - **Onset GWAS:**
    ```
    --phenoCol=event --eventTimeCol=age
    ```
    - `event`: Phecode status (case/control)
    - `age`: Diagnosis age for cases, or death/last visit age for controls.
  - **Progression GWAS:**
    ```
    --phenoCol=secondEvent --eventTimeCol=secondTime
    ```
    - `secondEvent`: Development of T1 status
    - `secondTime`: Time from T0 to T1 for cases, or death/last visit time from T0 for controls. 
- **Imputation quality:**
  - If genotype is imputed, please use `--is_imputed_data=TRUE --minInfo=0.3` when running SAIGE-GATE step2, so that output file has an `imputationInfo` column and variants with quality score<0.3 will be filtered out. 
  - If genotype is not imputed, `--is_imputed_data=FALSE` by default, and `MissingRate` column will be in the output file.
- **Parallelization:**
  - To speed up SAIGE-GATE step 2, split chromosomes into chunks using `--rangestoIncludeFile` and then merge the output files from each chunk. Ensure that the chromosome format in the range file matches the format in the genotype file (e.g., both should use "chr22" or both should use "22").
- **Chromosome X:**
  - For `ALL` analysis, only autosomal chromosomes should be included in the analysis.
  - For sex-stratified `F/M` analysis, chrX should be added to the analysis.
