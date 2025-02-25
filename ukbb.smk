from os.path import join
import os
import numpy as np
import pandas as pd
import sys


configfile: "config.yaml"

"""
Rules:
Please note phecode 185 is male only, only analysis on male samples should be run.
Please note only GWAS with N_case>=50 is to be run. Please adjust to your biobank.

1. onset_saige_step1: 
- The first step of SAIGE is to fit a null model using the sparse GRM. This is done separately for each population and F/M/ALL. 
- The output is a .rda file containing the fitted model and a .varianceRatio.txt file containing the variance ratio estimates.
- The input phenotype file should contain the columns "event" and "age" for the onset survival analysis. "event" is the phecode status, and "age" is the diagonosis age for cases, or current/death age for controls.
- The covariates used in the model are gPC, birthyear, and sex if analysis is on ALL.

2. onset_saige_step2
- The second step of SAIGE is to perform the association test using the fitted null model from step 1. This is done separately for each chromosome.
- The input is the .rda file and .varianceRatio.txt file from step 1, as well as the genotype data in BGEN format. VCF or plink format is also supported. See SAIGE documentation for details.
- Specify --is_imputed_data=TRUE and --minInfo=0.3 if genotype is imputed.
- Option to split into chunks with --rangestoIncludeFile.

3. progress_saige_step1
- The input phenotype file should contain the columns "secondEvent" and "secondTime" for the progress survival analysis. "secondEvent" is the T0 to T1 status, and "secondTime" is the progression time for cases, or current/death age from T0 diagnosis for controls.
- The covariates used in the model are gPC, diagAge, birthyear, and sex if analysis is on ALL.

4. progress_saige_step2
- The input is the .rda file and .varianceRatio.txt file from step 1, as well as the genotype data in BGEN format. VCF or plink format is also supported. See SAIGE documentation for details.
- Specify --is_imputed_data=TRUE and --minInfo=0.3 if genotype is imputed.
- Option to split into chunks with --rangestoIncludeFile.
 
"""

rule all:
    input:
        # Zip t0 and t1 together first, then expand population, sex, and chr independently
        expand(
            "saige_out_progress/{t0_t1[0]}_to_{t0_t1[1]}_{population}_{sex}_{chr}.txt.index",
            t0_t1=list(zip(config["progress"]["t0"], config["progress"]["t1"])),
            population=['EUR', 'AFR', 'CSA', 'EAS', 'AMR', 'MID'],
            sex=['F', 'M', 'ALL'],
            chr=[str(c) for c in range(1, 23)] + ["X"]  # Include chrX
        ),
        expand(
            "saige_out_onset/{t0}_{population}_{sex}_{chr}.txt.index",
            t0=config["onset"]["t0"],
            population=['EUR', 'AFR', 'CSA', 'EAS', 'AMR', 'MID'],
            sex=['F', 'M', 'ALL'],
            chr=[str(c) for c in range(1, 23)] + ["X"]  # Include chrX
        )



rule onset_saige_step1:
    input:
        sparsegrm="SparseGRM",
        grmid="sampleIDs",
        pheno="onset/pheno_{t0}_{population}_{sex}.txt",
        plinkbed="path/plink.bed",
    output:
        "saige_out_onset/{t0}_{population}_{sex}.rda",
        "saige_out_onset/{t0}_{population}_{sex}.varianceRatio.txt",
    log:
        "log/cluster_logs/ukbb.saige1.{t0}_{population}_{sex}.out",
    params:
        plinkpath="path/plink",
        covar=lambda wildcards: "PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,sex,birthyear" if wildcards.sex == "ALL" else "PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,birthyear",
    shell:
        """
        mkdir -p saige_out_onset/

        Rscript SAIGE/extdata/step1_fitNULLGLMM.R \
            --sparseGRMFile={input.sparsegrm} \
            --sparseGRMSampleIDFile={input.grmid} \
            --useSparseGRMtoFitNULL=TRUE \
            --plinkFile={params.plinkpath} \
            --phenoFile={input.pheno} \
            --phenoCol=event \
            --eventTimeCol=age \
            --covarColList={params.covar} \
            --sampleIDColinphenoFile=s \
            --traitType=survival \
            --outputPrefix=saige_out_onset/{wildcards.t0}_{wildcards.population}_{wildcards.sex} \
            --nThreads=1 \
            --traceCVcutoff=0.0025 \
            --ratioCVcutoff=0.001 \
            --pcgforUhatforSurvAnalysis=TRUE \
            --eventTimeBinSize=1 \
            > {log} 2>&1
        """

rule onset_saige_step2:
    input:
        bgen="ukb_imp_chr{chr}_v3.bgen",
        rda="saige_out_onset/{t0}_{population}_{sex}.rda",
        varratio="saige_out_onset/{t0}_{population}_{sex}.varianceRatio.txt",
        sample="sample.txt",
    output:
        "saige_out_onset/{t0}_{population}_{sex}_{chr}.txt",
        "saige_out_onset/{t0}_{population}_{sex}_{chr}.txt.index",
    log:
        "log/cluster_logs/ukbb.saige2.{t0}_{population}_{sex}_{chr}.txt.out",
    shell:
        """
        Rscript SAIGE/extdata/step2_SPAtests.R \
            --bgenFile={input.bgen} \
            --bgenFileIndex={input.bgen}.bgi \
            --sampleFile={input.sample} \
            --GMMATmodelFile={input.rda} \
            --varianceRatioFile={input.varratio} \
            --minMAC=20 \
            --LOCO=FALSE \
            --AlleleOrder=ref-first \
            --SAIGEOutputFile={output[0]} \
            --is_imputed_data=TRUE \
            --minInfo=0.3 \
            > {log} 2>&1
        """

        
rule progress_saige_step1:
    input:
        sparsegrm="SparseGRM",
        grmid="sampleIDs",
        pheno="progress/pheno_{t0}_to_{t1}_{population}_{sex}.txt",
        plinkbed="path/plink.bed",
    output:
        "saige_out_progress/{t0}_to_{t1}_{population}_{sex}.rda",
        "saige_out_progress/{t0}_to_{t1}_{population}_{sex}.varianceRatio.txt",
    log:
        "log/cluster_logs/ukbb.saige1.{t0}_{t1}_{population}_{sex}.out",
    params:
        plinkpath="path/plink",
        covar=lambda wildcards: "PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,sex,diagAge,birthyear" if wildcards.sex == "ALL" else "PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,diagAge,birthyear",
    shell:
        """
        mkdir -p saige_out_progress/

        Rscript SAIGE/extdata/step1_fitNULLGLMM.R \
            --sparseGRMFile={input.sparsegrm} \
            --sparseGRMSampleIDFile={input.grmid} \
            --useSparseGRMtoFitNULL=TRUE \
            --plinkFile={params.plinkpath} \
            --phenoFile={input.pheno} \
            --phenoCol=secondEvent \
            --eventTimeCol=secondTime \
            --covarColList={params.covar} \
            --sampleIDColinphenoFile=s \
            --traitType=survival \
            --outputPrefix=saige_out_progress/{wildcards.t0}_to_{wildcards.t1}_{wildcards.population}_{wildcards.sex} \
            --nThreads=1 \
            --traceCVcutoff=0.0025 \
            --ratioCVcutoff=0.001 \
            --pcgforUhatforSurvAnalysis=TRUE \
            --eventTimeBinSize=1 \
            > {log} 2>&1
        """

rule progress_saige_step2:
    input:
        bgen="ukb_imp_chr{chr}_v3.bgen",
        rda="saige_out_progress/{t0}_to_{t1}_{population}_{sex}.rda",
        varratio="saige_out_progress/{t0}_to_{t1}_{population}_{sex}.varianceRatio.txt",
        sample="sample.txt",
    output:
        "saige_out_progress/{t0}_to_{t1}_{population}_{sex}_{chr}.txt",
        "saige_out_progress/{t0}_to_{t1}_{population}_{sex}_{chr}.txt.index",
    log:
        "log/cluster_logs/ukbb.saige2.{t0}_{t1}_{population}_{sex}_{chr}.txt.out",
    shell:
        """
        Rscript SAIGE/extdata/step2_SPAtests.R \
            --bgenFile={input.bgen} \
            --bgenFileIndex={input.bgen}.bgi \
            --sampleFile={input.sample} \
            --GMMATmodelFile={input.rda} \
            --varianceRatioFile={input.varratio} \
            --minMAC=20 \
            --LOCO=FALSE \
            --AlleleOrder=ref-first \
            --SAIGEOutputFile={output[0]} \
            --is_imputed_data=TRUE \
            --minInfo=0.3 \
            > {log} 2>&1
        """


