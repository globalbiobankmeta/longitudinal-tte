from os.path import join
import os
import numpy as np
import pandas as pd
import sys


configfile: "config.yaml"

rule all:
    input:
        [
            f"ukbb/saige_out/{t0}_to_{t1}_{chr}.index"
            for t0, t1 in zip(config["pheno"]["t0"], config["pheno"]["t1"])
            for chr in np.arange(1, 23, 1)
        ]
        
rule saige_step1:
    input:
        sparsegrm="SparseGRM",
        grmid="sampleIDs",
        pheno="ukbb/pheno_{t0}_to_{t1}.txt",
        plinkbed="path/plink.bed",
    output:
        "ukbb/saige_out/{t0}_to_{t1}.rda",
        "ukbb/saige_out/{t0}_to_{t1}.varianceRatio.txt",
    log:
        "log/cluster_logs/ukbb.saige1.{t0}_{t1}.out",
    params:
        plinkpath="path/plink",
    shell:
        """
        mkdir -p ukbb/saige_out/

        Rscript SAIGE/extdata/step1_fitNULLGLMM.R \
            --sparseGRMFile={input.sparsegrm} \
            --sparseGRMSampleIDFile={input.grmid} \
            --useSparseGRMtoFitNULL=TRUE \
            --plinkFile={params.plinkpath} \
            --phenoFile={input.pheno} \
            --phenoCol=secondEvent \
            --eventTimeCol=secondTime \
            --covarColList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,sex,diagAge,birthyear \
            --sampleIDColinphenoFile=s \
            --traitType=survival \
            --outputPrefix=ukbb/saige_out/{wildcards.t0}_to_{wildcards.t1} \
            --nThreads=1 \
            --traceCVcutoff=0.0025 \
            --ratioCVcutoff=0.001 \
            --pcgforUhatforSurvAnalysis=TRUE \
            --eventTimeBinSize=1 \
            > {log} 2>&1
        """

rule saige_step2:
    input:
        bgen="ukb_imp_chr{chr}_v3.bgen",
        rda="ukbb/saige_out/{t0}_to_{t1}.rda",
        varratio="ukbb/saige_out/{t0}_to_{t1}.varianceRatio.txt",
        sample="sample.txt",
    output:
        "ukbb/saige_out/{t0}_to_{t1}_{chr}",
        "ukbb/saige_out/{t0}_to_{t1}_{chr}.index",
    log:
        "log/cluster_logs/ukbb.saige2.{t0}_{t1}_{chr}.out",
    resources:
        mem_gb=10,
        time_min=480,
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
            --SAIGEOutputFile=ukbb/saige_out/{wildcards.t0}_to_{wildcards.t1}_{wildcards.chr} \
            > {log} 2>&1
        """