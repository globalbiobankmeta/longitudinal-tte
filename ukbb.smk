from os.path import join
import os
import numpy as np
import pandas as pd
import sys


configfile: "config.yaml"

rule all:
    input:
        # Zip t0 and t1 together first, then expand pop, sex, and chr independently
        expand(
            "ukbb/saige_out/{t0}_to_{t1}_{pop}_{sex}_{chr}.index",
            t0=[x[0] for x in zip(config["progress"]["t0"], config["progress"]["t1"])],
            t1=[x[1] for x in zip(config["progress"]["t0"], config["progress"]["t1"])],
            pop=['EUR', 'AFR', 'CSA', 'EAS', 'AMR', 'MID'],
            sex=['F', 'M', 'ALL'],
            chr=[str(c) for c in range(1, 23)] + ["X"]  # Include chrX
        ),
        expand(
            "ukbb/saige_out/{t0}_{pop}_{sex}_{chr}.index",
            t0=config["onset"]["t0"],
            pop=['EUR', 'AFR', 'CSA', 'EAS', 'AMR', 'MID'],
            sex=['F', 'M', 'ALL'],
            chr=[str(c) for c in range(1, 23)] + ["X"]  # Include chrX
        )



rule onset_saige_step1:
    input:
        sparsegrm="SparseGRM",
        grmid="sampleIDs",
        pheno="ukbb/onset/pheno_{t0}_{pop}_{sex}.txt",
        plinkbed="path/plink.bed",
    output:
        "ukbb/saige_out/{t0}_{pop}_{sex}.rda",
        "ukbb/saige_out/{t0}_{pop}_{sex}.varianceRatio.txt",
    log:
        "log/cluster_logs/ukbb.saige1.{t0}_{pop}_{sex}.out",
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
            --phenoCol=event \
            --eventTimeCol=age \
            --covarColList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,sex,birthyear \
            --sampleIDColinphenoFile=s \
            --traitType=survival \
            --outputPrefix=ukbb/saige_out/{wildcards.t0}_{wildcards.pop}_{wildcards.sex} \
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
        rda="ukbb/saige_out/{t0}_{pop}_{sex}.rda",
        varratio="ukbb/saige_out/{t0}_{pop}_{sex}.varianceRatio.txt",
        sample="sample.txt",
    output:
        "ukbb/saige_out/{t0}_{pop}_{sex}_{chr}",
        "ukbb/saige_out/{t0}_{pop}_{sex}_{chr}.index",
    log:
        "log/cluster_logs/ukbb.saige2.{t0}_{pop}_{sex}_{chr}.out",
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
            --SAIGEOutputFile=ukbb/saige_out/{wildcards.t0}_{wildcards.pop}_{wildcards.sex}_{wildcards.chr} \
            > {log} 2>&1
        """

        
rule progress_saige_step1:
    input:
        sparsegrm="SparseGRM",
        grmid="sampleIDs",
        pheno="ukbb/progress/pheno_{t0}_to_{t1}_{pop}_{sex}.txt",
        plinkbed="path/plink.bed",
    output:
        "ukbb/saige_out/{t0}_to_{t1}_{pop}_{sex}.rda",
        "ukbb/saige_out/{t0}_to_{t1}_{pop}_{sex}.varianceRatio.txt",
    log:
        "log/cluster_logs/ukbb.saige1.{t0}_{t1}_{pop}_{sex}.out",
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
            --outputPrefix=ukbb/saige_out/{wildcards.t0}_to_{wildcards.t1}_{wildcards.pop}_{wildcards.sex} \
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
        rda="ukbb/saige_out/{t0}_to_{t1}_{pop}_{sex}.rda",
        varratio="ukbb/saige_out/{t0}_to_{t1}_{pop}_{sex}.varianceRatio.txt",
        sample="sample.txt",
    output:
        "ukbb/saige_out/{t0}_to_{t1}_{pop}_{sex}_{chr}",
        "ukbb/saige_out/{t0}_to_{t1}_{pop}_{sex}_{chr}.index",
    log:
        "log/cluster_logs/ukbb.saige2.{t0}_{t1}_{pop}_{sex}_{chr}.out",
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
            --SAIGEOutputFile=ukbb/saige_out/{wildcards.t0}_to_{wildcards.t1}_{wildcards.chr}_{wildcards.pop}_{wildcards.sex} \
            > {log} 2>&1
        """


