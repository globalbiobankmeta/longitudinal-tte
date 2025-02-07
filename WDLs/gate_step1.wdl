version 1.0

task step1 {
    input {
        File sparsegrmfile
        File sparsegrmsamplefile
        File bedfile
        File bimfile
        File famfile
        File phenofile
        String covariates
        String pheno
        Int cpu
        Int mem
        Float traceCVcutoff
        Float ratioCVcutoff
        String pcgforUhatforSurvAnalysis
    }

    command <<<
        step1_fitNULLGLMM.R \
            --sparseGRMFile=~{sparsegrmfile} \
            --sparseGRMSampleIDFile=~{sparsegrmsamplefile} \
            --useSparseGRMtoFitNULL=TRUE \
            --plinkFile=~{sub(bedfile, "\\.bed$", "")} \
            --phenoFile=~{phenofile} \
            --phenoCol=secondEvent \
            --eventTimeCol=TTE \
            --eventTimeBinSize=0.083 \
            --covarColList=~{covariates} \
            --sampleIDColinphenoFile=s \
            --traitType=survival \
            --outputPrefix=~{pheno} \
            --nThreads=1 \
            --traceCVcutoff=0.0025 \
            --ratioCVcutoff=0.001 \
            --pcgforUhatforSurvAnalysis=TRUE
    >>>
    
    output {
        File modelfile = pheno + ".rda"
        File varianceratio = pheno + ".varianceRatio.txt"
    }
    
    runtime {
        docker: '""" + os.environ["ARTIFACT_REGISTRY_DOCKER_REPO"] + """/wzhou88/saige:1.4.0'
        cpu: "~{cpu}"
        memory: "~{mem} GB"
        disks: "local-disk 120 HDD"
    }
}

workflow gate_step1 {
    input {
        File sparsegrmfile
        File sparsegrmsamplefile
        File bedfile
        File bimfile
        File famfile
        File phenofile
        String covariates
        String pheno
        Int cpu
        Int mem
        Float traceCVcutoff
        Float ratioCVcutoff
        String pcgforUhatforSurvAnalysis
    }
        
    call step1 {
        input:
            sparsegrmfile = sparsegrmfile,
            sparsegrmsamplefile = sparsegrmsamplefile,
            bedfile = bedfile,
            bimfile = bimfile,
            famfile = famfile,
            phenofile = phenofile,
            covariates = covariates,
            pheno = pheno,
            cpu = cpu,
            mem = mem,
            traceCVcutoff = traceCVcutoff,
            ratioCVcutoff = ratioCVcutoff,
            pcgforUhatforSurvAnalysis = pcgforUhatforSurvAnalysis
    }
}