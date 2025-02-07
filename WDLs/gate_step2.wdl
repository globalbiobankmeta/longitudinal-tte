version 1.0

task step2 {
    input {
        File bgenfile
        File bgenfileindex
        File samplefile
        File modelfile
        File varianceratiofile
        Int chr
    }

    command <<<
        export MKL_NUM_THREADS=1;export MKL_DYNAMIC=false; export OMP_NUM_THREADS=1; export OMP_DYNAMIC=false;
        
        step2_SPAtests.R \
            --bgenFile=~{bgenfile} \
            --bgenFileIndex=~{bgenfileindex} \
            --sampleFile=~{samplefile} \
            --GMMATmodelFile=~{modelfile} \
            --varianceRatioFile=~{varianceratiofile} \
            --minMAC=20 \
            --LOCO=FALSE \
            --AlleleOrder=ref-first \
            --SAIGEOutputFile="~{sub(basename(modelfile), "\\.rda$", "")}_~{chr}.GATE_results.txt"
    >>>
    
    output {
        File results = "~{sub(basename(modelfile), "\\.rda$", "")}_~{chr}.GATE_results.txt"
    }
    
    runtime {
        docker: '""" + os.environ["ARTIFACT_REGISTRY_DOCKER_REPO"] + """/wzhou88/saige:1.4.0'
        cpu: "1"
        memory: "4 GB"
        disks: "local-disk " + ceil(size(bgenfile, "G") + 35) + " HDD"
    }
}

workflow gate_step2 {
    input {
        String bgenbase
        File samplefile
        File modelfile
        File varianceratiofile
    }
    
    # Create an array with indices from 1 to 22
    Array[Int] chromosome_indices = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]

    # Scatter block to run the task for each chromosome/input file
    scatter (i in chromosome_indices) {
        call step2 {
            input:
                bgenfile = bgenbase + i + ".bgen",
                bgenfileindex = bgenbase + i + ".bgen.bgi",
                samplefile = samplefile,
                modelfile = modelfile,
                varianceratiofile = varianceratiofile,
                chr = i
        }
    }
}