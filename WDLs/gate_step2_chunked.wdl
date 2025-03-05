version 1.0

task step2 {
    input {
        File bgenfile
        File bgenfileindex
        File samplefile
        File modelfile
        File varianceratiofile
        File rangefile
        Int chr
        Int chunk
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
            --rangestoIncludeFile=~{rangefile} \
            --SAIGEOutputFile="~{sub(basename(modelfile), "\\.rda$", "")}_~{chr}_~{chunk}.GATE_results.txt"
    >>>
    
    output {
        File results = "~{sub(basename(modelfile), "\\.rda$", "")}_~{chr}_~{chunk}.GATE_results.txt"
    }
    
    runtime {
        docker: '""" + os.environ["ARTIFACT_REGISTRY_DOCKER_REPO"] + """/wzhou88/saige:1.4.4'
        cpu: "1"
        memory: "4 GB"
        disks: "local-disk " + ceil(size(bgenfile, "G") + 35) + " HDD"
    }
}

workflow gate_step2 {
    input {
        String bgenbase
        String rangebase
        File samplefile
        File modelfile
        File varianceratiofile
    }
    
    Array[Int] chromosome_indices = [21,22]
    Array[Int] chunk_indices = [1,2,3,4]
    
    # Create the array of all pairs
    Array[Pair[Int, Int]] runs = cross(chromosome_indices, chunk_indices)

    # Scatter block to run the task for each chromosome/input file
    scatter (run in runs) {
        call step2 {
            input:
                bgenfile = bgenbase + run.left + ".bgen",
                bgenfileindex = bgenbase + run.left + ".bgen.bgi",
                samplefile = samplefile,
                modelfile = modelfile,
                varianceratiofile = varianceratiofile,
                rangefile = rangebase + "chr" + run.left + "_chunk" + run.right + ".txt",
                chr = run.left,
                chunk = run.right
        }
    }
}