version 1.0

task step2 {
    input {
        File bgenfile
        File bgenfileindex
        File samplefile
        File modelfile
        File varianceratiofile
        Array[File] rangefiles
    }

    command <<<
        python3 <<EOF
        import os
        import subprocess
        import time
        processes = set()
        
        cmd_prefix = '''
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
        '''
        
        processes = set()

        for file in '~{sep=" " rangefiles}'.split(' '):
            cmd = cmd_prefix + '--rangestoIncludeFile=' + file
            cmd += ' --SAIGEOutputFile=' + '~{sub(basename(modelfile), "\\.rda$", "")}_' + os.path.basename(file) + '.GATE_results.txt'
            logfilepath = 'GATE_log_' + '~{sub(basename(modelfile), "\\.rda$", "")}_' + os.path.basename(file) + '.txt'

            logfile = open(logfilepath, 'w')

            # Start all processes
            p = subprocess.Popen(cmd, shell=True, stdout=logfile)
            processes.add(p)

            print(time.strftime("%Y/%m/%d %H:%M:%S") + ' ' + str(len(processes)) + ' processes started', flush=True)

        # Now, wait for all processes to finish
        while processes:
            time.sleep(60)
            for p in list(processes):  # Iterate over a copy to allow modification
                p_poll = p.poll()
                if p_poll is not None:
                    if p_poll > 0:
                        raise Exception('subprocess returned ' + str(p_poll))
                    processes.remove(p)  # Remove completed process

            print(time.strftime("%Y/%m/%d %H:%M:%S") + ' ' + str(len(processes)) + ' processes still running', flush=True)

        EOF
    >>>
    
    output {
        Array[File] out = glob("*.GATE_results.txt")
        Array[File] logs = glob("GATE_log_*.txt")
    }
    
    runtime {
        docker: '""" + os.environ["ARTIFACT_REGISTRY_DOCKER_REPO"] + """/wzhou88/saige:1.4.4'
        cpu: "4"
        memory: "16 GB"
        disks: "local-disk " + ceil(size(bgenfile, "G") + 35) + " HDD"
        zones: "us-central1-b"
        preemptible: 4
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
    
    # Create an array with indices from 1 to 22
    Array[Int] chromosomes = [22]

    # Scatter block to run the task for each chromosome/input file
    scatter (chrom in chromosomes) {
        Array[String] rangefiles = [rangebase + "chr" + chrom + "_chunk1.txt",
                                    rangebase + "chr" + chrom + "_chunk2.txt",
                                    rangebase + "chr" + chrom + "_chunk3.txt",
                                    rangebase + "chr" + chrom + "_chunk4.txt"]
        call step2 {
            input:
                bgenfile = bgenbase + chrom + ".bgen",
                bgenfileindex = bgenbase + chrom + ".bgen.bgi",
                samplefile = samplefile,
                modelfile = modelfile,
                varianceratiofile = varianceratiofile,
                rangefiles = rangefiles
        }
    }
}