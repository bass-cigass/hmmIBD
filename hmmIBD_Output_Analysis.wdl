version 1.0

## 
# WORKFLOW DEFINITION
workflow hmmIBD{
  input {
    File hmm_File
    File hmm_fract_File
    
  }

  call run_Pileup{
    input: 
    hmm_File = hmm_File,
    hmm_fract_File = hmm_fract_File,
  }
  
  # Outputs that will be retained when execution is complete
    output {
    File result = run_Pileup.result

  }
}
# Align a pair of FASTQs and output a bam file
task run_Pileup {
    input {
    File hmm_File
    File hmm_fract_File
     
    }
    command {
    set -euxo pipefail #if any of the command fails then the entire worfklow fails
    python /py/pileup.py ~{hmm_File} ~{hmm_fract_File} "output/result_plot.pdf"
    }
    
    runtime {
    docker: "basscigass/hmmibd:1.0.8"
    memory: 8+ " GiB"
    disks: "local-disk 50 HDD"
    cpu: 4
    preemptible: 0
    }
    
    output {
    File result = "output/result_plot.pdf"
    
    }
}