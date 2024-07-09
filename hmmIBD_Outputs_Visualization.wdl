version 1.0

## 
# WORKFLOW DEFINITION
workflow hmmIBD_Outputs_Visualization{
  input {
    File hmm_File
    File hmm_fract_File
    File locus
    
  }

  call run_Pileup{
    input: 
    hmm_File = hmm_File,
    hmm_fract_File = hmm_fract_File,
    locus_gene = locus
  }
  call plot_IBD{
    input: 
    hmm_fract_File = hmm_fract_File,
  }
  
  # Outputs that will be retained when execution is complete
    output {
    File pileupPlot = run_Pileup.result
    File IBDplot = plot_IBD.plot

  }
}
# Align a pair of FASTQs and output a bam file
task run_Pileup {
    input {
    File hmm_File
    File hmm_fract_File
    File locus_gene
     
    }
    command {
    set -euxo pipefail #if any of the command fails then the entire worfklow fails
    mkdir -p 'output'
    python /py/pileup_st.py ~{hmm_File} ~{hmm_fract_File} "output/result_plot.pdf" ~{locus_gene}
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
task plot_IBD {
    input {
    File hmm_fract_File
     
    }
    command {
    set -euxo pipefail #if any of the command fails then the entire worfklow fails
    mkdir - 'results'
    python /py/plot_ibd.py ~{hmm_fract_File} 
    }
    
    runtime {
    docker: "basscigass/hmmibd:1.0.8"
    memory: 8+ " GiB"
    disks: "local-disk 50 HDD"
    cpu: 4
    preemptible: 0
    }
    
    output {
    File plot = "results/plot_ibd.pdf"
    
    }
}