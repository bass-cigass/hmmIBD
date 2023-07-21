version 1.0

## Input:vcfFile, output_pfx

## Outputs :
## 
# WORKFLOW DEFINITION
workflow hmmIBD{
  input {
    #String sample_name 
    File genotype_data
    File? freqData
    String output_pfx = "hmm" 
  }
  #Calling task remove_sweeps 
 call remove_sweeps {
  input:
    infile = genotype_data,
    output_file = "sweeps_removed_seq.txt"
  }
  if (defined(freqData))  { 
      call remove_sweeps as sweeps_removed_freq {
    input:
      infile = freqData,
      output_file = "sweeps_removed_freq.txt"

    }
  }
#Calling task run_hmmIBD 
  call run_hmmIBD {
    input:
    data = remove_sweeps.sweeps_removed,
    freqData = if defined(freqData) then sweeps_removed_freq.sweeps_removed else freqData,
    output_pfx = output_pfx
    }
  
  # Outputs that will be retained when execution is complete
    output {
    File hmm_fract = run_hmmIBD.hmm_fract
    File hmm_file = run_hmmIBD.hmm_file

  }
}

######## Task remove_sweeps will remove the known resistances from the sequences list
task remove_sweeps {
    input {
    # Command parameters
    File? infile
    String output_file
    }
    command {
    set -euxo pipefail #if any of the command fails then the entire worfklow fails
    mkdir -p 'seq'

    python /py/remove_sweeps.py ~{infile} "seq/"+output_file
    }
    #runtime configuration
    runtime {
    docker: "basscigass/hmmibd:1.0.8"
    memory: 16+ " GiB"
    disks: "local-disk 50 HDD"
    cpu: 4
    preemptible: 0
    }
    
    output {
    File sweeps_removed = "seq/"+output_file
    }
}

#### Task run_hmmIBD will run the hmmIBD scipt with 2 required arguments: the vcf file (data) and the prefix of the output files (output_pfx)
# One optional argument can be provided: freqData
task run_hmmIBD {
    input {
    # Command parameters
    File data
    File? freqData
    String output_pfx 
    }
    command {
    set -euxo pipefail #if any of the command fails then the entire worfklow fails
    hmmIBD -i ~{data} ~{'-f '+freqData} -o ~{output_pfx}
    }
    
    runtime {
    docker: "basscigass/hmmibd:1.0.8"
    memory: 16+ " GiB"
    disks: "local-disk 50 HDD"
    cpu: 4
    preemptible: 0
    }
    
    output {
    File hmm_fract = output_pfx +".hmm_fract.txt"
    File hmm_file = output_pfx +".hmm.txt"
    }
}

