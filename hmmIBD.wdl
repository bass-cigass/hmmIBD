version 1.0

## Input:
## Outputs :

## 
# WORKFLOW DEFINITION
workflow hmmIBD{
  input {
    #String sample_name 
    File vcfFile
    File? freqData
    String output_pfx = "hmm"
    
  }

  call prepareData{
    input: 
    vcf = vcfFile
  }

  call run_hmmIBD {
    input:
    data = prepareData.gendata,
    freqData = freqData,
    output_pfx = output_pfx
    }
  
  # Outputs that will be retained when execution is complete
    output {
    File out1 = run_hmmIBD.output1
    File out2 = run_hmmIBD.output2

  }
}
# Align a pair of FASTQs and output a bam file
task run_hmmIBD {
    input {
    # Command parameters
    File data
    File? freqData
    String output_pfx 
    }
    command {
    set -euxo pipefail #if any of the command fails then the entire worfklow fails
    ./hmmIBD -i ~{data} ~{'-f '+freqData} -o ~{output_pfx}
    }
    
    runtime {
    docker: "basscigass/hmmibd:1.0.8"
    memory: 10+ " GiB"
    disks: "local-disk 50 HDD"
    cpu_cores: 8
    preemptible: 0
    }
    
    output {
    File output1 = output_pfx +".hmm_fract.txt"
    File output2 = output_pfx +".hmm.txt"
    }
}
task prepareData {
    input {
    # Command parameters
    File vcf
    }
    command {
      set -euxo pipefail #if any of the command fails then the entire worfklow fails
      mkdir 'output'
      mkdir 'seq'
      mkdir 'results'
      python /py/vcf2het.py ~{vcf}
      python /py/vcf2hmm.py ~{vcf}
      python /py/hetrate.py output/samp_het.txt
    }
    
    runtime {
    docker: "basscigass/hmmibd:1.0.8"
    memory: 32+ " GiB"
    disks: "local-disk 50 HDD"
    preemptible: 0
    }
    
    output {
    File samp_het = "output/samp_het.txt"
    File all_mono_samples = "output/all_mono_samples.txt"
    File all_poly_samples = "output/all_poly_samples.txt"
    File bad_mono_samples = "output/bad_mono_samples.txt"
    File good_mono_samples = "output/good_mono_samples.txt"
    File good_poly_samples = "output/good_poly_samples.txt"
    File good_samples = "output/good_samples.txt"
    File gendata = "seq/seq.txt"
    File allele = "seq/allele.txt"
    }
}
