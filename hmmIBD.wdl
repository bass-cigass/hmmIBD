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
    Boolean onlyGoodSamples = true
    
  }
  
  call prepareData{
    input: 
    vcf = vcfFile,
    onlyGoodSamples = onlyGoodSamples
  }

  call run_hmmIBD {
    input:
    data = prepareData.gendata,
    freqData = freqData,
    output_pfx = output_pfx
    }
  
  # Outputs that will be retained when execution is complete
    output {
    File hmm_fract = run_hmmIBD.hmm_fract
    File hmm_file = run_hmmIBD.hmm_file
    File samp_het = prepareData.samp_het
    File all_mono_samples = prepareData.all_mono_samples
    File all_poly_samples = prepareData.all_poly_samples
    File bad_mono_samples = prepareData.bad_mono_samples
    File good_mono_samples = prepareData.good_mono_samples
    File good_poly_samples = prepareData.good_poly_samples
    File good_samples = prepareData.good_samples
    File gendata = prepareData.gendata
    File allele = prepareData.allele
    File hetrate = prepareData.hetrate

  }
}
# Task run_hmmIBD will run the hmmIBD scipt with 2 required arguments: the vcf file (data) and the prefix of the output files (output_pfx)
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
task prepareData {
    input {
    # Command parameters
    File vcf
    Boolean onlyGoodSamples = true
    }
    
    command {
      set -euxo pipefail #if any of the command fails then the entire worfklow fails
      mkdir -p 'output'
      mkdir -p 'seq'
      mkdir -p 'results'

      python /py/vcf2het.py ~{vcf}
      python /py/hetrate.py output/samp_het.txt
      python /py/vcf2hmm.py ~{vcf} ~{false="" true = "output/good_mono_samples.txt" onlyGoodSamples}
    
    }
    
    runtime {
    docker: "basscigass/hmmibd:1.0.8"
    memory: 32+ " GiB"
    cpu: 8
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
    File hetrate = "results/hetrate.pdf"
    }
}
