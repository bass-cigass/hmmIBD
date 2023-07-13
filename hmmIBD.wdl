version 1.0

## Input:
## Outputs :

## 
# WORKFLOW DEFINITION
workflow hmmIBD{
  input {
    #String sample_name 
    File vcfFile
    File? snplist_preferred
    File? freqData
    String output_pfx = "hmm"
    
  }

  call prepareData{
    input: 
    vcf = vcfFile,
    snplist = snplist_preferred
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
    File output1 = output_pfx +".hmm_fract.txt"
    File output2 = output_pfx +".hmm.txt"
    }
}
task prepareData {
    input {
    # Command parameters
    File vcf
    File? snplist
    }
    String filename = 'output/cleaned_vcf.vcf'
    Boolean hasSnp = defined(snplist)
    command {
      set -euxo pipefail #if any of the command fails then the entire worfklow fails
      mkdir 'output'
      mkdir 'seq'
      mkdir 'results'

      if ~{hasSnp} ; then
        python /py/clean_vcf.py ~{vcf} ~{filename} ~{snplist}
        python /py/vcf2het.py ~{filename}
        python /py/hetrate.py output/samp_het.txt
        python /py/vcf2hmm.py ~{filename} output/good_mono_samples.txt
      else
        python /py/vcf2het.py ~{vcf}
        python /py/hetrate.py output/samp_het.txt
        python /py/vcf2hmm.py ~{vcf} output/good_mono_samples.txt
        touch ~{filename}
      fi
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
    File cleaned_vcf = "output/cleaned_vcf.vcf"
    }
}
