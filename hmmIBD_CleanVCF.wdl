version 1.0

## 
# WORKFLOW DEFINITION
workflow hmmIBD_CleanVCF{
  input {
    File vcfFile
     File snplist_preferred    
  }
 #Task definition
  call clean{
    input: 
    vcf = vcfFile,
    snplist = snplist_preferred
  }
  
  # Outputs that will be retained when execution is complete
    output {
    File result = clean.cleaned_vcf

  }
}
# clean task will run the script clean_vcf.py using 2 required arguments: the VCF File to clean and a list of prefered SNPs
task clean {
    input {
    # Command parameters
    File vcf
    File snplist
    }   

    command {
      set -euxo pipefail #if any of the command fails then the entire worfklow fails
      mkdir -p 'output'

      python /py/clean_vcf.py ~{vcf} 'output/cleaned_vcf.vcf' ~{snplist}
      
    }
    
    runtime {
    docker: "basscigass/hmmibd:1.0.8"
    memory: 32+ " GiB"
    cpu: 8
    disks: "local-disk 50 HDD"
    preemptible: 0
    }
    
    output {
    File cleaned_vcf = "output/cleaned_vcf.vcf"
    }
}
