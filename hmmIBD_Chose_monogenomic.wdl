version 1.0

# WORKFLOW DEFINITION
workflow hmmIBD_Chose_monogenomic{
  input {
    #String sample_name 
    File vcfFile
    Float? het_thresh
    Boolean onlyGoodSamples = true
    
  }
  call DecideMonogenomic{
    input: 
    vcf = vcfFile,
    onlyGoodSamples = onlyGoodSamples,
    het_thresh = het_thresh
  }

  # Outputs that will be retained when execution is complete
    output {
    File samp_het = DecideMonogenomic.samp_het
    File all_mono_samples = DecideMonogenomic.all_mono_samples
    File all_poly_samples = DecideMonogenomic.all_poly_samples
    File bad_mono_samples = DecideMonogenomic.bad_mono_samples
    File good_mono_samples = DecideMonogenomic.good_mono_samples
    File good_poly_samples = DecideMonogenomic.good_poly_samples
    File good_samples = DecideMonogenomic.good_samples
    File gendata = DecideMonogenomic.gendata
    File allele = DecideMonogenomic.allele
    File hetrate = DecideMonogenomic.hetrate

  }
}
# clean task will run the script clean_vcf.py using 2 required arguments: the VCF File to clean and a list of prefered SNPs
task DecideMonogenomic {
    input {
    # Command parameters
    File vcf
    Boolean onlyGoodSamples = true
    Float? het_thresh
    }
    
    command {
      set -euxo pipefail #if any of the command fails then the entire worfklow fails
      mkdir -p 'output'
      mkdir -p 'seq'
      mkdir -p 'results'

      python /py/vcf2het.py ~{vcf}
      python /py/hetrate.py ~{het_thresh} output/samp_het.txt 
      python /py/vcf2hmm.py ~{vcf} ~{false="" true = "output/good_mono_samples.txt" onlyGoodSamples}
    
    }
    
    runtime {
    docker: "basscigass/hmmibd:1.0.8"
    memory: 64+ " GiB"
    cpu: 16
    disks: "local-disk 100 HDD"
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