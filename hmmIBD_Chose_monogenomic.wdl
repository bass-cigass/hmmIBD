version 1.0

# WORKFLOW DEFINITION
workflow hmmIBD_Chose_monogenomic{
  input {
    #String sample_name 
    File vcfFile
    Float? het_thresh
    Boolean onlyGoodSamples = true
    String year = 'all'
    
  }
  call DecideMonogenomic{
    input: 
    vcf = vcfFile,
    onlyGoodSamples = onlyGoodSamples,
    het_thresh = het_thresh,
    target_year = year
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
    File genotype_data = DecideMonogenomic.genotype_data
    File allele = DecideMonogenomic.allele
    File freq = DecideMonogenomic.freq
    File seq = DecideMonogenomic.seq
    File hetrate = DecideMonogenomic.hetrate
    File thinnedSites = DecideMonogenomic.thinnedSite

  }
}
# clean task will run the script clean_vcf.py using 2 required arguments: the VCF File to clean and a list of prefered SNPs
task DecideMonogenomic {
    input {
    # Command parameters
    File vcf
    Boolean onlyGoodSamples = true
    String target_year
    Float? het_thresh
    String outdir = ' -s output/~{target_year}_good_mono_samples.txt'
    
    }
    
    command {
      set -euxo pipefail #if any of the command fails then the entire worfklow fails
      mkdir -p 'output'
      mkdir -p 'seq'
      mkdir -p 'results'
      mkdir -p 'hmmInput'

      python /py/vcf2het.py ~{vcf} ~{target_year}
      python /py/hetrate.py ~{target_year} ~{het_thresh}
      #uncomment and debug the line below for the case of all sample
      #python /py/vcf2hmm.py ~{vcf} 'seq/out' ~{false='' true='output dir' onlyGoodSamples}
      python /py/vcf2hmm_st.py ~{vcf} 'seq/out' ~{outdir} 
      python /py/thin_sites.py "seq/out_freq.txt" "seq/thinned_Site.txt"
      python /py/thin_seq.py "seq/thinned_Site.txt" "seq/out_seq.txt" "hmmInput/thin_seq.txt"
    
    }
    
    runtime {
    docker: "basscigass/hmmibd:1.1.9"
    memory: 64+ " GiB"
    cpu: 16
    disks: "local-disk 100 HDD"
    preemptible: 0
    }
    
    output {
    File samp_het = "output/~{target_year}_samp_het.txt"
    File all_mono_samples = "output/~{target_year}_mono_samples.txt"
    File all_poly_samples = "output/~{target_year}_poly_samples.txt"
    File bad_mono_samples = "output/~{target_year}_bad_mono_samples.txt"
    File good_mono_samples = "output/~{target_year}_good_mono_samples.txt"
    File good_poly_samples = "output/~{target_year}_good_poly_samples.txt"
    File good_samples = "output/~{target_year}_good_samples.txt"
    File genotype_data = "hmmInput/thin_seq.txt"
    File seq = "seq/out_seq.txt"
    File allele = "seq/out_allele.txt"
    File freq = "seq/out_freq.txt"
    File hetrate = "results/~{target_year}_hetrate.pdf"
    File thinnedSite = "seq/thinned_Site.txt"
    }
}