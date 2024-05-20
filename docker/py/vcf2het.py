#!/usr/local/bin/python3

from collections import Counter
import gzip
import re
import sys

def main() :
  kill_indel = True
  # Mininum genotyping call rate to accept variant
  #min_call = 0.80
  # Option: minimum number of alt allele copies to keep variant
  #  (1 = kill monomorphic, 2 = kill singletons)
  #min_copy = 0
  min_depth = 5
  
  if len(sys.argv) < 2 : sys.exit('Usage: vcf2het.py <input vcf file name>')
  infile = sys.argv[1]
  outfile = 'output/' + 'samp_het.txt'
    
  if re.search(r'\.gz$', infile) :
    inf = gzip.open(infile, 'rt')
  else :
    inf = open(infile, 'r')

  if len(sys.argv) == 3 :
    if sys.argv[2]=='False':
      kill_indel = False
  print('Kiling Indel :', kill_indel)    
  samples = []
  nindel = 0
  
  callcount_persamp = Counter()
  hetcount_persamp = Counter()
  depthsum_persamp = Counter()
  hidepthcount_persamp = Counter()
  triedsum_persamp = Counter()
  refcount_persamp = Counter()

  refcount_perdepth = Counter()
  homcount_perdepth = Counter()
  hetcount_perdepth = Counter()
  
  for line in inf :
    if re.match(r'\#CHROM', line) :
      samples = line.rstrip().split('\t')[9:]
      continue
    elif re.match(r'\#', line) :
      continue
      
    # data line
    pieces = line.rstrip().split('\t')
    chrom = pieces[0]
    # Handle falciparum chrom names
    match = re.search(r'^Pf3D7_(\d+)_v3', chrom)
    if match :
      chrom = int(match.group(1))
    else:
      continue
    pos = int(pieces[1])
    ref_all = pieces[3]
    alt_alls = pieces[4].split(',')
    kill = False
    if not re.match(r'^[ACGT\.]$', ref_all) :
      nindel += 1
      if kill_indel : continue
    for alt_all in alt_alls :
      if kill_indel and not re.match(r'^[ACGT\.]$', alt_all) :
        kill = True
    if kill :
      nindel += 1
      if kill_indel : continue

    formats = pieces[8].split(':')
    genotypes = pieces[9:]
    outline = '{:d}\t{:d}'.format(chrom, pos)
    nassay = 0
    for isamp, samp in enumerate(samples) :
      genotype = genotypes[isamp].split(':')
      call = ''
      depth = -1
      for iform, form in enumerate(formats) :
        if form == 'GT' :
          call = genotype[iform]
        elif form == 'DP' :
          if len(genotype) <= iform : depth = 0
          elif genotype[iform] == '.' : depth = 0 #line 83
          else : depth = int(genotype[iform])
          depthsum_persamp[samp] += depth
          triedsum_persamp[samp] += 1
          if depth > 5 : hidepthcount_persamp[samp] += 1
      allele = '-2'
      if depth == -1 : print('uh')
      g1 = g2 = '-2'
      if re.match(r'\d+', call) and depth >= min_depth :
        match = re.match(r'(\d+)[\/|](\d+)', call)
        g1 = match.group(1)
        g2 = match.group(2)
        callcount_persamp[samp] += 1
        if g1 != g2 :
          allele = '-1'
          hetcount_persamp[samp] += 1
          hetcount_perdepth[depth] += 1
        else :
          allele = g1
          homcount_perdepth[depth] += 1
          if allele == '0' :
            refcount_persamp[samp] += 1
            refcount_perdepth[samp] += 1

  print('N indel:', nindel)
  with open(outfile, 'w') as outf :
    print("samp\tN_hetcall\tN_call\tN_snp\tmean_depth\tdepth_gt5x\tnonref_fract", file=outf);
    for samp in samples :
      if callcount_persamp[samp] == 0 or callcount_persamp[samp] == hetcount_persamp[samp] :
        nonref_fract = 0
      else :
        nonref_fract = 1 - refcount_persamp[samp] / (callcount_persamp[samp] - hetcount_persamp[samp])
      print('{:s}\t{:d}\t{:d}\t{:d}\t{:.2f}\t{:.3f}\t{:.4f}'.
             format(samp, hetcount_persamp[samp], callcount_persamp[samp], triedsum_persamp[samp],
                    depthsum_persamp[samp] / triedsum_persamp[samp],
	            hidepthcount_persamp[samp] / triedsum_persamp[samp], nonref_fract), file=outf)
    
  
main()