#!/usr/local/bin/python3

from collections import Counter
import gzip
import re
import sys

def main() :
  kill_indel = True
  # Mininum genotyping call rate to accept variant
  min_call = 0.80
  # Option: minimum number of alt allele copies to keep variant
  #  (1 = kill monomorphic, 2 = kill singletons)
  min_depth = 5

  if len(sys.argv) < 3 : sys.exit('Usage: vcf2het.py <input vcf file name> <year or all>: '+len(sys.argv)+ ' parameters are given')
  infile = sys.argv[1]
  target_year = sys.argv[2]
  outfile = 'output/' + target_year + '_samp_het.txt' #'output/'+tag+'_samp_het.txt'
  depthfile = 'output/' + target_year + '_hetdepth.txt'


  
  if re.search(r'\.gz$', infile) :
    inf = gzip.open(infile, 'rt')
  else :
    inf = open(infile, 'r')

  good_samps = set()
  samples = []
  nindel = n_lowcall = n_lowfreq = 0
  
  callcount_persamp = Counter()
  hetcount_persamp = Counter()
  depthsum_persamp = Counter()
  hidepthcount_persamp = Counter()
  triedsum_persamp = Counter()
  refcount_persamp = Counter()

  refcount_perdepth = Counter()
  homcount_perdepth = Counter()
  hetcount_perdepth = Counter()
  
  if True :
    for iline, line in enumerate(inf) :
      if re.match(r'\#CHROM', line) :
        samples = line.rstrip().split('\t')[9:]
        for samp in samples :
          subp = samp.split('_')
          if len(subp) == 1 :
            subp = samp.split('-')
          if len(subp) < 3 or subp[0] != 'SEN' :
            continue      # non-standard Senegal name
          if target_year != 'all' :
            year = subp[1]
            print(year)
            if year != target_year : continue
            if samp in good_samps : print('Duplicate sample', samp)
          good_samps.add(samp)
        print(target_year, len(good_samps), 'good samples')
        continue
      elif re.match(r'\#', line) : continue
      
      # data line
      pieces = line.rstrip().split('\t')
      chrom = str(pieces[0])
      # Handle falciparum chrom names
      match = re.search(r'^Pf3D7_(\d+)_v3', chrom)
      if match :
        chrom = int(match.group(1))
      pos = int(pieces[1])
      ref_all = pieces[3]
      alt_alls = pieces[4].split(',')
      kill = False
      if not re.match(r'^[ACGT\.]$', ref_all) :
        nindel += 1
        #print('ref', ref_all)
        if kill_indel : continue
      for alt_all in alt_alls :
        if kill_indel and not re.match(r'^[ACGT\.]$', alt_all) :
          #print('alt', alt_all)
          kill = True
      if kill :
        nindel += 1
        if kill_indel : continue

      formats = pieces[8].split(':')
      genotypes = pieces[9:]
      nassay = ncall = 0 # nassay = number of samples that could have a call, ncall= number of sample that has been called
      for isamp, samp in enumerate(samples) :
        if samp not in good_samps : continue
        if isamp >= len(genotypes) : print(iline, chrom, pos, isamp, len(genotypes))
        genotype = genotypes[isamp].split(':')
        call = ''
        depth = -1
        for iform, form in enumerate(formats) :
          if form == 'GT' :
            call = genotype[iform]
          elif form == 'DP' :
            if len(genotype) <= iform : depth = 0
            elif genotype[iform] == '.' : depth = 0
            else : depth = int(genotype[iform])
            depthsum_persamp[samp] += depth
            triedsum_persamp[samp] += 1
            if depth > 5 : hidepthcount_persamp[samp] += 1
        allele = '-2'
        if depth == -1 : print('uh')
        g1 = g2 = '-2'
        if re.match(r'[0-9]+', call) and depth >= min_depth :
          match = re.match(r'([0-9]+)[\/|]([0-9]+)', call)
          g1 = match.group(1)
          g2 = match.group(2)
          callcount_persamp[samp] += 1
          ncall += 1
          if g1 != g2 :
            # het
            allele = '-1'
            hetcount_persamp[samp] += 1
            hetcount_perdepth[depth] += 1
          else :
            #homo
            allele = g1
            homcount_perdepth[depth] += 1
            if allele == '0' :
              refcount_persamp[samp] += 1
              refcount_perdepth[depth] += 1
        nassay += 1
      if nassay == 0 or ncall / nassay < min_call :
        n_lowcall += 1
        continue
              
  print('N indel killed:', nindel)
  print('N dropped for low call rate', n_lowcall)
  with open(outfile, 'w') as outf :
    print("samp\tN_hetcall\tN_call\tN_snp\tmean_depth\tdepth_gt5x\tnonref_fract", file=outf);
    for samp in samples :
      if samp not in good_samps : continue
      if callcount_persamp[samp] == 0 or callcount_persamp[samp] == hetcount_persamp[samp] :
        nonref_fract = 0
      else :
        nonref_fract = 1 - refcount_persamp[samp] / (callcount_persamp[samp] - hetcount_persamp[samp])
      print('{:s}\t{:d}\t{:d}\t{:d}\t{:.2f}\t{:.3f}\t{:.4f}'.
             format(samp, hetcount_persamp[samp], callcount_persamp[samp], triedsum_persamp[samp],
                    depthsum_persamp[samp] / triedsum_persamp[samp],
	            hidepthcount_persamp[samp] / triedsum_persamp[samp], nonref_fract), file=outf)
  with open(depthfile, 'w') as depf :
    print("depth\tN_hom\tN_nonref\tN_hom\tnonref_fract\thet_fract\n", file=depf)
    for dep in range(150) :
      print(dep, homcount_perdepth[dep], refcount_perdepth[dep], hetcount_perdepth[dep])
      nonref = homcount_perdepth[dep] - refcount_perdepth[dep]
      if hetcount_perdepth[dep] + homcount_perdepth[dep] == 0: het_fract = 0
      else: het_fract = hetcount_perdepth[dep] / (hetcount_perdepth[dep] + homcount_perdepth[dep])
      if homcount_perdepth[dep] == 0: nonref_fract = 0
      else: nonref_fract = nonref / homcount_perdepth[dep]
      print('{:d}\t{:d}\t{:d}\t{:d}\t{:.4f}\t{:.4f}'.
            format(dep, homcount_perdepth[dep], nonref, hetcount_perdepth[dep], nonref_fract, het_fract), file=depf)
  
main()
