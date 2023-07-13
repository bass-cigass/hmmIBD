#!/usr/local/bin/python3

# Script to output genotype data per sample formatted for hmmIBD

from collections import Counter
import gzip
import re
import sys

def main() :
  kill_indel = False
  # Mininum genotyping call rate to accept variant
  min_call = 0.80
  # Option: minimum number of alt allele copies to keep variant
  #  (1 = kill monomorphic, 2 = kill singletons)
  min_copy = 0
  min_depth = 5      # minimum read depth to accept a genotype
  
  if len(sys.argv) != 2 and len(sys.argv) != 3: sys.exit('Usage: vcf2hmm.py <input vcf file name> [<file of samples to use>]')
  infile = sys.argv[1]
  seqfile = 'seq/seq.txt'
  allfile = 'seq/allele.txt'

  all_samps = True
  good_samps = set()
  if len(sys.argv) == 3 :
    samp_file = sys.argv[2]
    all_samps = False
    sampf = open(samp_file, 'r')
    for line in sampf :
      samp = line.rstrip()
      good_samps.add(samp)

  # Decide whether the file is compressed or not, and open it accordingly
  if re.search(r'\.gz$', infile) :
    inf = gzip.open(infile, 'rt')
  else :
    inf = open(infile, 'r')

  samples = []
  depthsum = Counter()    # by sample
  ndepth = Counter()
  nhet = Counter()
  ncall = Counter()
  nocall = Counter()
  nhomo = 0
  nindel = 0
  n_fail_filter = n_pass_filter = 0
  
  nline = 0
  nsamp_used = 0
  with open(seqfile, 'w') as seqf, open(allfile, 'w') as allf :
    for line in inf :
      if re.match(r'\#CHROM', line) :
        samples = line.rstrip().split('\t')[9:]
        print(len(samples), 'samples')
        print('chrom\tpos', sep='', end='', file=seqf)
        for samp in samples :
          if all_samps or samp in good_samps :
            nsamp_used += 1
            print('\t', samp, sep='', end='', file=seqf)
        print('', file=seqf)
        print(nsamp_used, 'samples used')
        continue
      elif re.match(r'\#', line) :
        continue
      
      # data line
      nline += 1
      pieces = line.rstrip().split('\t')
      chrom = pieces[0]
      filter_state = pieces[6]
      if filter_state != 'PASS' :
        n_fail_filter += 1
        continue
      n_pass_filter += 1
      # Handle falciparum chrom names
      match = re.search(r'^Pf3D7_(\d+)_v3', chrom)
      if match :
        chrom = int(match.group(1))
      else :
        # Mito and apicoplast or something is wrong
        continue
      if chrom > 14 : print('high chrom number', chrom)
      pos = int(pieces[1])
      ref_all = pieces[3]
      alt_alls = pieces[4].split(',')
      kill = False
      if not re.match(r'^[ACGT\.]$', ref_all) :
        nindel += 1
        if kill_indel : continue
      for alt_all in alt_alls :
        if not re.match(r'[ACGT\.]$', alt_all) :
          kill = True
      if kill :
        nindel += 1
        if kill_indel : continue

      formats = pieces[8].split(':')
      genotypes = pieces[9:]
      outline = '{:d}\t{:d}'.format(chrom, pos)
      nassay = 0
      
      for isamp, samp in enumerate(samples) :
        if not all_samps and samp not in good_samps : continue
        genotype = genotypes[isamp].split(':')
        call = ''
        for iform, form in enumerate(formats) :
          if form == 'GT' :
            call = genotype[iform]
          elif form == 'DP' :
            if genotype[iform] == '.' : depth = 0
            else : depth = int(genotype[iform])
            depthsum[samp] += depth
            ndepth[samp] += 1
        allele = '-2'
        if re.match(r'\d+', call) and depth >= min_depth :
          match = re.match(r'(\d+)[\/|](\d+)', call)
          g1 = match.group(1)
          g2 = match.group(2)
          ncall[samp] += 1
          if g1 != g2 :
            allele = '-1'
            nhet[samp] += 1
          else :
            allele = g1
        else :
          #nocall
          allele = '-1'
        nassay += 1
        outline += ('\t' + allele)
      print(outline, file=seqf)
      print('{:d}\t{:d}'.format(chrom, pos), ref_all, '\t'.join(alt_alls), sep='\t', file=allf)
  print('N indel', nindel)
  print('N failed/passed filter', n_fail_filter, n_pass_filter)
main()
