#!/usr/local/bin/python3

# Script to rename samples in a vcf to the standard format (underscores, not hyphens) and restrict sites to the
#  preferred set
import gzip
import re
import sys

def main() :
  
  if len(sys.argv) != 4 : sys.exit('Usage: clean_vcf.py <input vcf file name> <output vcf file name> <file of snps to use>')
  infile = sys.argv[1]
  outfile = sys.argv[2]
  goodfile = sys.argv[3]
  #goodfile = '/Users/sfs/work/malaria/pf3k_2/output/snplist_preferred_b.txt'

  good_sites = set()
  with open(goodfile, 'r') as goodf :
    for line in goodf :
      pieces = line.rstrip().split('\t')
      good_sites.add( (int(pieces[0]), int(pieces[1])) ) # chrom, pos as ints
    
  
  # Decide whether the file is compressed or not, and open it accordingly
  if re.search(r'\.gz$', infile) :
    inf = gzip.open(infile, 'rt')
  else :
    inf = open(infile, 'r')

  samples = []
  
  nline = 0
  init_format = False
  with open(outfile, 'w') as outf :
    for line in inf :
      if re.match(r'\#CHROM', line) :
        pieces = line.rstrip().split('\t')
        samples = pieces[9:]
        print(len(samples), 'samples')
        print(*pieces[:9], sep='\t', end='', file=outf)
        for samp in samples :
          samp = samp.replace('-', '_')
          print('\t', samp, sep='', end='', file=outf)
        print('', file=outf)
        continue
      elif re.match(r'\#', line) :
        print(line, end='', file=outf)
        continue
      
      # data line
      line = line.rstrip()
      nline += 1
      pieces = line.split('\t')
      chrom_str = pieces[0]
      filter_state = pieces[6]
      if 'PASS' not in filter_state :
        continue
      # Handle falciparum chrom names
      match = re.search(r'^Pf3D7[_\-](\d+)[_\-]v3', chrom_str)
      if match :
        chrom = int(match.group(1))
      else :
        continue
      pos = int(pieces[1])
      if (chrom, pos) not in good_sites : continue
      genotypes = pieces[9:]

      # Fix the 0/0 depth=0 awkwardness
      outline = '\t'.join(pieces[:9])
      for isamp, samp in enumerate(samples) :
        field = genotypes[isamp]
        genotype = field.split(':')
      
        if not init_format :
          init_format = True
          formats = pieces[8].split(':')
          for iform, form in enumerate(formats) :
            if form == 'GT' :
              gt_idx = iform
            elif form == 'DP' :
              dp_idx = iform
        
        if genotype[dp_idx] == '.' :
          depth = 0
        else :
          depth = int(genotype[dp_idx])
        if depth == 0 :
          field = re.sub('^0[\/|]0', './.', field)
        outline += ('\t' + field)
      print(outline, file=outf)
      
main()
