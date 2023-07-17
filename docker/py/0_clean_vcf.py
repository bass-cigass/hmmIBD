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
      good_sites.add( (pieces[0], pieces[1]) ) # chrom, pos as strings
    
  
  # Decide whether the file is compressed or not, and open it accordingly
  if re.search(r'\.gz$', infile) :
    inf = gzip.open(infile, 'rt')
  else :
    inf = open(infile, 'r')

  samples = []
  
  nline = 0
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
      nline += 1
      pieces = line.rstrip().split('\t')
      chrom = pieces[0]
      filter_state = pieces[6]
      if filter_state != 'PASS' :
        continue
      # Handle falciparum chrom names
      match = re.search(r'^Pf3D7_(\d+)_v3', chrom)
      if match :
        chrom = match.group(1)
      else :
        continue
      pos = pieces[1]
      if (chrom, pos) not in good_sites : continue
      print(line, end='', file=outf)
      
main()
