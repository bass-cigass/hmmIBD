#!/usr/local/bin/python3

# Script to restrict sites to the preferred set and fix missing data bug (thanks, GATK!)
import gzip
import re
import sys

def main() :
  
  if len(sys.argv) != 4 : sys.exit('Usage: clean_vcf.py <input vcf file name> <output vcf file name> <file of snps to use>')
  infile = sys.argv[1]
  outfile = sys.argv[2]
  goodfile = sys.argv[3]
  #goodfile = '/Users/sfs/work/malaria/pf3k_2/output/snplist_preferred_b.txt'

  filter_fields = set()
  n_bad_sites = 0
  
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
  npass = ngood = nfiltered = n_missing_depth = 0
  init_format = False
  with open(outfile, 'w') as outf :
    for line in inf :
      if re.match(r'\#CHROM', line) :
        pieces = line.rstrip().split('\t')
        samples = pieces[9:]
        print(len(samples), 'samples')
        print(*pieces[:9], sep='\t', end='', file=outf)
        for samp in samples :
          print('\t', samp, sep='', end='', file=outf)
        print('', file=outf)
        my_count = {x : samples.count(x) for x in samples}
        for samp in my_count.keys() :
          if my_count[samp] > 1 : print('duplicate sample', samp, my_count[samp])
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
      filter_fields.add(filter_state)
      pos = int(pieces[1])
      if filter_state == 'PASS' :
        npass += 1
    # Handle falciparum chrom names
      match = re.search(r'^Pf3D7_(\d+)_v3', chrom_str)
      if match :
        chrom = int(match.group(1))
      else :
        continue
      if nline%40000 == 0 :print(nline, 'nbad', n_bad_sites, 'nfiltered', nfiltered, chrom, pos, filter_fields)
      if (chrom, pos) not in good_sites :
        n_bad_sites += 1
        continue
      #if filter_state != 'PASS' and filter_state != 'Core' :
      if filter_state != 'PASS' :
        nfiltered += 1
        continue
      ngood += 1
      genotypes = pieces[9:]
  
      # Fix the 0/0 depth=0 awkwardness
      outline = '\t'.join(pieces[:9])
      for isamp, samp in enumerate(samples) :
        field = genotypes[isamp]
        genotype = field.split(':')
      
        if not init_format or True:
          init_format = True
          formats = pieces[8].split(':')
          for iform, form in enumerate(formats) :
            if form == 'GT' :
              gt_idx = iform
            elif form == 'DP' :
              dp_idx = iform

        if len(genotype) <= dp_idx :
          #print(chrom, pos, field, genotype, isamp, samp)
          n_missing_depth += 1
          depth = 0
        elif genotype[dp_idx] == '.' :
          depth = 0
        else :
          depth = int(genotype[dp_idx])
        if depth == 0 :
          field = re.sub('^0[\/|]0', './.', field)
        outline += ('\t' + field)
      print(outline, file=outf)
  print('passing, good', npass, ngood)
  print('preferred sites that fail QC', nfiltered)
  print('filter fields encountered:', filter_fields)
  print('n sites not in preferred list:', n_bad_sites)
  print('genotypes with no depth field:', n_missing_depth)
  
main()