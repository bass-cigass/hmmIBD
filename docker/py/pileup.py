#!/usr/local/bin/python3

from collections import Counter
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys

def main() :
  
  if len(sys.argv) != 4 : sys.exit('usage: pileup.py <hmm fraction file> <hmm file> <pdf output file>')
  infile = sys.argv[1]
  fractfile = sys.argv[2]
  pdf_file = sys.argv[3]

  pp = PdfPages(pdf_file)
  cols = sns.color_palette("colorblind", n_colors=6)
  
  binsize = 1000
  nchrom = 14

  rels = set()
  with open(fractfile, 'r') as fractf :
    head = fractf.readline().rstrip().split()
    idx = {col : i for i, col in enumerate(head)}
    for line in fractf :
      pieces = line.rstrip().split()
      samp1 = pieces[idx['sample1']]
      samp2 = pieces[idx['sample2']]
      fIBD = float(pieces[idx['fract_sites_IBD']])
      if fIBD > 0.98 : continue
      if fIBD > 0.05 : 
        rels.add( (samp1, samp2) )
  
  ntot_nonrel = ntot_rel = 0
  chr_ends = Counter()
  with open(infile, 'r') as inf :
    head = inf.readline().rstrip().split()
    idx = {col: i for i, col in enumerate(head)}
    for line in inf :
      pieces = line.rstrip().split()
      chr = int(pieces[idx['chr']])
      end = int(pieces[idx['end']])
      if end > chr_ends[chr] : chr_ends[chr] = end
    inf.seek(0)

    depth = [ [0]*(chr_ends[x]//binsize+1) for x in range(nchrom+1) ]
    depth_rel = [ [0]*(chr_ends[x]//binsize+1) for x in range(nchrom+1) ]

    inf.readline()
    olds1 = olds2 = ''
    for line in inf :
      pieces = line.rstrip().split()
      diff = pieces[idx['different']]
      s1 = pieces[idx['sample1']]
      s2 = pieces[idx['sample2']]
      if s1 != olds1 or s2 != olds2 :
        if (s1, s2) in rels : ntot_rel += 1
        else : ntot_nonrel += 1

      olds1 = s1
      olds2 = s2
      if diff == '1' : continue
      chr = int(pieces[idx['chr']])
      end = int(pieces[idx['end']])
      start = int(pieces[idx['start']])
      bin_end = end // binsize
      bin_start = start // binsize
      if (s1, s2) in rels :
        for i in range(bin_start, bin_end+1) :
          depth_rel[chr][i] += 1
      else :
        for i in range(bin_start, bin_end+1) :
          depth[chr][i] += 1

  print('relative, non-relative pairs', ntot_nonrel, ntot_rel)

  for chrom in range(1, 15) :
    xs = [(x+.5)*binsize/1000000 for x in range(len(depth_rel[chrom])) ]

    fig, ax = plt.subplots()
    if ntot_rel > 0 :
      frac_depth_rel = np.array(depth_rel[chrom] , dtype=int) / ntot_rel
      ax.scatter(xs, frac_depth_rel, s=3, color=cols[1], label='Relatives')
    if ntot_nonrel > 0 :
      frac_depth = np.array(depth[chrom] , dtype=int) / ntot_nonrel
      ax.scatter(xs, frac_depth, s=3, color=cols[0], label='Non-relatives')
    ax.set_xlabel('Position on chromosome ' + str(chrom) + ' (Mb)')
    ax.set_ylabel('Fraction of times shared')
    ax.set_title('Sharing pileup, chromosome ' + str(chrom))
    ax.set_ylim(top=1.0)
    ax.legend()
    fig.savefig(pp, format='pdf')

  pp.close()

main()
