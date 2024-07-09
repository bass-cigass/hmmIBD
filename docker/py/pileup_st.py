#!/usr/local/bin/python3

from collections import Counter
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys

def main() :
  
  if len(sys.argv) != 6 : sys.exit('usage: pileup.py <file tag>')
  min_year = 2019
  print(f'Only analyzing year {min_year} and later')
  tag = sys.argv[1]
  infile = sys.argv[2]
  fractfile = sys.argv[3]
  pdf_file = sys.argv[4]
  gene_file = sys.argv[5]
  #infile = 'output/' + tag + '.hmm.txt'
  #fractfile = 'output/' + tag + '.hmm_fract.txt'
  #gene_file = '../sensurv/data_fixed/gene_locs.txt'
  #pdf_file = 'results/' + tag + '_pileup.pdf'
  outfile = f'output/{tag}.pileup.txt'
  pp = PdfPages(pdf_file)
  cols = sns.color_palette("colorblind", n_colors=7)
  binsize = 1000
  nchrom = 14
  year_idx = 1
  samp_sep = '_'
  outf = open(outfile, 'w')
  print('chrom\tpos\tfraction',file=outf)

  gene_start = {}
  gene_end = {}
  gene_name = {}
  with open(gene_file, 'r') as genef :
    for line in genef :
      pieces = line.rstrip().split()
      chrom = int(pieces[0])
      start = int(pieces[1]) / 1000000
      end = int(pieces[2]) / 1000000
      name = pieces[3]
      if chrom not in gene_start :
        gene_start[chrom] = []
        gene_end[chrom] = []
        gene_name[chrom] = []
      gene_start[chrom].append(start)
      gene_end[chrom].append(end)
      gene_name[chrom].append(name)

  pairs = set()
  clones = set()
  
  with open(fractfile, 'r') as fractf :
    head = fractf.readline().rstrip().split()
    idx = {col : i for i, col in enumerate(head)}
    for line in fractf :
      pieces = line.rstrip().split()
      samp1 = pieces[idx['sample1']]
      samp2 = pieces[idx['sample2']]
      if 'SEN' not in samp1 or 'SEN' not in samp2 : continue
      subp = samp1.split(samp_sep)
      year1 = int(subp[year_idx])
      subp = samp2.split(samp_sep)
      year2 = int(subp[year_idx])
      if year1 < min_year or year2 < min_year : continue
      fIBD = float(pieces[idx['fract_sites_IBD']])
      if fIBD > 0.98 :
        clones.add( (samp1, samp2) )
        continue
      pairs.add( (samp1, samp2) )
  
  ntot = 0
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
    print(chr_ends)

    depth = [ [0]*(chr_ends[x]//binsize+1) for x in range(nchrom+1) ]

    inf.readline()
    olds1 = olds2 = ''
    for line in inf :
      pieces = line.rstrip().split()
      diff = pieces[idx['different']]
      s1 = pieces[idx['sample1']]
      s2 = pieces[idx['sample2']]
      if s1 != olds1 or s2 != olds2 :
        if (s1, s2) in pairs : ntot += 1
      olds1 = s1
      olds2 = s2
      if (s1, s2) not in pairs : continue
      if diff == '1' : continue
      chr = int(pieces[idx['chr']])
      end = int(pieces[idx['end']])
      start = int(pieces[idx['start']])
      bin_end = end // binsize
      bin_start = start // binsize
      for i in range(bin_start, bin_end+1) :
        depth[chr][i] += 1
  print('pairs', ntot)


  for chrom in range(1, 15) :
    xs = [(x+0)*binsize/1000000 for x in range(len(depth[chrom])) ]
    xc = [(x+.5)*binsize/1000000 for x in range(len(depth[chrom])) ]
    xe = [(x+1)*binsize/1000000 for x in range(len(depth[chrom])) ]

    frac_depth = np.array(depth[chrom] , dtype=int) / ntot
    for ibin, pos in enumerate(xc) :
      print('{:d}\t{:.1f}\t{:.4f}'. format(chrom, pos * 1000, frac_depth[ibin]), file=outf)

    fig, ax = plt.subplots()
    init1 = False
    for i in range(len(xc)) :
      ax.plot((xs[i],xe[i]), (frac_depth[i], frac_depth[i]), color=cols[0], label='pileup', marker='', lw=1.5)
      if init1 :
        ax.plot((xe[i-1],xs[i]), (frac_depth[i-1], frac_depth[i]), color=cols[0], marker='', lw=.2)
      if frac_depth[i] > 0 : init1 = True
    ax.set_xlabel('Position on chromosome ' + str(chrom) + ' (Mb)')
    ax.set_ylabel('Fraction of times shared')
    ax.set_title('Sharing pileup, chromosome ' + str(chrom))
    #ax.set_ylim(top=0.2)
    lines = []
    labels = []
    lines.append(mlines.Line2D([], [], color=cols[0], marker=''))
    labels.append('Sharing')
    #  Add a line for candidate genes
    if chrom in gene_start :
      for igene in range(len(gene_start[chrom])) :
        #x = (gene_start[chrom][igene]/1000000, gene_end[chrom][igene]/1000000)
        #y = (0,0)
        #ax.plot(x, y, color=cols[2], label=gene_name[chrom][igene])
        lines.append(mlines.Line2D([], [], color=cols[1+igene], marker=''))
        labels.append(gene_name[chrom][igene])
        xpos = gene_start[chrom][igene] + 0.5 * (gene_end[chrom][igene] - gene_start[chrom][igene])
        color = cols[igene+1]
        plt.axvline(xpos, color=color, ls='dotted')
    ax.legend(handles=lines, labels=labels)
    fig.savefig(pp, format='pdf')
    plt.close('all')
    
  pp.close()

main()
