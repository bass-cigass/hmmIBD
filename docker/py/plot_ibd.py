#!/usr/local/bin/python3
# Look at distribution of pairs enriched for low real IBD, to figure out acceptance region. Or try, anyway.

from collections import defaultdict
from itertools import combinations 
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import seaborn as sns
import sys

def main() :
  if len(sys.argv) != 2 : sys.exit('usage: probe_ibd.py <hmm fract file>')
  infile = sys.argv[1]
  pdf_file = 'results/plot_ibd.pdf'
  cols = sns.color_palette("colorblind", n_colors=3)
  site_name_map = {}
  site_name_map['SES'] = 'DBL'
  site_name_map['MAD'] = 'TBA'

  pp = PdfPages(pdf_file)

  ibd_by_site = defaultdict(list)
  ibd_all = []
  ngen_all = []
  with open(infile, 'r') as inf :
    idx = {}
    head = inf.readline().rstrip().split()
    idx = {col : i for i, col in enumerate(head)}
    for line in inf :
      pieces = line.rstrip().split()
      samp1 = pieces[ idx['sample1'] ]
      samp2 = pieces[ idx['sample2'] ]
      subp = samp1.split('_')
      site1 = subp[1]
      if site1 in site_name_map : site1 = site_name_map[site1]
      subp = samp2.split('_')
      site2 = subp[1]
      if site2 in site_name_map : site2 = site_name_map[site2]
      
      ngen = float(pieces[ idx['N_generation'] ])
      ntrans = float(pieces[ idx['N_state_transition'] ])
      fIBD = float(pieces[ idx['fract_sites_IBD'] ])

      ibd_all.append(fIBD)
      if site1 == site2 :
        ibd_by_site[site1].append(fIBD)
      ngen_all.append(ngen)
    # end data reading
  # end open file

  fig, ax = plt.subplots()
  ax.hist(ibd_all, rwidth=0.85, bins=np.arange(0, 1.01, .01))
  ax.set_xlabel('Fraction of genome IBD')
  ax.set_ylabel('Number of sample pairs')
  ax.set_title('Distribution of fraction IBD (all sites)')  
  ax.set_yscale('log')
  fig.savefig(pp, format='pdf', bbox_inches='tight')

  sites = ibd_by_site.keys()
  for site in sites :
    fig, ax = plt.subplots()
    ax.hist(ibd_by_site[site], rwidth=0.85, bins=np.arange(0, 1.01, .01))
    ax.set_xlabel('Fraction of genome IBD')
    ax.set_ylabel('Number of sample pairs')
    ax.set_title('Distribution of fraction IBD within ' + site)  
    ax.set_yscale('log')
    fig.savefig(pp, format='pdf', bbox_inches='tight')
    
  pp.close()

main()
