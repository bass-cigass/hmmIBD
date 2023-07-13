#!/usr/local/bin/python3

from collections import Counter
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import re
import seaborn as sns
import sys

def main() :

  het_thresh = .002   # current definition of polygenomic
  cover_thresh = 0.25   # current def of good sample -- threshold on fraction of genome >= 5x cover
  cols = sns.color_palette("colorblind", n_colors=3)
  
  goutfile = 'output/good_mono_samples.txt'
  goutf = open(goutfile, 'w')

  aoutfile = 'output/all_mono_samples.txt'
  aoutf = open(aoutfile, 'w')

  bmoutfile = 'output/bad_mono_samples.txt'
  bmoutf = open(bmoutfile, 'w')

  apoutfile = 'output/all_poly_samples.txt'
  apoutf = open(apoutfile, 'w')
  
  poutfile = 'output/good_poly_samples.txt'
  poutf = open(poutfile, 'w')

  pall_outfile = 'output/good_samples.txt'
  pall_outf = open(pall_outfile, 'w')
  
  hetfile = 'output/samp_het.txt'
  hetf = open(hetfile, 'r')
  head = hetf.readline().rstrip().split()
  idx = {}
  for i in range(len(head)) :
    idx[head[i]] = i

  nsamp = 0
  nsamp_mono = 0
  nsamp_good = 0
  nsamp_mono_good = 0
  nsamp_mono_bad = 0
  nsamp_poly_all = 0
  nsamp_poly_good = 0
  hetrate = {}
  hetrate_good = {}
  for line in hetf :
    pieces = line.rstrip().split()
    samp = pieces[0]
    if '-' in samp:
      year = samp.split('-')[2]
      site = samp.split('-')[1]
    else :
      year = samp.split('_')[2]
      site = samp.split('_')[1]
    
    nhet = int(pieces[ idx['N_hetcall'] ])
    ntot = int(pieces[ idx['N_call'] ])
    if ntot == 0 : continue
    f_cover = float(pieces[ idx['depth_gt5x'] ])
    het_fract = nhet / ntot 
    samp = pieces[ idx['samp'] ]
    if re.search(r'Control', samp) or re.search(r'Ctrl', samp) :
      print('skipping', samp)
      continue
    hetrate[samp] = het_fract
    nsamp += 1
    if het_fract <= het_thresh : 
      print(samp, file=aoutf)
      nsamp_mono += 1
    if f_cover >= cover_thresh :
      hetrate_good[samp] = het_fract
      nsamp_good += 1
      # coverage > threshold and either poly or mono
      print(samp, file=pall_outf)
      
      if het_fract <= het_thresh : 
        # monogenomic and coverage > threshold
        print(samp, file=goutf)
        nsamp_mono_good += 1
    else :    # bad samples
      if het_fract <= het_thresh : 
        # monogenomic and coverage < threshold
        print(samp, file=bmoutf)
        nsamp_mono_bad += 1
    if het_fract > het_thresh : 
      print(samp, file=apoutf)
      nsamp_poly_all += 1
      if f_cover >= cover_thresh :
        nsamp_poly_good += 1
        print(samp, file=poutf)
    
  print('N samples:', nsamp)
  print('N mono samples:', nsamp_mono)
  print('N good samples:', nsamp_good)
  print('N good mono samples:', nsamp_mono_good)
  print('N good poly samples:', nsamp_poly_good)
  print('N all poly samples:', nsamp_poly_all)
  print('N bad mono samples:', nsamp_mono_bad)

  pdf_file = 'results/hetrate.pdf'
  pp = PdfPages(pdf_file)
  fig, ax = plt.subplots()
  hets = [float(x) for x in hetrate.values()]
  print(len(hets))
  ax.hist(x=hets, bins=np.arange(0,0.03,.0006), rwidth=0.85, color=cols[0])
  ax.set_xlabel('Fraction of heterozygous calls')
  ax.set_ylabel('Number of samples')
  ax.set_title('Distribution of het call rate per sample (all samples)')
#  ax.set_ylim(top=40)
  fig.savefig(pp, format='pdf', bbox_inches='tight')

  fig1, ax1 = plt.subplots()
  hets_good = [float(x) for x in hetrate_good.values()]
  print(len(hets_good))
  ax1.hist(x=hets_good, bins=np.arange(0,0.03,.0005), rwidth=0.85, color=cols[0])
  ax1.set_xlabel('Het call rate')
  ax1.set_ylabel('Number of samples')
  ax1.set_title('Distribution of het call rate per sample (>50% high coverage)')
  fig1.savefig(pp, format='pdf', bbox_inches='tight')

  pp.close()

def rand_jitter(l):
  stdev = .01*(max(l)-min(l))
#  stdev = .15
  return l + np.random.randn(len(l)) * stdev



main()
