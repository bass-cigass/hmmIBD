#!/usr/local/bin/python3

from collections import Counter
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import re
import seaborn as sns
import sys

def main() :

  if len(sys.argv) != 2 :
    sys.exit('usage: hetrate.py takes at least 1 arg: gives the <tag>')
  tag = sys.argv[1]
  het_thresh = .002   # current definition of polygenomic (try 0.06 for more monogenomics sample) (nhet/ncal)
  if len(sys.argv) == 3 : 
    het_thresh = float(sys.argv[1])
    print('runing heterate with threshold '+str(het_thresh))
  cover_thresh = 0.75   # current def of good sample -- threshold on fraction of genome >= 5x cover
  cols = sns.color_palette("colorblind", n_colors=3)
  #steve adds
  year_idx = 1
  samp_sep = '_'
  do_year = True
  if tag == '2019' :
    do_year = False
  
  #end steve adds
  
  goutfile = 'output/'+tag+'good_mono_samples.txt'
  goutf = open(goutfile, 'w')

  aoutfile = 'output/'+tag+'_mono_samples.txt'
  aoutf = open(aoutfile, 'w')

  bmoutfile = 'output/'+tag+'bad_mono_samples.txt'
  bmoutf = open(bmoutfile, 'w')

  apoutfile = 'output/'+tag+'_poly_samples.txt'
  apoutf = open(apoutfile, 'w')
  
  poutfile = 'output/'+tag+'good_poly_samples.txt'
  poutf = open(poutfile, 'w')

  pall_outfile = 'output/'+tag+'good_samples.txt'
  pall_outf = open(pall_outfile, 'w')
  
  hetfile = 'output/'+tag+'_samp_het.txt'
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
  hetrate_good = []
  mean_depth_good = []
  hetrate_odd = []
  mean_depth_odd = []
  fract5x_good = []
  fract5x_all = []
  depth5x_list = []
  year_list = []
  het_list = []
  samp_list = []
  for line in hetf :
    pieces = line.rstrip().split()
    samp = pieces[0]
    if tag != '2019' and 'SEN' not in samp : continue
    subp = samp.split(samp_sep)
    ntot = int(pieces[ idx['N_call'] ])
    if ntot == 0 : continue
    if do_year : 
      year = int(subp[year_idx])
      year_list.append(year)
    nhet = int(pieces[ idx['N_hetcall'] ])
    nsnp = int(pieces[ idx['N_snp']])
    f_cover = float(pieces[ idx['depth_gt5x'] ])
    het_fract = nhet / ntot 
    depth5x_list.append(f_cover)
    het_list.append(het_fract)
    samp_list.append(samp)
    fract5x_all.append(f_cover)
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
      hetrate_good.append(het_fract)
      fract5x_good.append(f_cover)
      mean_depth_good.append( float(pieces[ idx['mean_depth'] ]) )
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

  print('N SNPs:', nsnp)
  print('N samples:', nsamp)
  print('N mono samples:', nsamp_mono)
  print('N good samples:', nsamp_good)
  print('N good mono samples:', nsamp_mono_good)
  print('N good poly samples:', nsamp_poly_good)
  print('N all poly samples:', nsamp_poly_all)
  print('N bad mono samples:', nsamp_mono_bad)

  pdf_file = 'results/' + tag + '_hetrate.pdf'
  pp = PdfPages(pdf_file)

  hets = [float(x) for x in hetrate.values()]
  fig, ax = plt.subplots()
  ax.hist(x=hets, bins=np.arange(0,0.03,.0006), rwidth=0.85, color=cols[0])
  ax.set_xlabel('Fraction of heterozygous calls')
  ax.set_ylabel('Number of samples')
  ax.set_title('Distribution of het call rate per sample (all samples)')
#  ax.set_ylim(top=40)
  fig.savefig(pp, format='pdf', bbox_inches='tight')

  fig1, ax1 = plt.subplots()
  ax1.hist(x=hetrate_good, bins=np.arange(0,0.03,.0005), rwidth=0.85, color=cols[0])
  ax1.set_xlabel('Het call rate')
  ax1.set_ylabel('Number of samples')
  ax1.set_title('Distribution of het call rate per sample (good samples)')
  fig1.savefig(pp, format='pdf', bbox_inches='tight')

  fig1, ax1 = plt.subplots()
  ax1.hist(x=mean_depth_good, rwidth=0.85, color=cols[0], bins=np.arange(0, 200, 5))
  ax1.set_xlabel('Mean depth')
  ax1.set_ylabel('Number of samples')
  ax1.set_title('Mean read depth (good samples)')
  fig1.savefig(pp, format='pdf', bbox_inches='tight')

  fig1, ax1 = plt.subplots()
  ax1.hist(x=fract5x_all, rwidth=0.85, color=cols[1], bins=np.arange(0, 1, .02))
  ax1.set_xlabel('Fraction of sites with read depth > 5')
  ax1.set_ylabel('Number of samples')
  ax1.set_title('Fraction good coverage')
  fig1.savefig(pp, format='pdf', bbox_inches='tight')

  if do_year :
    fig, ax = plt.subplots()
    ax.scatter(rand_jitter(np.array(year_list), 0.15), depth5x_list, color=cols[0], marker='.', s=3)
    ax.set_xlabel('Sample year')
    ax.set_ylabel('Fraction of genome with >5x depth')
    ax.set_title('Coverage vs year')
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    fig.savefig(pp, format='pdf', bbox_inches='tight')
    
  fig, ax = plt.subplots()
  ax.scatter(depth5x_list, het_list, color=cols[0], marker='.', s=3)
  ax.set_ylabel('Het rate')
  ax.set_xlabel('Fraction of genome with >5x depth')
  ax.set_title('Het rate vs coverage')
  ax.set_ylim(top=0.02, bottom=0)
  fig.savefig(pp, format='pdf', bbox_inches='tight')

  fig, ax = plt.subplots()
  ax.scatter(depth5x_list, het_list, color=cols[0], marker='.', s=3)
  ax.set_ylabel('Het rate')
  ax.set_xlabel('Fraction of genome with >5x depth')
  ax.set_title('Het rate vs coverage (zoomed in)')
  ax.set_ylim(top=0.007, bottom = 0)
  for isamp, samp in enumerate(samp_list) :
    if het_list[isamp] < 0.0005 and depth5x_list[isamp] > 0.6 and depth5x_list[isamp] < 0.8 :
      # print(samp, het_list[isamp], depth5x_list[isamp])
      #f isamp %3 == 0 :
      ax.annotate(samp, (depth5x_list[isamp], het_list[isamp]), textcoords='offset points', xytext=(0,0), fontsize=1)
  fig.savefig(pp, format='pdf', bbox_inches='tight')

  fig, ax = plt.subplots()
  ax.scatter(depth5x_list, het_list, color=cols[0], marker='.', s=3)
  ax.set_ylabel('Het rate')
  ax.set_xlabel('Fraction of genome with >5x depth')
  ax.set_title('Het rate vs coverage (zoomed in)')
  ax.set_ylim(top=0.007, bottom = 0)
  for isamp, samp in enumerate(samp_list) :
    pass
  fig.savefig(pp, format='pdf', bbox_inches='tight')
  
  pp.close()

def rand_jitter(l):
  stdev = .01*(max(l)-min(l))
#  stdev = .15
  return l + np.random.randn(len(l)) * stdev



main()
