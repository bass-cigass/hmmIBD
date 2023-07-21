#!/usr/local/bin/python3

import sys

def main() :
  if len(sys.argv) != 4 : sys.exit('usage: thin_seq.py <input thinning file> <input seq or freq file> <output file>')
  good = set()
  with open(sys.argv[1], 'r') as inf :
    for line in inf :
      pieces = line.rstrip().split()
      chrom = pieces[0]
      pos = pieces[1]
      good.add( (chrom, pos) )

  with open(sys.argv[2], 'r') as seqf, open(sys.argv[3],'w') as outf :
    head = seqf.readline().rstrip()
    print(head, file=outf)
    for line in seqf :
      pieces = line.rstrip().split()
      chrom = pieces[0]
      pos = pieces[1]
      if (chrom, pos) in good : print(line.rstrip(), file=outf)
main()
