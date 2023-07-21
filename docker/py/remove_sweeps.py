#!/usr/local/bin/python3

import sys

def main() :
  sweeps = {}
  sweeps[4] = []
  sweeps[4].append( (600_000,850_000) )
  sweeps[6] = []
  sweeps[6].append( (1_050_000,1_400_000) )
  sweeps[7] = []
  sweeps[7].append( (0,700_000) )
  sweeps[8] = []
  sweeps[8].append( (450_000,800_000) )

  
  if len(sys.argv) != 3 : sys.exit('usage: remove_sweeps.py <input seq> <output file>')

  with open(sys.argv[1], 'r') as seqf, open(sys.argv[2],'w') as outf :
    head = seqf.readline().rstrip()
    print(head, file=outf)
    for line in seqf :
      pieces = line.rstrip().split()
      chrom = int(pieces[0])
      skip = False
      if chrom in sweeps :
        pos = int(pieces[1])
        for sweep in sweeps[chrom] :
          if pos >= sweep[0] and pos <= sweep[1] :
            skip = True
            break
      if skip : continue
      print(line.rstrip(), file=outf)

main()
