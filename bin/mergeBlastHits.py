#!/usr/bin/env python3

'''

Written by Nick Toda

'''

import pandas as pd
import argparse
import os.path

def main_func(args):
  '''
  Main function body. 
  Takes in an args object containing the command line arguments
  passed to the program. 
  '''
  mytab = pd.read_csv(args.infile, sep='\t',header=None)
  if args.method != "hit":
    mytab = mytab.sort_values(by=[0,1,8,9])
  j=0
  while j<(len(mytab)-1):
    if mytab.iloc[j,0]==mytab.iloc[j+1,0] and mytab.iloc[j,1]==mytab.iloc[j+1,1]: 
      # Option to merge exact hits with start and stop nearly identical
      if args.method == "exact":
        if (abs(mytab.iloc[j,8]-mytab.iloc[j+1,8])<=args.distance) and (abs(mytab.iloc[j,9]-mytab.iloc[j+1,9])<=args.distance):
          if mytab.iloc[j,10] <= mytab.iloc[j+1,10]:
            mytab=mytab.drop(mytab.index[j+1])
          else:
            mytab=mytab.drop(mytab.index[j])
        elif mytab.iloc[j,8] < mytab.iloc[j+1,8] and mytab.iloc[j,9] > mytab.iloc[j+1,9]:
          mytab=mytab.drop(mytab.index[j+1])
        elif mytab.iloc[j,8] > mytab.iloc[j+1,8] and mytab.iloc[j,9] < mytab.iloc[j+1,9]:
          mytab=mytab.drop(mytab.index[j])
        else:
          j+=1
      # Option to merge hits that are next to each other, all coordinates must have stop greater than start
      # does not combine overlaps unless completely contained
      elif args.method == "flanking_positive":
        if abs(mytab.iloc[j+1,8]-mytab.iloc[j,9])<=args.distance:
          mytab.iloc[j,9] = mytab.iloc[j+1,9]
          mytab.iloc[j,7] = mytab.iloc[j+1,7]
          # Update evalue to the most significant
          if mytab.iloc[j,10] <= mytab.iloc[j+1,10]:
            mytab.iloc[j,10] = mytab.iloc[j+1,10]
          else:
            mytab.iloc[j+1,10] = mytab.iloc[j,10]
          mytab=mytab.drop(mytab.index[j+1])
        # One or the other hit is completely contained by the other
        elif mytab.iloc[j,8] < mytab.iloc[j+1,8] and mytab.iloc[j,9] > mytab.iloc[j+1,9]:
          mytab=mytab.drop(mytab.index[j+1])
        elif mytab.iloc[j,8] > mytab.iloc[j+1,8] and mytab.iloc[j,9] < mytab.iloc[j+1,9]:
          mytab=mytab.drop(mytab.index[j])
        else:
          j+=1
      else:
        j+=1
    else:
      j+=1
  mytab.to_csv(args.outfile,sep="\t",header=False,index=False)
  return


parser = argparse.ArgumentParser(description='Merge split blast hits into single')
parser.add_argument('infile', metavar='blasthits', 
                    help='Blast output file')
parser.add_argument('outfile', metavar='output', 
                    help='Output file name')
parser.add_argument('method', metavar='method', 
                    help='''Type of merge to do: overlap, exact, flanking, hit. 
                      overlap: combine hits that overlap
                      exact: combine hits with start and stop positions within distance
                      flanking: combine hits within distance
                      hit: combine blast hits that have multiple adjacent hits within distance
                      Except for hit, different queries are combined if they overlap''')
parser.add_argument('distance', metavar='distance', 
                    help='Distance to merge across', type=int)

args = parser.parse_args()

'''
Error handling. File must exist.
'''
if not args.method in ['overlap', 'exact', 'flanking', 'flanking_positive']:
  parser.print_help()
  raise Exception(args.method," not a recognized method. Must be overlap, exact, flanking, hit")
if not os.path.isfile(args.infile):
  parser.print_help()
  raise Exception(args.infile," does not exist")
  
main_func(args)



