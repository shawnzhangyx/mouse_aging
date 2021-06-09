#!/bin/python

import argparse
from multiprocessing import Pool

parser = argparse.ArgumentParser(description='filter bam based on QNAMES')
parser.add_argument('--bam-prefix', type=str, dest="bam", help='bam file')
parser.add_argument('--statH', type=str, dest="statH", help='input statH matrix')
parser.add_argument('-o', '--outPrefix', type=str, dest="outPrefix", help='output prefix')
NCPU = 30
args = parser.parse_args()

import numpy as np
import pysam
from time import perf_counter as pc

def run():
  """ Run standard NMF on rank """
  start_time = pc()
  """ init input files """
  bamf = args.bam
  statHf = args.statH
  outPrefix = args.outPrefix
  print("filter out bam files")
  generate_bams(bamf, statHf, outPrefix)
  end_time = pc()
  print('Used (secs): ', end_time - start_time)

def generate_bams(bamf, statHf, prefix):
  o_stat_H = np.genfromtxt(statHf, dtype=None, names=True)
  cluster = np.max(o_stat_H['cluster']).astype(int) + 1
  p = Pool(NCPU)
  for idx in np.unique(o_stat_H['cluster']):
    for age in ["03","10","18"]: 
      for rep in ["rep1","rep2"]:
        p.apply_async(generate_bam_worker, (bamf, o_stat_H, idx, age, rep, prefix))
  p.close()
  p.join()

def generate_bam_worker(bamf, o_stat_H, cluster,age,rep, prefix):
  name = bamf + "." + age + "." + rep + ".dedup.filter.bam"
  #print(name)
  bamF = pysam.AlignmentFile(name)
  qnames = list(o_stat_H[np.where( (o_stat_H['cluster']==cluster) & 
      ( o_stat_H['stage'] == int(age) ) & 
      ( o_stat_H['rep'] == str.encode(rep) )
      )]['barcode'].astype(str))
  qnames_set = set(qnames)
  bam_fname = prefix + "." + "metacell_" + str(cluster) + "." + age + "." + rep + ".bam"
  print("For metaCell =", cluster, "The filtered bam is writing to:", bam_fname)
  obam = pysam.AlignmentFile(bam_fname, "wb", template=bamF)
  for b in bamF.fetch(until_eof=True):
    if b.query_name.split(':')[0] in qnames_set:
      obam.write(b)
  obam.close()
  bamF.close()
  print("metaCell =", cluster,"Stage= ",stage, "Rep= ",rep," witting finished.")

 

if __name__ == "__main__":
  """filter bam based on QNAMES"""
  run()