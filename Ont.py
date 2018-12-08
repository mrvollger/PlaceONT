#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser(description="")
parser.add_argument("ref", nargs="?", help="input ref file")
parser.add_argument("reads",nargs="?", help="reads file")
parser.add_argument("--fofn", help="fofn file", default=None)
parser.add_argument("-t","--threads", help="threads", type=int, default=1)
parser.add_argument('-d', action="store_true", default=False)
args = parser.parse_args()

import os
import re
import mappy as mp

ref = mp.Aligner(args.ref, preset="map-ont", bw=50000, best_n=1, n_threads=args.threads)
if not ref: raise Exception("ERROR: failed to load/build index")

def makeAlns(reads):
	hits = []
	for name, seq, qual in mp.fastx_read(reads):
		for hit in ref.map(seq):
			print(hit)

hits = []
if(args.fofn is not None):
	for f in open(args.fofn).readlines():
		f=f.srtip()
		makeAlns(f)

if(args.reads is not None):
		makeAlns(args.reads)


