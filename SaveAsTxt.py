#!/usr/bin/python3
import pandas as pd
import os
import argparse

ap = argparse.ArgumentParser()

ap.add_argument('--input', help='path of input file')
ap.add_argument('--output', help='path of output file')
args = ap.parse_args()
df = pd.read_parquet(args.input)
df.to_csv(args.output, sep='\t')