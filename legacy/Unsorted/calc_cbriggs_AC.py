#!/usr/bin/env python3

import sys
import pandas as pd

arg_in = sys.argv[1]

def determine_header(file_name):
	headcount = 0
	header_line = []
	format_in = 0
	with open(file_name+".csv", "w") as new_f:
		with open(file_name, "r") as f:
			for line in f:
				if line.startswith("##"):
					headcount += 1
					new_f.write(line)
				if line.startswith("#CHROM"):
					header_line = line.split()
					format_in = header_line.index("FORMAT")
					header_len = len(header_line)
	return(headcount,format_in,header_len)

skiplines,format_col,last_col = determine_header(arg_in)
print(format_col,last_col)
df = pd.read_csv(arg_in,skiprows=skiplines,header=0,sep="\t", usecols=range(format_col+1,last_col))
new_ACinfo = df.apply(pd.to_numeric, errors='coerce').sum(1)
df.drop(df.index, inplace=True)
df = pd.read_csv(arg_in,skiprows=skiplines,header=0,sep="\t")
df["INFO"] = 'AC=' + new_ACinfo.astype(int).astype(str)
##contig=<ID=I>
new_f = open(arg_in+".csv","a")
for i in df["#CHROM"].unique():
	new_f.write("##contig=<ID="+i+">\n")
new_f.write("##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Total number of alternate alleles in called genotypes\">\n")
df.to_csv(new_f,sep="\t",index=False)
new_f.close()
