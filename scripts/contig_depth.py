#!/usr/bin/python3

'''
Python script to calculate and plot contig read depth

usage: python3 contig_depth.py samtools_depth_file

output: plots read depth over contigs

'''
import sys
import plotly.express as px
import numpy as np
import scipy.stats as stats
import time
import math

# Get number reads for progress statistics
ct1 = 0
print("prepping")
with open(sys.argv[1], 'r') as file:
	for line in file:
		ct1 += 1

print(ct1)

print("running")
start = time.time()
ct2 = 0
ave_cov = []
max_pos = []
current_contig = None
covs = []

# Coallate coverage statistics
with open(sys.argv[1], "r") as file:

	for line in file:

		ct2 += 1
		if ct2 % 50000 == 0:
			curr = time.time()
			print(f"{(100*float(ct2) / ct1):.1f}% - estimated {math.floor(((curr - start) / (ct2/ct1) - (curr - start)))} seconds remaining", end="\r")	


		contig, pos, cov = line.strip().split()
		pos = int(pos)
		cov = int(cov)

		if contig != current_contig:
			if len(covs) > 0:
				ave_cov.append(sum(covs) / len(covs))
				max_pos.append(pos)
				covs = []
			current_contig = contig

		else:
			covs.append(cov)

# Calculate and print statistics to console
bins = np.linspace(0, 1000, 201)
statistic, bin_edges, binnumber = stats.binned_statistic(
    x=ave_cov, values=max_pos, statistic='sum', bins=bins)
print(ave_cov)
print(max_pos)
print(statistic)
print(bin_edges)
print(binnumber)


# Plot results
fig = px.bar(y=statistic, x=bin_edges[0:200])
fig.write_html("cov.html")
