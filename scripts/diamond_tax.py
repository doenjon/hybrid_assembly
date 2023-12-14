#!/usr/bin/python3

'''
Python script summarize the information in a diamond blast file

usage: python3 diamond_tax.py <diamond_blast_output>

output: pdf summary

Save Entrez username as api_key to environment variables ENTREZ_USER and ENTREZ_API

'''

import sys
import time
import traceback 
from Bio import Entrez
from Bio import SeqIO
import matplotlib.pyplot as plt 

Entrez.email = os.environ.get('ENTREZ_USER')
Entrez.api_key = os.environ.get('ENTREZ_API')

# Get diamond hits file from command line
input_file = sys.argv[1]

# Go through hits file, get species counts
taxa = {}
with open(input_file) as file:
	for line in file:
		_, tax_id, *_ = line.strip().split("\t")

		if tax_id == "":
			tax_id = "unclassified"

		if tax_id in taxa.keys():
			taxa[tax_id]["count"] += 1
		else:
			taxa[tax_id] = {}
			taxa[tax_id]["count"] = 1

# Get NCBI information about the speceis in the hit file
tax_ids = list(taxa.keys())
handle = Entrez.efetch(db="taxonomy", id=tax_ids)
response = Entrez.read(handle)

# Annotate species with more taxonomic information
for r in response:
	try:
		tax_id = r["TaxId"]
		lineage = r["LineageEx"]
		sci_name = r["ScientificName"]
		ranks = {}
		ranks_of_interest = ["kingdom", "phylum", "class", "order", "family", "genus"]
		for l in lineage:
			if l["Rank"].lower() in ranks_of_interest:
				ranks[l["Rank"].lower()] = l["ScientificName"]
		taxa[tax_id]["rank"] = ranks 
		taxa[tax_id]["ScientificName"] = sci_name
	except KeyError as e:
		pass
	except Exception as e:
		taxa[tax_id]["rank"] = "unknown" 
		taxa[tax_id]["ScientificName"] = "unknown"

rem_keys = []
for key in taxa.keys():
	try:
		taxa[key]["rank"]
		taxa[key]["ScientificName"]
	except:
		taxa[key]["rank"] = "unknown"
		taxa[key]["ScientificName"] = "unknown"


###############################################################################################
# Print a pdf of top 10 species in sample

species_count = [(taxa[tax_id]["ScientificName"], taxa[tax_id]["count"]) for tax_id in taxa.keys()]
species_count.sort(key = lambda x: x[1], reverse = True)
for l in species_count:
	print(l)
print()

x = [x[0] for x in species_count][:10]
y = [y[1] for y in species_count][:10]

fig = plt.figure(figsize = (10, 8))
plt.bar(x, y, color ='maroon', width = 0.4)
plt.xticks(rotation='vertical')
plt.tight_layout()
 
plt.xlabel("Species")
plt.ylabel("Sequene Count")
plt.savefig("Species_count.pdf")

###############################################################################################
# Print a pdf of top 10 phylum in sample

cluster_rank = "phylum"
cluster_count = {}
for tax in taxa.keys():
	try:
		rank_classification = taxa[tax]["rank"][cluster_rank]
		if rank_classification in cluster_count.keys():
			cluster_count[rank_classification] += taxa[tax]["count"]
		else:
			cluster_count[rank_classification] = taxa[tax]["count"]
	except:
		print(f"lost {tax}: {taxa[tax]}")

cluster_count = [(rank, cluster_count[rank]) for rank in cluster_count.keys()]
cluster_count.sort(key = lambda x: x[1], reverse = True)
for l in cluster_count:
	print(l)
print()
x = [x[0] for x in cluster_count][:10]
y = [y[1] for y in cluster_count][:10]

fig = plt.figure(figsize = (10, 8))
plt.bar(x, y, color ='maroon', width = 0.4)
plt.xticks(rotation='vertical')
plt.tight_layout()
 
plt.xlabel(cluster_rank)
plt.ylabel("Sequene Count")
plt.savefig("tax_count.pdf")

