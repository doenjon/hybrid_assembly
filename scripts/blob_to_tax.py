#!/usr/bin/python3

'''
Python script to convert blobtools output to taxonomic classification of contigs

usage: python3 blob_to_tax.py --blob <blob.json> [--blast <blast_file> ]

output: plots of number and classification of genes found on contigs

Save Entrez username as api_key to environment variables ENTREZ_USER and ENTREZ_API

'''

import sys
import json
import time
import traceback 
import argparse
from Bio import Entrez
import os 

Entrez.email = os.environ.get('ENTREZ_USER')
Entrez.api_key = os.environ.get('ENTREZ_API')

parser = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('--blob',
					help='blobtools json file',
					type=str,
					required=True)

parser.add_argument('--blast',
					help='blast file to annotate contigs',
					type=str,
					required=False)

args = parser.parse_args()

# Read diamond hits
hits = {}	# contig - annotations

if args.blast is not None:

	# for progress meter
	line_count = 0
	with open(args.blast, "r") as d_file:
		for line in d_file:
			line_count += 1

	# Process blast hits	
	processed = 0
	with open(args.blast, "r") as d_file:
		for line in d_file:
			processed += 1
			contig, _, _, _, qstart, qend, sseqid, stitle, *_ = line.strip().split("\t")

			qstart = int(qstart)
			qend = int(qend)

			tmp = qstart
			qstart = min(qstart, qend)
			qend = max (tmp, qend)


			# Don't process overlapping hits
			process = False
			if contig not in hits.keys():
				process = True
			else:
				process = True
				for hit in hits[contig]:
					if (hit["start"] <= qstart <= hit["end"]) or (hit["start"] <= qend <= hit["end"]):
						process = False
						break
			
			if process:
				try:
					annotations = [stitle, f" qrange:({qstart} to {qend})"]

					if contig not in hits.keys():
						hits[contig] = []	
					hits[contig].append({"annotation": ";".join(annotations), "start":qstart, "end":qend, "sseqid":sseqid})

				except Exception as e:
					traceback.print_exc()
					print(e) 

		count_accessions = 0
		for contig in hits.keys():
			# Annotate just the top hit for each contig
			if not any(keyword in hits[contig][0]["annotation"].lower() for keyword in ["retro", "gag", "pol", "integrase", "copia"]):
				# print(f"get {hits[contig][0]}")
				count_accessions += 1

		print(f"annotating top {count_accessions} accessions")
		processed = 0
		for contig in hits.keys():
			# Annotate just the top hit for each contig
			if not any(keyword in hits[contig][0]["annotation"].lower() for keyword in ["retro", "gag", "pol", "integrase", "copia"]):
				processed += 1
				try:
					sseqid = hits[contig][0]["sseqid"]
					handle = Entrez.efetch(db="protein", id=sseqid, RetMax=1, retmode="xml")
					response = Entrez.read(handle)

					feature_table = response[0]['GBSeq_feature-table']
					for feat in feature_table:

							feat_annotations = []
							if feat["GBFeature_key"] == "Protein":
								for i in range(len(feat["GBFeature_quals"])):
									if feat["GBFeature_quals"][i]["GBQualifier_name"] in ["product","name"]:
										feat_annotations.append(feat["GBFeature_quals"][i]["GBQualifier_value"])
							elif feat["GBFeature_key"] == "Region":
								for i in range(len(feat["GBFeature_quals"])):
									if feat["GBFeature_quals"][i]["GBQualifier_name"] in ["region_name", "note"]:
										feat_annotations.append(feat["GBFeature_quals"][i]["GBQualifier_value"])
							if feat_annotations:
								hits[contig][0]["annotation"] += "/".join(feat_annotations)

					print(f"{( 100 * processed / count_accessions):.2f}%  {hits[contig][0]['annotation'][:120]}", end="\n", flush=True)
				
				except Exception as e:
					print(e)
					traceback.print_exc()


# PARSE BLOBFILE - write output file
with open("blob_taxonomy.tsv", "w") as out_file:

	header = ",".join(["Name", "length", "coverage", "species", "genus", "family", "order", "phylum", "superkingdom"])
	out_file.write(header + "\n")
	with open(args.blob, "r") as blob_file:
		blob = json.load(blob_file)

		blobs = blob["dict_of_blobs"]

		for key in blobs.keys():
			contig = blobs[key]

			name = contig["name"]
			length = str(contig["length"])
			coverage = str(contig["read_cov"]["bam0"])
			tax = contig["taxonomy"]["bestsum"]

			tax_list = []
			for val in ["species", "genus", "family", "order", "phylum", "superkingdom"]:
				tax_list.append(tax[val]["tax"])

			if name in hits.keys():
				annotation = [hits[name][i]["annotation"] for i in range(len(hits[name]))]
				annotation_text = f"{len(hits[name])} genes found "
				for a in hits[name]:
					annotation_text += f"\t{a['annotation']}"
				out_line = "\t".join([name, length, coverage, *tax_list, annotation_text])
			else:
				out_line = "\t".join([name, length, coverage, *tax_list])

			out_file.write(out_line + "\n")

