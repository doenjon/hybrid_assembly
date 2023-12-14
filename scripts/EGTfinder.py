###############################################################################
#
#   EGTfinder - A Targeted and simple approach to finding HGT events between 
#               an endosymbiont and host organisms
#
# 	Example usage:
# 		python3 EGTfinder.py build_db --endo CP000815 --endo_adj NZ_CH724159.1,NC_005042.1,NC_019675.1 \
#			--host_adj resources/GCA_000512085.1_Reti_assembly1.0_cds_from_genomic.fasta,\
#				resources/GCA_000698435.2_ASM69843v2_cds_from_genomic.fasta 
# 		python3 EGTfinder.py blast --host resources/Paulinella_transcriptome.fasta --evalue 1e-5
###############################################################################

import sys
import time
import re
import pickle
import os.path
import argparse
import traceback
import subprocess 
from datetime import datetime
import numpy as np
from Bio import SeqIO
from Bio import SearchIO
from Bio import Entrez
from Bio.Blast.Applications import NcbiblastxCommandline
import matplotlib.pyplot as plt 

Entrez.email = "jdoenier@stanford.edu"	# Need to fill in before running...
#Entrez.api_key = "" # bad...

class EGTfinder():
	""" 
	A class to organize implementation of blast database building search
	"""

	def __init__(self, args):
		"""
		simple init

		Args
			args - argparse object from command line options
		"""

		self.args = args

		# Load parameters needed for using the blast database if running in blast mode
		if self.args.command == "blast":
			self._load_db()

		self.threshold = 0	# Default summed bitscore require to call a hit a HGT event

	def build_db(self):
		"""
		make a database from a set of fasta files
		"""

		# Collect fasta files/accession numbers from command line input
		fastas = []
		fastas.append(self.args.endo)
		fastas.extend(self.args.endo_adj.split(","))
		if len(self.args.host_adj) > 0:
			fastas.extend(self.args.host_adj.split(","))
		
		# Download resources if given an NCBI accession
		resource_locs = []
		for fasta in fastas:
			if ".fasta" not in fasta: # file is an accession
				fasta = self._download_NCBI_fasta(fasta, self.args.force_download)
			else:
				fasta = self.args.download_location + fasta
			resource_locs.append(fasta)

		# Concatinate fasta files
		build_fastas = self.args.download_location + "all_fasta.fasta"
		with open(build_fastas, 'w') as outfile:
			for file in resource_locs:
				with open(file) as infile:
					outfile.write(infile.read())

		# Build database - inspired from https://www.biostars.org/p/249968/
		print("Building blast database...", end="\t", flush=True)
		blastdb_cmd = f'makeblastdb -in {build_fastas} -dbtype prot -title custom_EGT_db -out {self.args.db}'
		DB_process = subprocess.Popen(blastdb_cmd,
							  shell=True,
							  stdin=subprocess.PIPE,
							  stdout=subprocess.PIPE,
							  stderr=subprocess.PIPE)
		DB_process.wait()

		blastdb_cmd = f'diamond makedb --in {build_fastas} -d {self.args.db}'
		DB_process = subprocess.Popen(blastdb_cmd,
							  shell=True,
							  stdin=subprocess.PIPE,
							  stdout=subprocess.PIPE,
							  stderr=subprocess.PIPE)
		DB_process.wait()

		print("Done")	# With making database

		# Save annotations for subsequent blast searches
		pickle_file = self.args.db + ".pickle"
		endo_adj = self.args.endo_adj.split(",")
		host_adj = self.args.host_adj.split(",")
		with open(pickle_file, "wb") as file:
			pickle.dump({"endo": self.args.endo, "endo_adj": endo_adj, "host_adj": host_adj}, file)

	def _load_db(self):
		"""
		private method to load parameters required to correctly use the database
		"""
		pickle_file = self.args.db + ".pickle"
		annotations = pickle.load(open(pickle_file, "rb"))

		self.endo = annotations["endo"]
		self.endo_adj = annotations["endo_adj"]

		if len(annotations["host_adj"]) <= 1 and annotations["host_adj"][0] == "":
			self.host_adj = [] # It's an empty space...
		else:
			self.host_adj = annotations["host_adj"]

		print(f"db_params: ")
		print(f"\tendosymbiont: {self.endo}")
		print(f"\tendosymbiont related organisms: {', '.join(self.endo_adj)}")
		print(f"\thost related organisms: {', '.join(self.host_adj)}")

	def blast(self, fasta, skip_blast=False):
		"""
		blast all sequences in a fasta file against the custom blast database

		Args:
			fasta - a string NCBI accession 

		Returns:
			contigs - a list of contig names that were determined as be HGT events.
		"""
		# If given an accession number, download from NCBI
		if "fasta" not in fasta: 
			fasta = self._download_NCBI_fasta(fasta)

		print("Performing blast... ", end="\t", flush=True)
		blast_output = f"tmp_blast_{datetime.now().strftime('%s')}.xml"	# File to store blast results in

		# In some cases, it is useful to performing the blast and instead parse blast
		# results from a previous search.
		if not skip_blast:
			blast_cli = NcbiblastxCommandline(query=fasta, 
												db=self.args.db, 
												evalue=self.args.evalue, 
												outfmt=5,
												out=blast_output)
		else:
			print("skipping blast...", end="\t", flush=True)


		blastdb_cmd = f'diamond blastx -q {fasta} -d {self.args.db} -o {blast_output} -f 5 -F 15 -b4 -c2 --sensitive'
		DB_process = subprocess.Popen(blastdb_cmd,
							  shell=True)
		DB_process.communicate()
		DB_process.wait()

		print("Done") # performing blast

		total_hits = 0
		contigs = []

		# File to write results to. Could also be moved to stdout in the future
		with open(f"tmp_output.txt", "w") as out_file:	

			# Open blast result file to analyse hits.
			# for query in SearchIO.parse("tmp_blast_1622762558.xml", "blast-xml"):
			for query in SearchIO.parse(blast_output, "blast-xml"):
				print(query)

				if len(query.hits) > 0:	# The query has blast hits

					for place in range(0, query.seq_len, 500):
						
						# parameters to anlayse blast results for a single query
						endo_count = 0
						endo_adj_count = 0
						host_adj_count = 0

						# Go through each hit for a query and determin the organism the hit came from
						hits = []
						for hit in query.hits:
							if place  < hit.hsps[0].hit_end < place + 500:
								hits.append(hit)
								if self.endo in hit.id:
									endo_count += 1
								if "XP_0" in hit.id:
									host_adj_count += 1
								for endo_adj in self.endo_adj:
									if endo_adj in hit.id:
										endo_adj_count += 1
										
						# If the query only has hits from the endosybiont adjacent organisms over the threshold, 
						# report the query as a potential HGT
						if endo_count == 0 and endo_adj_count > 0 and host_adj_count == 0:
							for hit in hits:
								print(dir(hit))
								hit_string = []
								hit_string2 = []
								hit_string.append(str(hit.hsps[0].query_start))
								hit_string.append(str(hit.hsps[0].query_end))
								hit_string.append(str(hit.hsps[0].bitscore))
								hit_string.append(str(hit.hsps[0].evalue))
								hit_string2.append(hit.id)
								hit_string2.append(hit.description)

								print(hit.query_id, file=out_file, end=", ")
								print(", ".join(hit_string), file=out_file, end=", ")

								print(f"{endo_count}, {host_adj_count}, {endo_adj_count}", file=out_file, end=", ")
								total_hits += 1
								contig_num = re.findall(r'\d+', query.id)[0]
								contigs.append(contig_num)

								print(", ".join(hit_string2), file=out_file)


		print(f"total: {total_hits}")

		# Return the HGT genes from the host genome.
		return contigs

	def benchmark(self):
		"""
		Compare results to high confidence published results. (Nowick et. a. 2011)
		"""
		# load known hits
		known_hits = []
		with open("known_hits.txt") as file:
			for line in file:
				known_hits.append(line.strip().split()[0])

		# Populate blast file
		self.blast(self.args.host)

		results = {}
		# Compare results to published HGTs for a varity of thresholds
		for threshold in range(0, 25, 1001):
			self.threshold = threshold
			contigs = self.blast(self.args.host, skip_blast=True)
			contig_count = 0
			for contig in contigs:
				if contig in known_hits:
					print(contig)
					contig_count += 1
			results[threshold] = {"known_hits":contig_count, "total": len(contigs)}

		threshold = list(results.keys())
		known_hit_percent = [val["known_hits"]/len(known_hits) for val in list(results.values())] 
		total_hits = [val["total"] for val in list(results.values())]
		fraction_known_over_total = [val["known_hits"]/(val["total"]+0.001) for val in list(results.values())]
	
		# Print a graph of the data
		plt.style.use('seaborn-poster')

		fig, axs = plt.subplots(3, sharex=True)
		fig.set_size_inches(10, 16)

		axs[0].plot(threshold, total_hits, label="total hits")
		axs[0].set_ylabel("Total hits")
		axs[1].plot(threshold, known_hit_percent, label="Percent of known hits recovered")
		axs[1].set_ylabel("Fraction of known\nHGT events recovered")
		axs[1].set_ylim(0,1)
		axs[2].plot(threshold, fraction_known_over_total, label="Fraction of total hits that are known")
		axs[2].set_ylabel("Fraction of hits that\nare known HGT events")
		axs[2].set_xlabel("Threshold (bitscore)")
		axs[2].set_xlim(0,1000)

		plt.savefig(f"benchmark_.pdf")


	def _download_NCBI_fasta(self, accession, force=False):
		"""
		A private function to download a fasta file from NCBI servers

		Args:
			accession - a string of the NCBI accesssion value to download
			force - download the file, even if it is already downloaded

		Returns:
			fasta_file_name - string location/file name
		"""
		# location/name to save file
		fasta_file_name = f"{self.args.download_location}{accession}.fasta"
		print(f"Trying to get {accession}...", end="\t", flush=True)
		
		# Don't redownload a file if it already exists
		if os.path.exists(fasta_file_name) and not force:
			print(f"{accession} already downloaded")
		else:
			print(f"Downloading {accession}... ", end="\t", flush=True)

			try:
				# Download a genome assembly
				if "GC" in accession[:2]:
					handle = Entrez.elink(dbfrom="assembly", db="nucleotide", id=accession, maxret=1)
					response = Entrez.read(handle)
					
					# Collect nucore ids, download individuall
					nuc_ids = []
					for db in response[0]["LinkSetDb"]:
						if db["LinkName"] == "assembly_nuccore":
							nuc_ids = [nuc_id["Id"] for nuc_id in db["Link"]]
		
					for nuc_accession in nuc_ids:
						time.sleep(1) # be nice to NCBI servers
						print(f"\t\t Downloading {nuc_ids}")
						handle = Entrez.efetch(db="nucleotide", id=nuc_accession, rettype="fasta_cds_aa", maxret=1)
						record = handle.read()

						# append all records to the same file
						with open(fasta_file_name, "a") as file:
							file.write(record)

				# Download a nucore record
				else:
					handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta_cds_aa", maxret=1)
					record = handle.read()
					with open(fasta_file_name, "w") as file:
						file.write(record)
				print(f"Finished")	# Downloading accession

			except Exception as e:
				print(f"\nFailed to download {accession} from NCBI")
				traceback.print_exc()
				sys.exit(1)	# quit on failed download

		return fasta_file_name

if __name__ == "__main__":

	# PARAMETERS
	parser = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)

	parser.add_argument('--download_location',
						help='location to store downloads/database',
						default="resources/",
						type=str,
						required=False)

	# Parser for running blast search
	subparser = parser.add_subparsers(dest='command')   
	blast_parser = subparser.add_parser('blast', help='Perform blast search')
	blast_parser.add_argument('--host',
						help='Host organism cds fasta or NCBI accession to download fasta',
						type=str,
						required=True)
	blast_parser.add_argument('--db',
						help='location of blast database',
						type=str,
						required=False,
						default='EGTfinder.db')
	blast_parser.add_argument("--evalue",
						help='evalue theshold for blast search',
						type=float,
						required=False,
						default=1e-15)
	blast_parser.add_argument('--download_location',
						help='location to store downloads/database',
						default="resources/",
						type=str,
						required=False)

	# Parser for building databsae
	build_db_parser = subparser.add_parser('build_db', help='Build custom blast database to search against')
	build_db_parser.add_argument('--host_adj',
						help='Host-related cds fasta or NCBI accession to download fasta (comma seperated, no space)',
						type=str,
						required=False,
						default="")
	build_db_parser.add_argument('--endo',
						help='Endosymbiont cds fasta or NCBI accession to download fasta',
						type=str,
						required=True)
	build_db_parser.add_argument('--endo_adj',
						help='Endosymbiont-related cds fasta or NCBI accession to \
						download fasta(comma seperated, no space)',
						type=str,
						required=True)
	build_db_parser.add_argument('--db',
						help='location of blast database',
						type=str,
						required=False,
						default='EGTfinder.db')
	build_db_parser.add_argument('--force_download',
						help='Force file to be redownloaded if present',
						type=bool,
						required=False,
						default=False)
	build_db_parser.add_argument('--download_location',
						help='location to store downloads/database',
						default="resources/",
						type=str,
						required=False)

	# Parse arguments
	args = parser.parse_args()


	# If in build database mode, build databse
	if args.command == "build_db":
		finder = EGTfinder(args)
		finder.build_db()

	# if in blast mode, run blast
	elif args.command == "blast":
		finder = EGTfinder(args)
		results = finder.blast(args.host)
		# finder.benchmark()




