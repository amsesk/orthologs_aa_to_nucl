import sys
import os
import pandas as pd
import argparse

class blastdb_maker(object):
	def __init__(self, marker, ltp_list, prepend):
		self.marker = marker
		self.ltp_list = ltp_list
		self.genome_fasta_prepend_path = prepend

	def canonicalize_fasta_paths(self):
		return [os.path.join(self.genome_fasta_prepend_path, l) for l in self.ltp_list]

	def to_shell(self):
		full_paths = self.canonicalize_fasta_paths()
		print(f"cat {' '.join(full_paths)} > {marker}.missing.fasta")
		print(f"makeblastdb -in {marker}.missing.fasta -out {marker}.missing.fasta -title {marker}.missing.fasta -dbtype nucl")

	def title(self):
		return f"{marker}.missing.fasta"

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--isolates', help = "Path to isolate spread sheet with `LTP` and `genome_fasta` columns.")
	parser.add_argument('-o', '--ortholog_occupancy', help = "Path to three-column ortholog occupancy output from `orthologs_missing_from.R`.")
	parser.add_argument('-g', '--genome_dir', help = "Path to directory containing the genome fastas for the taxon set.")
	parser.add_argument('-f', '--ortholog_fasta_dir', help = "Path to directory containing to ortholog nucleotide fastas output by main.")
	args = parser.parse_args()

	ltp_to_asm = pd.read_excel(args.isolates)[['LTP', 'genome_fasta']].set_index('LTP').to_dict()['genome_fasta']

	with open(args.ortholog_occupancy) as occ:
		occ.readline()
		for line in occ:
			spl = [x.strip() for x in line.split('\t')]

			marker = spl[0]
			missing = spl[1].split(',')
			present = spl[2].split(',')

			if len(missing) == 1 and len(missing[0]) == 0:
				continue

			m = blastdb_maker(marker, [ltp_to_asm[ltp] for ltp in missing], args.genome_dir)
			m.to_shell()

			print(f"blastn -query {os.path.join(args.ortholog_fasta_dir, marker)}.fasta -db {m.title()} -evalue 1e-50 -num_threads 1 -outfmt 6 -out {marker}.fasta.blastout -max_target_seqs 1")

