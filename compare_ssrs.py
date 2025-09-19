from Bio import SeqIO
from collections import defaultdict
import os

# List of species directories
species_dirs = ["indica", "longipes", "persiciforma", "sylvatica"]
ssr_files = [os.path.join(d, "all_ssr_for_blast.fasta") for d in species_dirs]

ssr_dict = defaultdict(set)  # motif/flank -> set of species
ssr_seqs = defaultdict(list) # motif/flank -> list of (species, record)

for species, ssr_file in zip(species_dirs, ssr_files):
    if not os.path.exists(ssr_file):
        print(f"Warning: {ssr_file} not found.")
        continue
    for record in SeqIO.parse(ssr_file, "fasta"):
        # Use full sequence as key for strict comparison, or motif only for loose
        key = str(record.seq)
        ssr_dict[key].add(species)
        ssr_seqs[key].append((species, record))

num_species = len(species_dirs)
unique = [k for k, v in ssr_dict.items() if len(v) == 1]
common = [k for k, v in ssr_dict.items() if len(v) == num_species]
polymorphic = [k for k, v in ssr_dict.items() if 1 < len(v) < num_species]

print(f"Unique SSRs: {len(unique)}")
print(f"Common SSRs: {len(common)}")
print(f"Polymorphic SSRs: {len(polymorphic)}")

# Optionally, write results to files
with open("unique_ssrs.fasta", "w") as u, open("common_ssrs.fasta", "w") as c, open("polymorphic_ssrs.fasta", "w") as p:
    for k in unique:
        for species, record in ssr_seqs[k]:
            u.write(f">{species}|{record.description}\n{record.seq}\n")
    for k in common:
        for species, record in ssr_seqs[k]:
            c.write(f">{species}|{record.description}\n{record.seq}\n")
    for k in polymorphic:
        for species, record in ssr_seqs[k]:
            p.write(f">{species}|{record.description}\n{record.seq}\n")

print("Results written to unique_ssrs.fasta, common_ssrs.fasta, polymorphic_ssrs.fasta")
