from Bio import SeqIO

# Input and output file names
input_file = "all_ssr_flanks.fasta"       # your current FASTA
output_file = "all_ssr_for_blast.fasta"   # merged file for BLAST

# Read all FASTA records
records = list(SeqIO.parse(input_file, "fasta"))

with open(output_file, "w") as out:
    # Go through sequences in groups of 3 (upstream, SSR, downstream)
    for i in range(0, len(records), 3):
        if i + 2 < len(records):  # make sure we have a full set
            upstream = str(records[i].seq)
            ssr = str(records[i+1].seq)
            downstream = str(records[i+2].seq)

            # Extract SSR number from middle header
            header = records[i+1].description
            # Example: "NC_061916.1_Corydalis_trisecta_chloroplast,_complete_genome_SSR2_(A)12_6368_6379"
            # We just keep it as is + "_full"
            new_header = header.replace(",", "_").replace(" ", "_") + "_full"

            # Merge sequences
            combined_seq = upstream + ssr + downstream

            # Write to FASTA (single line sequence)
            out.write(f">{new_header}\n{combined_seq}\n")

print(f"✅ Merged file created: {output_file}")
