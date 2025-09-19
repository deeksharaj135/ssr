from Bio import SeqIO

# Input files
fasta_file = "indica.fasta"
misa_file = "indica.fasta.misa"

# Read the genome
record = SeqIO.read(fasta_file, "fasta")
sequence = record.seq
seq_length = len(sequence)

# Open output file
with open("all_ssr_flanks.fasta", "w") as out:
    with open(misa_file) as f:
        header = next(f)  # skip header line

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue  # skip bad lines

            # Parse columns (MISA format: ID, SSR_nr, SSR_type, Start, End, Motif, Repeats…)
            seq_name = parts[0]
            ssr_number = parts[1]
            ssr_type = parts[2]
            start = int(parts[5])
            end = int(parts[6])
            motif = parts[3]

            # Define flanking regions
            up_start = max(1, start - 200)
            up_end = start - 1
            down_start = end + 1
            down_end = min(seq_length, end + 200)

            # Extract sequences
            upstream = sequence[up_start-1:up_end] if up_end >= up_start else ""
            ssr_seq = sequence[start-1:end]
            downstream = sequence[down_start-1:down_end] if down_end >= down_start else ""

            # Write to FASTA
            out.write(f">{seq_name}_SSR{ssr_number}_{motif}_{start}_{end}_upstream200\n{upstream}\n")
            out.write(f">{seq_name}_SSR{ssr_number}_{motif}_{start}_{end}\n{ssr_seq}\n")
            out.write(f">{seq_name}_SSR{ssr_number}_{motif}_{start}_{end}_downstream200\n{downstream}\n")

print("✅ All SSRs with flanking regions saved in all_ssr_flanks.fasta")
