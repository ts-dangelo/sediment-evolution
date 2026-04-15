#!/usr/bin/env python3

import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


nucleotide_file = "/mnt/s1/labs/lindsay/epc/prodigal_mmseqs/all_genes.ffn"
core_alns_dir   = "/mnt/s1/labs/lindsay/epc/prodigal_mmseqs/aero_misag_gene_faas_sc"
core_ffns_dir   = "/mnt/s1/labs/lindsay/epc/prodigal_mmseqs/aero_misag_gene_faas_sc_ffn"

os.makedirs(core_ffns_dir, exist_ok=True)

print("Loading nucleotide sequences...")
nucleotide_sequences = {
    record.id: record.seq
    for record in SeqIO.parse(nucleotide_file, "fasta")
}
print(f"Loaded {len(nucleotide_sequences)} nucleotide sequences")

for aln_file in sorted(os.listdir(core_alns_dir)):
    if not aln_file.endswith(".faa"):
        continue

    aln_path = os.path.join(core_alns_dir, aln_file)
    output_ffn_path = os.path.join(
        core_ffns_dir,
        aln_file.replace(".faa", ".ffn")
    )

    amino_acid_records = list(SeqIO.parse(aln_path, "fasta"))
    updated_records = []

    for record in amino_acid_records:
        if record.id in nucleotide_sequences:
            updated_records.append(
                SeqRecord(
                    nucleotide_sequences[record.id],
                    id=record.id,
                    description=""
                )
            )
        else:
            protein_length = len(record.seq)
            updated_records.append(
                SeqRecord(
                    "-" * protein_length,
                    id=record.id,
                    description="Missing nucleotide sequence"
                )
            )
            print(
                f"WARNING: {record.id} not found in nucleotide DB "
                f"(alignment: {aln_file})"
            )

    with open(output_ffn_path, "w") as out_handle:
        SeqIO.write(updated_records, out_handle, "fasta")

print("Done creating FFN files.")
