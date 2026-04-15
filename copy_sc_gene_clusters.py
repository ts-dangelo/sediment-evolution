#!/usr/bin/env python3

import os
import sys
import shutil
from Bio import SeqIO

def genome_from_header(header):
    return header.split('.')[0]

def is_single_copy(fasta_path):
    genomes = []
    for rec in SeqIO.parse(fasta_path, "fasta"):
        genomes.append(genome_from_header(rec.id))
    return len(genomes) == len(set(genomes))

def filter_single_copy(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    for fname in os.listdir(input_dir):
        if not fname.endswith(".faa"):
            continue

        input_path = os.path.join(input_dir, fname)

        if is_single_copy(input_path):
            shutil.copy(input_path, os.path.join(output_dir, fname))
            print(f"KEEP: {fname}")
        else:
            print(f"SKIP: {fname}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python filter_single_copy_gene_clusters.py <input_dir> <output_dir>")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_dir = sys.argv[2]

    filter_single_copy(input_dir, output_dir)
