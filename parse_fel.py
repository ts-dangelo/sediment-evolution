#!/usr/bin/env python3
"""
parse_fel.py
Parse FEL stdout files from a HyPhy analysis directory.

Usage:
    python parse_fel.py -i ../hyphy/subset_rep_version/fel/ -o fel_results.csv
"""

import os
import re
import argparse
import pandas as pd


def parse_fel_stdout(stdout_path):
    results = []

    num_seqs = None
    len_seq = None
    dnds_ratio = None
    gene_cluster = os.path.basename(stdout_path).split(".")[0]

    with open(stdout_path) as f:
        lines = f.readlines()

    # Pass 1: metadata
    for line in lines:
        if "Loaded a multiple sequence alignment" in line:
            clean_line = re.sub(r"\*+", "", line)
            stats_match = re.search(r'(\d+)\s+sequences,\s+(\d+)\s+codons', clean_line)
            if stats_match:
                num_seqs = int(stats_match.group(1))
                len_seq = int(stats_match.group(2))

            gc_match = re.search(r'/([^/]+)\.codon\.fna', line)
            if gc_match:
                gene_cluster = gc_match.group(1)

        if "non-synonymous/synonymous rate ratio for" in line:
            dnds_match = re.search(r'=\s*([\d.]+)', line)
            if dnds_match:
                dnds_ratio = float(dnds_match.group(1))

    # Pass 2: codon table
    in_table = False
    parsed_codons = 0
    n_div = 0
    n_pur = 0

    for line in lines:
        if line.strip().startswith("|") and "Codon" in line:
            in_table = True
            continue

        if in_table and line.strip().startswith("|") and "Selection detected?" not in line:
            parts = [x.strip() for x in line.strip().split("|")[1:-1]]
            if len(parts) != 6:
                continue

            codon, partition, alpha, beta, lrt, sel_result = parts

            p_match = re.search(r"p\s*=\s*([\d.]+)", sel_result)
            if not p_match:
                continue

            p_value = float(p_match.group(1))

            if "Pos" in sel_result:
                site_class = "Diversifying"
                n_div += 1
            elif "Neg" in sel_result:
                site_class = "Purifying"
                n_pur += 1
            else:
                site_class = "Neutral"

            results.append({
                "gene_cluster": gene_cluster,
                "codon": int(codon),
                "partition": int(partition),
                "alpha": float(alpha),
                "beta": float(beta),
                "LRT": float(lrt),
                "p-value": p_value,
                "class": site_class,
                "num_seqs": num_seqs,
                "len_seq": len_seq,
                "dN/dS": dnds_ratio
            })
            parsed_codons += 1

        elif in_table and not line.strip().startswith("|"):
            break

    prop_div = n_div / len_seq if len_seq else 0
    prop_pur = n_pur / len_seq if len_seq else 0

    if parsed_codons == 0:
        return pd.DataFrame([{
            "gene_cluster": gene_cluster,
            "codon": None,
            "partition": None,
            "alpha": None,
            "beta": None,
            "LRT": None,
            "p-value": None,
            "class": "NA",
            "num_seqs": num_seqs,
            "len_seq": len_seq,
            "dN/dS": dnds_ratio,
            "prop_div": prop_div,
            "prop_pur": prop_pur
        }])

    df = pd.DataFrame(results)
    df["prop_div"] = prop_div
    df["prop_pur"] = prop_pur

    return df


def parse_fel_directory(directory):
    all_results = []

    json_files = [f for f in os.listdir(directory) if f.endswith(".fel.json")]
    gene_clusters_with_json = set(f.split(".")[0] for f in json_files)

    for file in os.listdir(directory):
        if file.endswith(".fel.stdout"):
            gene_cluster_name = file.split(".")[0]
            if gene_cluster_name in gene_clusters_with_json:
                full_path = os.path.join(directory, file)
                df = parse_fel_stdout(full_path)
                if not df.empty:
                    all_results.append(df)

    if all_results:
        return pd.concat(all_results, ignore_index=True)
    else:
        print("[!] No valid selection results found.")
        return pd.DataFrame()


def main():
    parser = argparse.ArgumentParser(
        description="Parse HyPhy FEL stdout files into a summary CSV."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Directory containing .fel.stdout and .fel.json files"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output CSV file path (e.g. fel_results.csv)"
    )
    args = parser.parse_args()

    if not os.path.isdir(args.input):
        print(f"[!] Input directory not found: {args.input}")
        raise SystemExit(1)

    print(f"Parsing FEL results from: {args.input}")
    fel_all = parse_fel_directory(args.input)

    if fel_all.empty:
        print("[!] No results to write.")
        raise SystemExit(1)

    fel_all.to_csv(args.output, index=False)
    print(f"Results written to: {args.output}")
    print(f"Total gene clusters: {fel_all['gene_cluster'].nunique()}")
    print(f"Total codon rows: {len(fel_all)}")


if __name__ == "__main__":
    main()
