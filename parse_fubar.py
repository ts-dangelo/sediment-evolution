#!/usr/bin/env python3
"""
parse_fubar.py
Parse FUBAR stdout files from a HyPhy analysis directory.

Usage:
    python parse_fubar.py -i ../hyphy/subset_rep_version/fubar/ -o fubar_results.csv
"""

import os
import re
import argparse
import pandas as pd


def parse_fubar_summary(file_path):
    with open(file_path, "r") as f:
        content = f.read()

    if "could not be opened for reading" in content:
        return None

    dN, dS, dN_dS = None, None, None
    for line in content.split("\n"):
        if "* synonymous rate =" in line:
            dS = float(line.split("=")[-1].strip())
        elif "* non-synonymous rate =" in line:
            dN = float(line.split("=")[-1].strip())

    if dN is not None and dS is not None and dS != 0:
        dN_dS = dN / dS

    if "## FUBAR inferred no sites under subject to positive selection" in content:
        positive = "No"
        num_sites = 0
    else:
        positive = "Yes"
        for line in content.split("\n"):
            if "## FUBAR inferred" in line and "sites subject to diversifying positive selection" in line:
                num_sites = int(line.split()[3])
                break
        else:
            num_sites = 0

    len_seq = None
    for line in content.split("\n"):
        if "Loaded a multiple sequence alignment with" in line:
            len_seq = int(line.split("**")[3].split()[0])

    prop_pos = num_sites / len_seq if len_seq else 0

    return {
        "dN": dN,
        "dS": dS,
        "dN/dS": dN_dS,
        "Positive": positive,
        "num_sites": num_sites,
        "len_seq": len_seq,
        "prop_pos": prop_pos
    }


def parse_fubar_codons(file_path):
    with open(file_path, "r") as f:
        content = f.read()

    if "could not be opened for reading" in content:
        return []

    num_seqs = None
    len_seq = None
    for line in content.split("\n"):
        if "Loaded a multiple sequence alignment with" in line:
            match = re.search(r"\*\*(\d+)\*\* sequences, \*\*(\d+)\*\* codons", line)
            if match:
                num_seqs = int(match.group(1))
                len_seq = int(match.group(2))
            break

    codon_data = []
    lines = content.splitlines()
    in_table = False
    for line in lines:
        if re.match(r"\|\s+Codon\s+\|", line):
            in_table = True
            continue
        if in_table:
            if line.strip().startswith("|:"):
                continue
            if not line.strip() or not line.startswith("|"):
                break

            fields = [field.strip() for field in line.strip('| \n').split('|')]
            if len(fields) != 5:
                continue

            try:
                codon = int(fields[0])
                partition = int(fields[1])
                alpha = float(fields[2])
                beta = float(fields[3])
                posterior = float(fields[4].split('=')[-1].strip())

                codon_data.append({
                    "codon": codon,
                    "partition": partition,
                    "alpha": alpha,
                    "beta": beta,
                    "posterior_prob": posterior,
                    "num_seqs": num_seqs,
                    "len_seq": len_seq
                })
            except ValueError:
                continue

    return codon_data


def process_fubar_directory(directory):
    all_rows = []

    for filename in os.listdir(directory):
        if not filename.endswith(".stdout"):
            continue

        file_path = os.path.join(directory, filename)
        gene_cluster = filename.replace(".fubar.stdout", "")

        summary = parse_fubar_summary(file_path)
        if summary is None:
            continue

        if summary["Positive"] == "Yes":
            codon_rows = parse_fubar_codons(file_path)
            for row in codon_rows:
                row.update(summary)
                row["gene_cluster"] = gene_cluster
                all_rows.append(row)
        else:
            row = {
                "gene_cluster": gene_cluster,
                "codon": None,
                "partition": None,
                "alpha": None,
                "beta": None,
                "posterior_prob": None,
                "num_seqs": None,
                "len_seq": summary["len_seq"],
                "dN": summary["dN"],
                "dS": summary["dS"],
                "dN/dS": summary["dN/dS"],
                "Positive": summary["Positive"],
                "num_sites": summary["num_sites"],
                "prop_pos": summary["prop_pos"]
            }
            all_rows.append(row)

    return pd.DataFrame(all_rows)


def main():
    parser = argparse.ArgumentParser(
        description="Parse HyPhy FUBAR stdout files into a summary CSV."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Directory containing .fubar.stdout files"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output CSV file path (e.g. fubar_results.csv)"
    )
    args = parser.parse_args()

    if not os.path.isdir(args.input):
        print(f"[!] Input directory not found: {args.input}")
        raise SystemExit(1)

    print(f"Parsing FUBAR results from: {args.input}")
    df = process_fubar_directory(args.input)

    if df.empty:
        print("[!] No results to write.")
        raise SystemExit(1)

    fubar_all = df[[
        "gene_cluster", "codon", "partition", "alpha", "beta", "posterior_prob",
        "num_seqs", "len_seq", "dN", "dS", "dN/dS",
        "Positive", "num_sites", "prop_pos"
    ]]

    fubar_all.to_csv(args.output, index=False)
    print(f"Results written to: {args.output}")
    print(f"Total gene clusters: {fubar_all['gene_cluster'].nunique()}")
    print(f"Total codon rows:    {len(fubar_all)}")


if __name__ == "__main__":
    main()
