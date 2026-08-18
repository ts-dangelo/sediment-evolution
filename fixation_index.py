import pandas as pd
import numpy as np
from Bio import SeqIO
from collections import Counter
import os
from scipy.stats import entropy

def read_fasta_to_dict(fasta_file):
    """read in aligned amino acid fastas"""
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        lib_id = record.id.split('.')[0]
        sequences[lib_id] = str(record.seq)
    return sequences

def get_aa_at_position(sequence, position):
    """get amino acid sequence at given position (hyphy results are 1-indexed)"""
    idx = position - 1
    if idx >= len(sequence):
        return None
    aa = sequence[idx]
    # Ignore gaps
    if aa == '-':
        return None
    return aa

def calculate_fst(aa_by_bioproject):
    """
    From methods:

    For codons under diversifying selection, the geographic distribution of individual alleles was summarized using the fixation index (Fst) 30,31. 
    
    Ht = 1 - sum_i (n_i / N)^2: where Ht is the total heterozygosity of the alleles in a given codon, where ni is the number of a given amino acid allele across the total of all samples (N)

    h_k = 1 - sum_i (n_ik / N_k)^2

    Hs = sum_k (N_k / N) * H_k

    Fst = (H_T - H_S) / H_T 
    
     where Ht is the total heterozygosity of the alleles in a given codon, where ni is the number of a given amino acid allele across the total of all samples (N) Hk is the within-location (sampling site i.e BioProject) heterozygosity where nik and Nk are allele and total counts within a given site and Hs is the average within-location heterozygosity of each allele in a given codon, normalized by the relative proportion of a given allele to the total sequence number in the codon (N). The fixation index (Fst) is on a 0-1 scale, with 1 meaning complete geographic segregation of alleles and 0 meaning no geographic signal. We report data for codons under diversifying selection that are present at a minimum of four of the nine global sites in this analysis (Supplemental Figure 8).

    """
    all_aas = set()
    for aas in aa_by_bioproject.values():
        all_aas.update(aas)

    # Calculate total frequencies (across all populations)
    total_counts = Counter()
    total_n = 0
    for aas in aa_by_bioproject.values():
        total_counts.update(aas)
        total_n += len(aas)

    if total_n == 0:
        return None

    # Calculate Ht (total heterozygosity)
    Ht = 1 - sum((count / total_n) ** 2 for count in total_counts.values())

    # Calculate Hs (weighted by population size)
    Hs = 0
    for aas in aa_by_bioproject.values():
        n_k = len(aas)
        if n_k == 0:
            continue
        counts = Counter(aas)
        h_k = 1 - sum((count / n_k) ** 2 for count in counts.values()) # within-location (bioproject) heterozygosity (sum = given allele, n_k = total aa counts within a given site)
        Hs += (n_k / total_n) * h_k   # proportional weighting part (dividing by total_n)

    # Calculate FST
    if Ht == 0:
        return 0

    fst = (Ht - Hs) / Ht
    return fst

def calculate_shannon_evenness(amino_acids):
    """Calculate Shannon evenness - Not used in the publication"""
    if len(amino_acids) == 0:
        return None
    
    counts = Counter(amino_acids)
    n_species = len(counts)
    
    if n_species <= 1:
        return 1.0  # Perfect evenness when only one species
    
    # Calculate Shannon entropy
    proportions = np.array(list(counts.values())) / len(amino_acids)
    shannon_entropy = entropy(proportions, base=np.e)
    
    # Normalize by max possible entropy (log of number of species)
    max_entropy = np.log(n_species)
    evenness = shannon_entropy / max_entropy if max_entropy > 0 else 0
    
    return evenness

def analyze_diversifying_codons(df, lib_to_bioproject, fasta_dir, min_bioprojects=4):
    """
    Analyze codons under diversifying selection
    Only includes positions where sequences come from at least min_bioprojects bioprojects
    Gaps in alignments not analyzed 
    """
    results = []
    
    # Filter for diversifying selection
    diversifying = df[df['class'] == 'Diversifying'].copy()
    
    for idx, row in diversifying.iterrows():
        gene_cluster = row['gene_cluster']
        codon_pos = int(row['codon'])
        
        # Construct fasta filename
        fasta_file = os.path.join(fasta_dir, f"{gene_cluster}.aln")
        
        if not os.path.exists(fasta_file):
            print(f"Warning: {fasta_file} not found")
            continue
        
        # Read sequences
        sequences = read_fasta_to_dict(fasta_file)
        
        # Extract amino acids at this position for each bioproject
        aa_by_bioproject = {}
        all_amino_acids = []
        
        for lib_id, seq in sequences.items():
            # Get bioproject for this library
            bioproject = lib_to_bioproject.get(lib_id)
            if not pd.notna(bioproject):  # ← CHANGED: Now filters out both None and nan
                continue
            
            # Get amino acid at position
            aa = get_aa_at_position(seq, codon_pos)
            if aa is None:  # Skip gaps and invalid positions
                continue
            
            # Store by bioproject
            if bioproject not in aa_by_bioproject:
                aa_by_bioproject[bioproject] = []
            aa_by_bioproject[bioproject].append(aa)
            all_amino_acids.append(aa)
        
        # Filter: require at least min_bioprojects
        if len(aa_by_bioproject) < min_bioprojects:
            continue
        
        if len(all_amino_acids) == 0:
            continue
        
        # Calculate statistics
        fst = calculate_fst(aa_by_bioproject)
        
        # Unique residues
        unique_residues = len(set(all_amino_acids))
        
        # Most frequent residue percentage
        aa_counts = Counter(all_amino_acids)
        most_common_aa, most_common_count = aa_counts.most_common(1)[0]
        most_frequent_pct = (most_common_count / len(all_amino_acids)) * 100
        
        # Shannon evenness
        shannon_even = calculate_shannon_evenness(all_amino_acids)
        
        results.append({
            'gene_cluster': gene_cluster,
            'codon': codon_pos,
            'fst': fst,
            'unique_residues': unique_residues,
            'most_frequent_residue': most_common_aa,
            'most_frequent_pct': most_frequent_pct,
            'shannon_evenness': shannon_even,
            'n_sequences': len(all_amino_acids),
            'n_bioprojects': len(aa_by_bioproject)
        })
    
    return pd.DataFrame(results)


fasta_dir = "/directory/with/single/copy/gene/cluster/faa/files'
results_df = analyze_diversifying_codons(fel_annots_df, lib_to_bioproject, fasta_dir, min_bioprojects=1) # data plotted is min_bioprojects=3
# lib_to_bioproject is a dictionary linking fasta header names to Bioprojects
# fel_annots_df is a dataframe of the Fixed Effects Likelihood results with the associated annotation data 
results_df.to_csv("diversifying_codon_statistics-weighted-All.csv")
