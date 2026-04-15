#!/usr/bin/env python3
import pysam
import re
import sys

def calculate_identity_from_md(read):
    """Calculate percent identity using MD tag"""
    if not read.has_tag('MD'):
        return None
    
    md = read.get_tag('MD')
    
    matches = sum(int(x) for x in re.findall(r'\d+', md))
    aligned_length = read.query_alignment_length
    
    if aligned_length > 0:
        return (matches / aligned_length) * 100
    return 0

def filter_bam_by_identity(input_bam, output_bam, min_identity=95):
  
    infile = pysam.AlignmentFile(input_bam, "rb", check_sq=False)
    outfile = pysam.AlignmentFile(output_bam, "wb", template=infile)
    
    total = 0
    kept = 0
    
    for read in infile:
        total += 1
        if read.is_unmapped:
            continue
            
        identity = calculate_identity_from_md(read)
        
        if identity is not None and identity >= min_identity:
            outfile.write(read)
            kept += 1
    
    infile.close()
    outfile.close()
    
    print(f"Total reads: {total}")
    print(f"Kept reads: {kept} ({kept/total*100:.2f}%)")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python filter_bam_identity.py input.bam output.bam [min_identity]")
        sys.exit(1)
    
    input_bam = sys.argv[1]
    output_bam = sys.argv[2]
    min_identity = float(sys.argv[3]) if len(sys.argv) > 3 else 95.0
    
    filter_bam_by_identity(input_bam, output_bam, min_identity)
