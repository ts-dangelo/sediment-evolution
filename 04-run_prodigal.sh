#!/bin/bash

set -euo pipefail

module use /mod/scgc/
module load prodigal/2.6.2
module load parallel/20180722

NUM_JOBS=12

WORKDIR="/home/tdangelo/epc/prodigal_mmseqs/"
INPUT_DIR="/home/tdangelo/epc/mmdb_analysis/misag_sags_mmdb951kf_mags/" ## Directory containing all MIMAG SAGs and the MAGs produced

FAA_DIR="${WORKDIR}/aero_misag_gene_faas"
FFN_DIR="${WORKDIR}/aero_misag_gene_ffns"
GFF_DIR="${WORKDIR}/aero_misag_gene_gffs"
GBK_DIR="${WORKDIR}/aero_misag_gene_gbks"
LOG_DIR="${WORKDIR}/aero_misag_gene_logs"

mkdir -p "$FAA_DIR" "$FFN_DIR" "$GFF_DIR" "$GBK_DIR" "$LOG_DIR"

cd "$WORKDIR"

run_prodigal() {
    fna_file="$1"

    base=$(basename "$fna_file")
    base="${base%%.*}"

    prodigal \
        -i "$fna_file" \
        -a "${FAA_DIR}/${base}.faa" \
        -d "${FFN_DIR}/${base}.ffn" \
        -o "${GFF_DIR}/${base}.gff" \
        -p meta \
        -q \
        > "${LOG_DIR}/${base}.log" 2>&1

    prodigal \
        -i "$fna_file" \
        -f gbk \
        -p meta \
        -q \
        > "${GBK_DIR}/${base}.gbk"
}


export -f run_prodigal
export FAA_DIR FFN_DIR GFF_DIR GBK_DIR LOG_DIR

find "$INPUT_DIR" -maxdepth 1 \( -name "*.fna" -o -name "*.fasta" -o -name "*.fa" \) \
| parallel -j "$NUM_JOBS" run_prodigal


echo "Job finished at $(date)"
