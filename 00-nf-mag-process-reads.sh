## how read QC was performed using the SOP from nf-core/mag
module load nextflow
module load singularity

nextflow run nf-core/mag \
    -resume \
    -profile charlie \
    --input samplesheet.csv \
    --outdir processed_reads \
    --skip_gtdb \
    --skip_spades \
    --skip_megahit \
    --skip_prodigal \
    --skip_prokka \
    --skip_binning


## pull the reads from the nf work dir and rename

mkdir -p qc_reads

find ./work \
  -type f \
  -name '*_run0_phix_removed.unmapped_*fastq.gz' \
  -exec cp -t ./qc_reads/. {} +

shopt -s nullglob

for f in ./qc_reads/*.fastq.gz; do
    base=$(basename "$f")

    if [[ "$base" =~ (R1|_1)\.fastq\.gz$ ]]; then
        read="_1"
    elif [[ "$base" =~ (R2|_2)\.fastq\.gz$ ]]; then
        read="_2"
    else
        continue
    fi

    stripped=${base/_R1/}
    stripped=${stripped/_R2/}
    stripped=${stripped/_1.fastq.gz/}
    stripped=${stripped/_2.fastq.gz/}
    stripped=${stripped%.fastq.gz}

    stripped=${stripped/_run0_phix_removed.unmapped/_QC}

    newname="qc_reads/${stripped}${read}.fastq.gz"

    if [[ "$f" != "$newname" ]]; then
        echo "  $base -> $(basename "$newname")"
        mv -n "$f" "$newname"
    fi
done
