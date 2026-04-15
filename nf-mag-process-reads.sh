## how read QC was performed using the SOP from nf-core/mag
module load nextflow
module load singularity

nextflow run nf-core/mag \
    -resume \
    -profile charlie \
    --input samplesheet.csv \
    --outdir processes-reads \
    --skip_gtdb \
    --skip_spades \
    --skip_megahit \
    --skip_prodigal \
    --skip_prokka \
    --skip_binning \
    -c nextflow.mod.config
