!/usr/bin/env bash

BAM_DIR="./bams/"
OUTFILE="bam_mean_depth_coverage.tsv"

if [[ ! -f "$OUTFILE" ]]; then
    echo -e "filename\tmean_depth\tpercent_genome_covered" > "$OUTFILE"
fi

declare -A DONE
tail -n +2 "$OUTFILE" | awk '{print $1}' | while read -r fname; do
    DONE["$fname"]=1
done

for bam in "$BAM_DIR"/*.bam; do
    [[ -e "$bam" ]] || continue

    fname=$(basename "$bam")

    if [[ ${DONE["$fname"]+_} ]]; then
        echo "Skipping $fname (already done)"
        continue
    fi

    read mean_depth percent <<< $(
        samtools depth -a "$bam" | \
        awk '{
                sum += $3
                total++
                if ($3 > 0) covered++
             }
             END {
                if (total > 0)
                    printf "%.6f %.4f", sum/total, (covered/total)*100
                else
                    print "0 0"
             }'
    )

    echo -e "$fname\t$mean_depth\t$percent" >> "$OUTFILE"
done
