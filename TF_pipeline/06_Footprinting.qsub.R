#!/bin/bash -l


echo "=========================================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID : $JOB_ID"
echo "Current job name : $JOB_NAME"
echo "=========================================================="


for file in "$file_dir"/*.shifted.sorted.bam; do
    sample=$(basename "$file" | sed 's/.shifted.sorted.bam//')
    
    samtools view -h $bam_dir/${sample}.shifted.sorted.bam \
    $(cut -f1 $main_dir/output/chroms_to_keep.txt) | samtools view \
    -b -o $output_dir/${sample}_filtered.bam

    samtools index $output_dir/${sample}_filtered.bam
    

done

for file in "$file_dir"/*_peaks.broadPeak; do
    sample=$(basename "$file" | sed 's/_peaks.broadPeak//')
    
    awk 'NR==FNR {chrom[$1]; next} $1 in chrom' $main_dir/output/chroms_to_keep.txt $bed_dir/${sample}_peaks.broadPeak > $output_dir/${sample}_filtered.bed


done


####RGT-HINT
samples=()

for file in "$file_dir"/*_merged.bam; do
    sample_name=$(basename "$file" | sed 's/_merged.bam//')
    
    samples+=("$sample_name")
done

sample=${samples[($SGE_TASK_ID-1)]}


rgt-hint footprinting --atac-seq --paired-end --output-location=$footprint_outdir --output-prefix=${sample} \
  $file_dir/${sample}_merged.bam \
  $file_dir/${sample}_merged.bed \
  --organism=hg38


rgt-motifanalysis matching --organism=hg38 --input-files $(ls $footprint_dir/*.bed)



####TOBIAS
TOBIAS ATACorrect --bam $bam_dir/${sample}.blacklist-filtered.bam \
      --genome /projectnb/tcwlab/RawData/Genome/hg38/Homo_sapiens_assembly38.fasta \
      --peaks $bed_dir/${sample}.bed \
      --blacklist /projectnb/tcwlab/RefData/ENCODE/Blacklist/lists/hg38-blacklist.v2.bed \
      --outdir $output_dir --cores 8

TOBIAS FootprintScores --signal $bw_dir/${sample}.blacklist-filtered_corrected.bw \
  --regions $bed_dir/${sample}.bed \
  --output $output_dir/${sample}_footprints.bw \
  --cores 8


echo "++++++++++++++++++"
echo "Done"
echo "=========================================================="
echo "Finished JOB #$JOB_ID on : $(date)"
echo "=========================================================="