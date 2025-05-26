#!/bin/bash -l


echo "=========================================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID : $JOB_ID"
echo "Current job name : $JOB_NAME"
echo "=========================================================="


module load samtools
module load picard
module load bedtools

genome_index_dir=$main_dir/output/genome_index

aln_sorted_dir=$main_dir/output/alignment
mkdir -p $aln_sorted_dir

#bowtie2-build $GENOME_DIR/Homo_sapiens_assembly38.fasta $genome_index_dir/hg38

SAMPLE="{}"
R1=$fastp_dir/${SAMPLE}_R1.trimmed.fastq.gz
R2=$fastp_dir/${SAMPLE}_R2.trimmed.fastq.gz

bowtie2 --local --very-sensitive --no-mixed --no-discordant -x $genome_index_dir/hg38 \
-1 "$R1" -2 "$R2" | samtools view -bS - > "$aln_sorted_dir/${SAMPLE}.bam"


for file in "$aln_sorted_dir"/*.bam; do
    sample=$(basename "$file" | sed 's/.bam//')
    
    #Sort the output bam files
    samtools sort $aln_dir/${sample}.bam -o $aln_sorted_dir/${sample}_sorted.bam 
    
    #Generate index
    samtools index $aln_sorted_dir/${sample}_sorted.bam

done


sample="{}"

samtools idxstats $aln_sorted_dir/${sample}_sorted.bam > $aln_qc_dir/${sample}_sorted.idxstats
samtools flagstat $aln_sorted_dir/${sample}_sorted.bam > $aln_qc_dir/${sample}_sorted.flagstat

#Remove mitochondrial reads
samtools view -h $aln_sorted_dir/${sample}_sorted.bam | grep -v chrM | \
samtools sort -O bam -o $aln_qc_dir/${sample}.rmChrM.bam -T $aln_sorted_dir

#Remove duplicated reads
picard MarkDuplicates QUIET=true \
  INPUT=$aln_qc_dir/${sample}.rmChrM.bam \
  OUTPUT=$aln_qc_dir/${sample}.marked.bam \
  METRICS_FILE=$aln_qc_dir/${sample}.dup.metrics \
  REMOVE_DUPLICATES=false \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=LENIENT \
  TMP_DIR= $aln_qc_dir
  
rm $aln_qc_dir/${sample}.rmChrM.bam

# -f 2 retains only mapped pairs
# -F 1548 exclues unmapped and duplicated reads
# -q 30 mapping score threshold of 30
samtools view -h -b -f 2 -F 1548 -q 30 $aln_qc_dir/${sample}.marked.bam | samtools sort -o $aln_qc_dir/${sample}.filtered.bam

rm $aln_qc_dir/${sample}.marked.bam

#Remove reads within the blacklist regions
bedtools intersect -nonamecheck -v -abam $aln_qc_dir/${sample}.filtered.bam -b /projectnb/tcwlab/RefData/ENCODE/Blacklist/lists/hg38-blacklist.v2.bed > $aln_qc_dir/${sample}.tmp.bam

#Sort and index the bam file
samtools sort -O bam -o $aln_qc_dir/${sample}.blacklist-filtered.bam $aln_qc_dir/${sample}.tmp.bam
samtools index $aln_qc_dir/${sample}.blacklist-filtered.bam

rm $aln_qc_dir/${sample}.filtered.bam
rm $aln_qc_dir/${sample}.tmp.bam


multiqc $alignment_dir -o $outdir


for file in "$aln_sorted_dir"/*filtered.bam; do
    sample=$(basename "$file" | sed 's/\.blacklist-filtered.bam//')
    
    # -AS Assume Sorted
    # --INCLUDE_DUPLICATES 
    # -M Read percentage threshold for FR, RF, and Tandem
    
    picard CollectInsertSizeMetrics -I $aln_sorted_dir/${sample}.blacklist-filtered.bam \
    -AS true \
    --INCLUDE_DUPLICATES false \
    -M 0.05 \
    -O $outdir/${sample}_size_metrics.txt \
    -H $outdir/${sample}_size_histogram.pdf
    
done


echo "++++++++++++++++++"
echo "Done"
echo "=========================================================="
echo "Finished JOB #$JOB_ID on : $(date)"
echo "=========================================================="
