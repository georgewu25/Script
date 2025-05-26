#!/bin/bash -l

echo "=========================================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID : $JOB_ID"
echo "Current job name : $JOB_NAME"
echo "=========================================================="



module load fastqc
module load fastp
module load multiqc

for fq in $raw_data_dir/*.fastq.gz; do
    fastqc -o $fastqc_dir "$fq"
done


for R1 in $raw_data_dir/*_R1_001.fastq.gz; do
    SAMPLE=$(basename $R1 | sed 's/_R1_001.fastq.gz//')

    # Define the corresponding R2 file
    R2="$raw_data_dir/${SAMPLE}_R2_001.fastq.gz"

    fastp -i "$R1" -I "$R2" -o "$fastqc_trimmed_dir/${SAMPLE}_R1.trimmed.fastq.gz" -O "$fastqc_trimmed_dir/${SAMPLE}_R2.trimmed.fastq.gz" \
        --detect_adapter_for_pe -j "$fastqc_trimmed_dir/${SAMPLE}.fastp.json" -h "$fastqc_trimmed_dir/${SAMPLE}.fastp.html"
        
done


fqc_file_count=`ls $fastqc_dir/*.html |wc -l`
echo "TOTAL FASTQC File Count=$fqc_file_count"

multiqc $fastqc_dir -o $outdir


fqc_file_count=`ls $fastqc_trimmed_dir/*.html |wc -l`
echo "TOTAL FASTQC File Count=$fqc_file_count"

multiqc $fastqc_trimmed_dir -o $outdir



echo "++++++++++++++++++"
echo "Done"
echo "=========================================================="
echo "Finished JOB #$JOB_ID on : $(date)"
echo "=========================================================="
