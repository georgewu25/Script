#!/bin/bash -l


echo "=========================================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID : $JOB_ID"
echo "Current job name : $JOB_NAME"
echo "=========================================================="

module load miniconda
conda activate /projectnb/tcwlab/LabMember/mwu/conda_envs/MACS2
conda activate /projectnb/tcwlab/LabMember/mwu/conda_envs/MACS3

###MACS2
for file in "$file_dir"/*.blacklist-filtered.bam; do
    sample=$(basename "$file" | sed 's/\.blacklist-filtered.bam//')
    
    #Convert to BED format
    macs2 randsample -i $bam_dir/${sample}.blacklist-filtered.bam -f BAMPE -p 100 -o $file_dir/${sample}.bed

done


macs2 callpeak -t $file_dir/${sample}.bed -n $sample -f BEDPE -g hs \
  --nomodel --shift -37 --extsize 73 \
  --keep-dup all \
  -q 0.01 \
  -B \
  --broad \
  --cutoff-analysis \
  --outdir $outdir

###MACS3
for file in "$file_dir"/*filtered.bam; do
    sample=$(basename "$file" | sed 's/\.blacklist-filtered.bam//')
    
    #Convert to BED format
    macs3 filterdup --keep-dup all -f BAMPE \
      -i $alignment_dir/${sample}.blacklist-filtered.bam -o $outdir/${sample}.bedpe
    
    macs3 hmmratac -i $outdir/${sample}.bedpe -f BEDPE -n ${sample} --outdir $outdir \
      l 5 -u 50 -c 2 \
      --hmm-type gaussian

done



echo "++++++++++++++++++"
echo "Done"
echo "=========================================================="
echo "Finished JOB #$JOB_ID on : $(date)"
echo "=========================================================="
