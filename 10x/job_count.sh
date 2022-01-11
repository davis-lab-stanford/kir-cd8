#!/usr/bin/env bash
#SBATCH --job-name=cellranger-count
#SBATCH --export=ALL
#SBATCH --account=mmdavis
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --time=2-00:00:00
#SBATCH --mem=128G
#SBATCH --qos=normal

. /data/10x/env.sh && cellranger count \
    --id="$sample_name" \
    --transcriptome=/data/10x/GRCh38_genome_ERCC_reference_genome \
    --fastqs="$fastq_path" \
    --sample="$sample_name" \
    --jobmode=local --localcores=16 --localmem=115;
