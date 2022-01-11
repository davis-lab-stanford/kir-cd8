#!/usr/bin/env bash
#SBATCH --job-name=cellranger-mkref
#SBATCH --account=mmdavis
#SBATCH --export=ALL
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --time=2-00:00:00
#SBATCH --mem=40G
#SBATCH --qos=normal
#SBATCH -o /data/10x/cellranger_mkref_job.out
#SBATCH -e /data/10x/cellranger_mkref_job.err

# run this with sbatch

cd /data/10x/ && . env.sh && cellranger mkref \
    --genome=GRCh38_genome_ERCC_reference_genome \
    --fasta=/data/10x/GRCh38.genome.ERCC.fa \
    --genes=/data/10x/GRCh38_v21_gene_and_ERCC.filtered.gtf;
