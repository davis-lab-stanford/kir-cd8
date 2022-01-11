#!/bin/sh

#SBATCH --account=mmdavis
#SBATCH --time=2-00:00:00
#SBATCH --job-name=STAR_align
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G

# set input_file1, input_file2, wDir in outside wrapper script

module load samtools/1.4;
module load python/2.7.3;
module load rsem/1.2.31;
module load STAR/2.6.1d;
module load gcc/9.2.0-centos_7;

base_name=$(basename "$input_file1")
echo ${base_name%%_R1_001.fastq.gz}

fa_file=/data/smartseq/GRCh38.genome.ERCC.fa
GTF_file=/data/smartseq/GRCh38_v21_gene_and_ERCC.gtf

gunzip -c $input_file1 > $wDir/Input_R1_001.fastq
gunzip -c $input_file2 > $wDir/Input_R2_001.fastq

## Start STAR alignment
STAR --genomeDir /data/smartseq/STAR_GRCh38_150bp --readFilesIn $wDir/Input_R1_001.fastq $wDir/Input_R2_001.fastq --runThreadN 16 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outReadsUnmapped Fastx

## get un-sorted bam file
samtools view -bT $fa_file $wDir/Aligned.out.sam > $wDir/Aligned_out.bam

## sort the bam file by read names for htseq-count
samtools sort -O bam -n -T $wDir/accepted_hits_ -o $wDir/accepted_hits_sorted.bam $wDir/Aligned_out.bam

## get sam file from read-name-sorted bam file
samtools view $wDir/accepted_hits_sorted.bam > $wDir/accepted_hits_sorted.sam

## Create counts
./htseq_count.sh --stranded=no --type=exon --idattr=gene_name --mode=intersection-nonempty $wDir/accepted_hits_sorted.sam $GTF_file > $wDir/output_counts_transcript_id.txt

rm $wDir/*.fastq
rm $wDir/*sam
rm $wDir/Aligned_out.bam
rm $wDir/accepted_hits_sorted.bam
rm -r $wDir/_STARtmp
rm $wDir/Unmapped.out.mate*
