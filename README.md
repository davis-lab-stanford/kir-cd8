# KIR+ CD8 T Cells

## Installation

```bash
# Install requirements
pip install --upgrade pip wheel
pip install -r requirements.txt
```

Install 10x cellranger 3.1.0 and configure in `env.sh`.

### Create Smartseq reference genome

```bash
# To account for ERCC spike-ins, we create a new reference that includes ERCC
# get human reference from https://www.gencodegenes.org/human/release_21.html
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_21/GRCh38.genome.fa.gz
gunzip GRCh38.genome.fa.gz
cat GRCh38.genome.fa ERCC_spike_in.fasta > /data/smartseq/GRCh38.genome.ERCC.fa

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_21/gencode.v21.annotation.gtf.gz
gunzip gencode.v21.annotation.gtf.gz
cat gencode.v21.annotation.gtf ERCC_spike_in.gtf > /data/smartseq/GRCh38_v21_gene_and_ERCC.gtf

# Make STAR index
# see https://www.biostars.org/p/221781/#221815
STAR --runThreadN 4 --runMode genomeGenerate \
    --genomeDir /data/smartseq/STAR_GRCh38_150bp \
    --genomeFastaFiles /data/smartseq/GRCh38.genome.ERCC.fa \
    --sjdbGTFfile /data/smartseq/GRCh38_v21_gene_and_ERCC.gtf \
    --sjdbOverhang 150;
```

### Create 10x reference

Create 10x reference to be consistent with above:

```bash
# Filter with cellranger mkgtf
# see https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references
source env.sh && cellranger mkgtf \
    /data/10x/GRCh38_v21_gene_and_ERCC.gtf \
    /data/10x/GRCh38_v21_gene_and_ERCC.filtered.gtf \
    --attribute=gene_type:protein_coding \
    --attribute=gene_type:lincRNA \
    --attribute=gene_type:antisense \
    --attribute=gene_type:IG_LV_gene \
    --attribute=gene_type:IG_V_gene \
    --attribute=gene_type:IG_V_pseudogene \
    --attribute=gene_type:IG_D_gene \
    --attribute=gene_type:IG_J_gene \
    --attribute=gene_type:IG_J_pseudogene \
    --attribute=gene_type:IG_C_gene \
    --attribute=gene_type:IG_C_pseudogene \
    --attribute=gene_type:TR_V_gene \
    --attribute=gene_type:TR_V_pseudogene \
    --attribute=gene_type:TR_D_gene \
    --attribute=gene_type:TR_J_gene \
    --attribute=gene_type:TR_J_pseudogene \
    --attribute=gene_type:TR_C_gene;

# Index with cellranger mkref
# see https://kb.10xgenomics.com/hc/en-us/articles/115003616023-Files-and-criteria-used-to-generate-10x-references-
sbatch mkref.sh;
```

## Smartseq2 data (`smartseq/`)

Run RNA-seq transcript quantification:

```bash
mkdir -p /data/smartseq2/logs/
mkdir -p /data/smartseq2/out/

# Start SLURM jobs
# Test first with --dry-run argument. Review proposed jobs with: cat jobs.json | python -m json.tool;
# For usage, see: ./enqueue_smartseq_jobs.py -h
python enqueue_smartseq_jobs.py \
    --log-path='/data/smartseq2/logs/' \
    --output-path='/data/smartseq2/out/' \
    --jobs-json-path='/data/smartseq2/jobs.json' \
    '/data/smartseq2/fastq/';
```

Reshape output into a cells x genes count matrix using `make_count_matrix.R`.

Run UMAP in Seurat, and export cell metadata using `merge_and_umap.R`.

Plot UMAP of all cells with `plot_umap.ipynb`.

## 10x data (`10x/`)

Run cellranger gene expression pipeline:

```bash
mkdir -p /data/10x/logs_cellranger_count/
mkdir -p /data/10x/out_cellranger_count/

# Datasets:
# dataset_peripheral_blood: Pan T cells from blood of healthy controls (HC) and MS patients

# Start SLURM jobs
# Test first with --dry-run argument. Review proposed jobs with: cat jobs_cellranger_count.json | python -m json.tool;
# For usage, see: ./enqueue_10x_jobs.py -h
python enqueue_10x_jobs.py \
    --log-path='/data/10x/logs_cellranger_count/' \
    --output-path='/data/10x/out_cellranger_count/' \
    --jobs-json-path='/data/10x/jobs_cellranger_count.json' \
    '/data/10x/dataset_peripheral_blood';
```

CD8+ T cells (from HC and MS) were identified based on expression of CD8A and CD8B transcripts and integrated with blood CD8+ T cells from COVID-19 patients (ArrayExpress: E-MTAB-9357) using Seurat: `hc_ms.R` and `hc_ms_covid19_integration.R`.

---
## Development

```bash
# Install pre-commit
pre-commit install

# Run lint
make lint
```
