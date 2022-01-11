setwd("/data/smartseq2")

library(Seurat)
library(dplyr)

df <- read.delim("counts.txt", sep = "\t", header = T, row.names = 1)
ERCC_genes <- grep(pattern = "^ERCC-", x = rownames(x = df), value = TRUE)
df <- df[!(rownames(df) %in% ERCC_genes), ]
seurat_combined_obj <- CreateSeuratObject(counts = df, project = "KIR")
metadata <- read.csv("metadata.csv", row.names = 1)
seurat_combined_obj <- AddMetaData(seurat_combined_obj, metadata = metadata)

## run seurat
# Calculates the percentage of counts originating from a set of features
seurat_combined_obj[["percent.mt"]] <- PercentageFeatureSet(object = seurat_combined_obj, pattern = "^MT-")

# Visualize QC metrics as a violin plot

VlnPlot(seurat_combined_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  pt.size = .1,
  ncol = 3
)

# Trim off cells based on the vln plot results
seurat_combined_obj <- subset(seurat_combined_obj,
  subset = nFeature_RNA > 800 & nFeature_RNA < 5000 & percent.mt < 15
)

# Normalize data
seurat_combined_obj <- NormalizeData(seurat_combined_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
seurat_combined_obj <- FindVariableFeatures(seurat_combined_obj,
  selection.method = "vst",
  nfeatures = 3000
)

### Scaling the data
## However, particularly for advanced users who would like to use this
## functionality, we strongly recommend the use of our new normalization
## workflow, sctransform. The method is described in our recent preprint,
## with a separate vignette using Seurat v3 here. As with ScaleData, the
## function SCTransform also includes a vars.to.regress parameter.

seurat_combined_obj <- ScaleData(seurat_combined_obj,
  features = rownames(seurat_combined_obj)
)

seurat_combined_obj <- RunPCA(object = seurat_combined_obj, verbose = TRUE)
seurat_combined_obj <- JackStraw(object = seurat_combined_obj, dims = 30, verbose = TRUE)
seurat_combined_obj <- ScoreJackStraw(object = seurat_combined_obj, dims = 1:30)
JackStrawPlot(object = seurat_combined_obj, dims = 1:30)
ElbowPlot(object = seurat_combined_obj, ndims = 30)

## run umap
seurat_combined_obj <- RunUMAP(
  object = seurat_combined_obj,
  dims = 1:20, reduction = "pca",
  min.dist = 0.1
)

seurat_combined_obj <- FindNeighbors(object = seurat_combined_obj, dims = seq(1, 20), verbose = TRUE)
seurat_combined_obj <- FindClusters(
  object = seurat_combined_obj,
  resolution = 0.5,
  verbose = TRUE
)
DimPlot(object = seurat_combined_obj, label = TRUE, reduction = "umap")

# Rename clusters
new.cluster.ids <- c("1", "4", "5", "2", "6", "3")
names(new.cluster.ids) <- levels(seurat_combined_obj)
seurat_combined_obj <- RenameIdents(seurat_combined_obj, new.cluster.ids)

saveRDS(seurat_combined_obj, file = "KIR_pos_ss2.rds")

tiff(filename = "UMAP.tiff", units = "cm", height = 15, width = 15, res = 300)
DimPlot(object = seurat_combined_obj, label = TRUE, reduction = "umap")
dev.off()

diff_markers <- FindAllMarkers(
  object = seurat_combined_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = .5
)

diff_markers <- diff_markers[diff_markers$p_val_adj < 0.05, ]
top_genes <- diff_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_logFC)
write.csv(diff_markers, "diff_markers.csv")
write.csv(top_genes, "top10genes.csv")

tiff(filename = "heatmap_top10genes.tiff", width = 25, height = 50, units = "cm", res = 600)
DoHeatmap(
  object = seurat_combined_obj,
  features = top_genes$gene,
  slot = "scale.data", size = 5
) + NoLegend()
dev.off()
