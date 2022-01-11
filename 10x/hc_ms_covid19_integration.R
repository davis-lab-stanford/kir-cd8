setwd("/data/10x")
# Set up seurat objects
library(dplyr)
library(Seurat)

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)
Convert("fo8_thru_in8_GEX_for_CD8.h5ad", dest = "h5seurat", overwrite = TRUE)
covid <- LoadH5Seurat("fo8_thru_in8_GEX_for_CD8.h5seurat", assays = "counts")
Idents(covid) <- "louvain"
FeaturePlot(covid, features = "KIR3DL1")
VlnPlot(covid, features = "KIR3DL1")
DimPlot(covid, reduction = "umap", label = T)
covid19 <- CreateSeuratObject(counts = covid@assays$RNA@counts, meta.data = covid@meta.data, min.features = 200, min.cells = 3)
table(covid19@meta.data$ICU_status)

# create new column in the metadata for Neg/Pos status of gene
Idents(covid19) <- "ICU_status"
table(Idents(covid19))
covid19$disease <- plyr::mapvalues(
  x = covid19$ICU_status,
  from = c("Healthy", "ICU", "Non-ICU"),
  to = c("HD", "COVID-19", "COVID-19")
)
table(covid19@meta.data$disease)
covid19[["percent.mt"]] <- PercentageFeatureSet(covid19, pattern = "^MT-")
VlnPlot(covid19, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
covid19 <- subset(covid19, subset = nFeature_RNA > 800 & nFeature_RNA < 3000 & percent.mt < 10)

### subset CD8+ T cells
cd19_row_index <- grep("CD19", rownames(covid19$RNA@data))
trdc_row_index <- grep("TRDC", rownames(covid19$RNA@data))
cd8a_row_index <- grep("CD8A", rownames(covid19$RNA@data))
cd8b_row_index <- grep("CD8B", rownames(covid19$RNA@data))
length(which((covid19$RNA@data[rownames(covid19$RNA@data)[cd19_row_index], ] > 0)))
length(which((covid19$RNA@data[rownames(covid19$RNA@data)[trdc_row_index], ] > 0)))
length(which((covid19$RNA@data[rownames(covid19$RNA@data)[cd8b_row_index], ] > 0) | (covid19$RNA@data[rownames(covid19$RNA@data)[cd8a_row_index], ] > 0)))
length(which((covid19$RNA@data[rownames(covid19$RNA@data)[cd19_row_index], ] == 0) &
  (covid19$RNA@data[rownames(covid19$RNA@data)[trdc_row_index], ] == 0) &
  (covid19$RNA@data[rownames(covid19$RNA@data)[cd8a_row_index], ] > 0) |
  (covid19$RNA@data[rownames(covid19$RNA@data)[cd8b_row_index], ] > 0)))
# create new column in the metadata for Neg/Pos status of gene
covid19@meta.data$cd8 <- "Neg"
covid19@meta.data$cd8[which((covid19$RNA@data[rownames(covid19$RNA@data)[cd19_row_index], ] == 0) &
  (covid19$RNA@data[rownames(covid19$RNA@data)[trdc_row_index], ] == 0) &
  (covid19$RNA@data[rownames(covid19$RNA@data)[cd8a_row_index], ] > 0) |
  (covid19$RNA@data[rownames(covid19$RNA@data)[cd8b_row_index], ] > 0))] <- "Pos"
# subset only the positive cells; change the ident first.
Idents(covid19) <- "cd8"
table(Idents(covid19))
covid19 <- subset(covid19, ident = "Pos")
table(Idents(covid19))

# Load data of CD8 T from HC and MS
hcms <- readRDS("hc_ms_cd8.rds")

# Normalization
covid19 <- NormalizeData(covid19)
hcms <- NormalizeData(hcms)
# Find Variable Features
covid19 <- FindVariableFeatures(covid19, selection.method = "vst", nfeatures = 2000)
hcms <- FindVariableFeatures(hcms, selection.method = "vst", nfeatures = 2000)

CD8.anchor <- FindIntegrationAnchors(
  object.list = list(covid19, hcms),
  dims = 1:20
)
CD8.integrated <- IntegrateData(anchorset = CD8.anchor, dims = 1:20, verbose = T)

### first generate data and scale data in RNA assay
DefaultAssay(CD8.integrated) <- "RNA"
CD8.integrated <- NormalizeData(object = CD8.integrated, normalization.method = "LogNormalize", scale.factor = 1e4)
CD8.integrated <- FindVariableFeatures(object = CD8.integrated, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
CD8.integrated <- ScaleData(CD8.integrated, verbose = FALSE)

## change to integrated assay
DefaultAssay(CD8.integrated) <- "integrated"
dpi <- 300
png(file = "qc.png", width = dpi * 16, height = dpi * 8, units = "px", res = dpi, type = "cairo")
VlnPlot(object = CD8.integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

png(file = "umi-gene.png", width = dpi * 6, height = dpi * 5, units = "px", res = dpi, type = "cairo")
FeatureScatter(object = CD8.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

# Run the standard workflow for visualization and clustering
DefaultAssay(CD8.integrated) <- "integrated"
CD8.integrated <- ScaleData(CD8.integrated, verbose = FALSE)
CD8.integrated <- RunPCA(CD8.integrated, verbose = T, npcs = 50)
CD8.integrated <- ProjectDim(object = CD8.integrated)
ElbowPlot(object = CD8.integrated, ndims = 50)


### cluster
CD8.integrated <- FindNeighbors(object = CD8.integrated, dims = seq(1, 20))
CD8.integrated <- FindClusters(object = CD8.integrated, resolution = 0.3)

### tsne and umap
CD8.integrated <- RunUMAP(CD8.integrated, reduction = "pca", dims = 1:20)

tiff(filename = "UMAP.tiff", units = "cm", height = 20, width = 20, res = 300)
DimPlot(object = CD8.integrated, reduction = "umap", label = T) + NoLegend()
dev.off()

DefaultAssay(CD8.integrated) <- "RNA"
FeaturePlot(CD8.integrated, features = "IGLC2", label = T, sort.cell = T)

## remove contamination of B cells
CD8.integrated <- subset(CD8.integrated, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8))
DefaultAssay(CD8.integrated) <- "integrated"
CD8.integrated <- ScaleData(CD8.integrated, verbose = FALSE)
CD8.integrated <- RunPCA(CD8.integrated, verbose = T, npcs = 50)
CD8.integrated <- ProjectDim(object = CD8.integrated)
ElbowPlot(object = CD8.integrated, ndims = 50)
CD8.integrated <- FindNeighbors(object = CD8.integrated, dims = seq(1, 20))
CD8.integrated <- FindClusters(object = CD8.integrated, resolution = 0.3)
CD8.integrated <- RunUMAP(CD8.integrated, reduction = "pca", dims = 1:20)
tiff(filename = "UMAP.tiff", units = "cm", height = 20, width = 20, res = 300)
baseplot <- DimPlot(object = CD8.integrated, reduction = "umap", label = T, label.size = 4) + NoLegend()
baseplot + FontSize(x.title = 20, y.title = 20, x.text = 20, y.text = 20)
dev.off()
saveRDS(CD8.integrated, "CD8_integrated.rds")

tiff(filename = "KIR2DL3_disease.tiff", units = "cm", height = 20, width = 60, res = 300)
FeaturePlot(CD8.integrated, features = "KIR2DL3", order = T, split.by = "disease")
dev.off()


# identify conserved cell type markers
library(metap)
DefaultAssay(CD8.integrated) <- "RNA"
CD8.integrated <- NormalizeData(object = CD8.integrated, normalization.method = "LogNormalize", scale.factor = 1e4)
CD8.integrated <- FindVariableFeatures(object = CD8.integrated, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
CD8.integrated <- ScaleData(CD8.integrated, verbose = FALSE)

CD8.integrated <- RenameIdents(CD8.integrated,
  `0` = "naive CD8", `1` = "effector CD8", `2` = "GZMK+ CD8",
  `3` = "memory CD8", `4` = "KIR+ effector CD8", `5` = "proliferating CD8", `6` = "MAIT", `7` = "IFN-stim CD8", `8` = "effector CD8"
)
naive_CD8.markers <- FindConservedMarkers(CD8.integrated, ident.1 = "naive CD8", grouping.var = "disease", verbose = TRUE)
effector_CD8.markers <- FindConservedMarkers(CD8.integrated, ident.1 = "effector CD8", grouping.var = "disease", verbose = TRUE)
GZMK_CD8.markers <- FindConservedMarkers(CD8.integrated, ident.1 = "GZMK+ CD8", grouping.var = "disease", verbose = TRUE)
memory_CD8.markers <- FindConservedMarkers(CD8.integrated, ident.1 = "memory CD8", grouping.var = "disease", verbose = TRUE)
KIR_effector_CD8.markers <- FindConservedMarkers(CD8.integrated, ident.1 = "KIR+ effector CD8", grouping.var = "disease", verbose = TRUE)
proliferating_CD8.markers <- FindConservedMarkers(CD8.integrated, ident.1 = "proliferating CD8", grouping.var = "disease", verbose = TRUE)
MAIT.markers <- FindConservedMarkers(CD8.integrated, ident.1 = "MAIT", grouping.var = "disease", verbose = TRUE)
IFN_stim_CD8.markers <- FindConservedMarkers(CD8.integrated, ident.1 = "IFN-stim CD8", grouping.var = "disease", verbose = TRUE)

write.csv(IFN_stim_CD8.markers, "IFN_stim_CD8.csv")

## pick KIR+ CD8
kir3dl1_row_index <- grep("KIR3DL1", rownames(CD8.integrated$RNA@data))
kir3dl2_row_index <- grep("KIR3DL2", rownames(CD8.integrated$RNA@data))
kir2dl3_row_index <- grep("KIR2DL3", rownames(CD8.integrated$RNA@data))
kir2dl1_row_index <- grep("KIR2DL1", rownames(CD8.integrated$RNA@data))

# how many cells are positive?
length(which((CD8.integrated$RNA@data[rownames(CD8.integrated$RNA@data)[kir3dl1_row_index], ] > 0)))
length(which((CD8.integrated$RNA@data[rownames(CD8.integrated$RNA@data)[kir3dl2_row_index], ] > 0)))
length(which((CD8.integrated$RNA@data[rownames(CD8.integrated$RNA@data)[kir2dl3_row_index], ] > 0)))
length(which((CD8.integrated$RNA@data[rownames(CD8.integrated$RNA@data)[kir2dl1_row_index], ] > 0)))
length(which((CD8.integrated$RNA@data[rownames(CD8.integrated$RNA@data)[kir3dl1_row_index], ] > 0) | (CD8.integrated$RNA@data[rownames(CD8.integrated$RNA@data)[kir3dl2_row_index], ] > 0) | (CD8.integrated$RNA@data[rownames(CD8.integrated$RNA@data)[kir2dl3_row_index], ] > 0) | (CD8.integrated$RNA@data[rownames(CD8.integrated$RNA@data)[kir2dl1_row_index], ] > 0)))
# create new column in the metadata for Neg/Pos status of gene
CD8.integrated@meta.data$kir <- "KIR-"
CD8.integrated@meta.data$kir[which((CD8.integrated$RNA@data[rownames(CD8.integrated$RNA@data)[kir3dl1_row_index], ] > 0) | (CD8.integrated$RNA@data[rownames(CD8.integrated$RNA@data)[kir3dl2_row_index], ] > 0) | (CD8.integrated$RNA@data[rownames(CD8.integrated$RNA@data)[kir2dl3_row_index], ] > 0) | (CD8.integrated$RNA@data[rownames(CD8.integrated$RNA@data)[kir2dl1_row_index], ] > 0))] <- "KIR+"
# subset only the positive cells; change the ident first.
Idents(CD8.integrated) <- "kir"
table(Idents(CD8.integrated))
kirpos_cd8 <- subset(cd8.combined, ident = "Pos")
table(Idents(kirpos_cd8))
table(kirpos_cd8@meta.data$disease)
tiff(filename = "KIR_disease.tiff", units = "cm", height = 20, width = 55, res = 300)
baseplot <- DimPlot(CD8.integrated, group.by = "kir", split.by = "disease", cols = c("blue", "orange"), order = "KIR+")
baseplot + FontSize(x.title = 20, y.title = 20, x.text = 20, y.text = 20)
dev.off()

## DEG analysis between KIR+ effector and KIR- effector CD8
library(Seurat)
CD8 <- readRDS("CD8_integrated.rds")
effector_vs_KIRCD8 <- FindConservedMarkers(CD8, ident.1 = "KIR+ effector CD8", ident.2 = "effector CD8", grouping.var = "disease", verbose = TRUE)
write.csv(effector_vs_KIRCD8, "KIR+ vs KIR- effector.csv")
