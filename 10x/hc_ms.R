setwd("/data/10x")
# Set up seurat objects
library(dplyr)
library(Seurat)

ReadData <- function(dataset_name, data_dir, disease_label) {
  data <- Read10X(data.dir = data_dir)
  seurat_obj <- CreateSeuratObject(counts = data, project = dataset_name, min.cells = 3, min.features = 200)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 800 & nFeature_RNA < 3000 & percent.mt < 10)
  seurat_obj@meta.data$disease <- disease_label

  return(seurat_obj)
}

MS9020 <- ReadData("MS9020", "out_cellranger_count/9020/filtered_feature_bc_matrix", "MS")
MS658 <- ReadData("MS658", "out_cellranger_count/658/filtered_feature_bc_matrix", "MS")
MS0109 <- ReadData("MS0109", "out_cellranger_count/01092019/filtered_feature_bc_matrix", "MS")
MS0816 <- ReadData("MS0816", "out_cellranger_count/08162018/filtered_feature_bc_matrix", "MS")
MS0910 <- ReadData("MS0910", "out_cellranger_count/09102018/filtered_feature_bc_matrix", "MS")
MS1113 <- ReadData("MS1113", "out_cellranger_count/11132018/filtered_feature_bc_matrix", "MS")
HC11 <- ReadData("HC11", "out_cellranger_count/HC11/filtered_feature_bc_matrix", "HD")
HC12 <- ReadData("HC12", "out_cellranger_count/HC12/filtered_feature_bc_matrix", "HD")
HC13 <- ReadData("HC13", "out_cellranger_count/HC13/filtered_feature_bc_matrix", "HD")
HC14 <- ReadData("HC14", "out_cellranger_count/HC14/filtered_feature_bc_matrix", "HD")

SubsetTCells <- function(dataset) {
  ## subset CD8+ T cells
  cd8a_row_index <- grep("CD8A", rownames(dataset$RNA@data))[1]
  cd8b_row_index <- grep("CD8B", rownames(dataset$RNA@data))
  trdc_row_index <- grep("TRDC", rownames(dataset$RNA@data))

  length(which((dataset$RNA@data[rownames(dataset$RNA@data)[cd8a_row_index], ] > 0)))
  length(which((dataset$RNA@data[rownames(dataset$RNA@data)[cd8b_row_index], ] > 0)))
  length(which((dataset$RNA@data[rownames(dataset$RNA@data)[cd8a_row_index], ] > 0) & (dataset$RNA@data[rownames(dataset$RNA@data)[cd8b_row_index], ] > 0) & (dataset$RNA@data[rownames(dataset$RNA@data)[trdc_row_index], ] == 0)))

  # create new column in the metadata for Neg/Pos status of gene
  dataset@meta.data$cd8 <- "Neg"
  dataset@meta.data$cd8[which((dataset$RNA@data[rownames(dataset$RNA@data)[cd8a_row_index], ] > 0) & (dataset$RNA@data[rownames(dataset$RNA@data)[cd8b_row_index], ] > 0) & (dataset$RNA@data[rownames(dataset$RNA@data)[trdc_row_index], ] == 0))] <- "Pos"

  # subset only the positive cells; change the ident first.
  Idents(dataset) <- "cd8"
  table(Idents(dataset))

  dataset_cd8 <- subset(dataset, ident = "Pos")
  table(Idents(dataset_cd8))

  return(dataset_cd8)
}

## subset CD8+ T cells
MS0109_cd8 <- SubsetTCells(MS0109)
MS658_cd8 <- SubsetTCells(MS658)
MS0816_cd8 <- SubsetTCells(MS0816)
MS0910_cd8 <- SubsetTCells(MS0910)
MS9020_cd8 <- SubsetTCells(MS9020)
MS1113_cd8 <- SubsetTCells(MS1113)
HC11_cd8 <- SubsetTCells(HC11)
HC12_cd8 <- SubsetTCells(HC12)
HC13_cd8 <- SubsetTCells(HC13)
HC14_cd8 <- SubsetTCells(HC14)

# Normalization
MS0109_cd8 <- NormalizeData(MS0109_cd8)
MS658_cd8 <- NormalizeData(MS658_cd8)
MS0816_cd8 <- NormalizeData(MS0816_cd8)
MS0910_cd8 <- NormalizeData(MS0910_cd8)
MS9020_cd8 <- NormalizeData(MS9020_cd8)
MS1113_cd8 <- NormalizeData(MS1113_cd8)
HC11_cd8 <- NormalizeData(HC11_cd8)
HC12_cd8 <- NormalizeData(HC12_cd8)
HC13_cd8 <- NormalizeData(HC13_cd8)
HC14_cd8 <- NormalizeData(HC14_cd8)

# Find Variable Features
MS0109_cd8 <- FindVariableFeatures(MS0109_cd8, selection.method = "vst", nfeatures = 2000)
MS658_cd8 <- FindVariableFeatures(MS658_cd8, selection.method = "vst", nfeatures = 2000)
MS0816_cd8 <- FindVariableFeatures(MS0816_cd8, selection.method = "vst", nfeatures = 2000)
MS0910_cd8 <- FindVariableFeatures(MS0910_cd8, selection.method = "vst", nfeatures = 2000)
MS9020_cd8 <- FindVariableFeatures(MS9020_cd8, selection.method = "vst", nfeatures = 2000)
MS1113_cd8 <- FindVariableFeatures(MS1113_cd8, selection.method = "vst", nfeatures = 2000)
HC11_cd8 <- FindVariableFeatures(HC11_cd8, selection.method = "vst", nfeatures = 2000)
HC12_cd8 <- FindVariableFeatures(HC12_cd8, selection.method = "vst", nfeatures = 2000)
HC13_cd8 <- FindVariableFeatures(HC13_cd8, selection.method = "vst", nfeatures = 2000)
HC14_cd8 <- FindVariableFeatures(HC14_cd8, selection.method = "vst", nfeatures = 2000)

# Perform integration
cd8.anchors <- FindIntegrationAnchors(object.list = list(MS0109_cd8, MS658_cd8, MS0816_cd8, MS0910_cd8, MS9020_cd8, MS1113_cd8, HC11_cd8, HC12_cd8, HC13_cd8, HC14_cd8), dims = 1:20)
cd8.combined <- IntegrateData(anchorset = cd8.anchors, dims = 1:20, verbose = T)

saveRDS(cd8.combined, file = "hc_ms_cd8.rds")
