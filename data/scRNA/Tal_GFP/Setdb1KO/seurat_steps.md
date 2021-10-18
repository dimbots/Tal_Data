## Load Libraries
library(dplyr)
library(Seurat)
library(patchwork)

## Seurat Commands

## Creat object A
```
Lgr5Cre_Setdb1KO_A.data=Read10X(data.dir="/media/dimbo/10T/data/talianidis_data/single_cell/Talianidis_GFP_2_fixed/Merge_Samples/seurat/21L004866/filtered_feature_bc_matrix/")
Lgr5Cre_Setdb1KO_A=CreateSeuratObject(counts=Lgr5Cre_Setdb1KO_A.data, project="Lgr5Cre_Setdb1KO_A")
```
## Creat object B
```
Lgr5Cre_Setdb1KO_B.data=Read10X(data.dir="/media/dimbo/10T/data/talianidis_data/single_cell/Talianidis_GFP_2_fixed/Merge_Samples/seurat/21L004870/filtered_feature_bc_matrix/")
Lgr5Cre_Setdb1KO_B=CreateSeuratObject(counts=Lgr5Cre_Setdb1KO_B.data, project="Lgr5Cre_Setdb1KO_B")
```
## Merge objects
```
Lgr5Cre_Setdb1KO=merge(Lgr5Cre_Setdb1KO_A, y=Lgr5Cre_Setdb1KO_B, add.cell.ids=c("Replicate_A","Replicate_B"), project="Lgr5Cre_Setdb1_KO")
```
## QC and selecting cells for further analysis
```
Lgr5Cre_Setdb1KO[["percent.mt"]] <- PercentageFeatureSet(Lgr5Cre_Setdb1KO, pattern = "^mt-")
```
## Visualize QC metrics as a violin plot
```
pdf("violin_plot_Setdb1KO", width=15, height=8)
VlnPlot(Lgr5Cre_Setdb1KO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pdf("violin_plot_Setdb1KO", width=15, height=8)
```
## Visualize feauter scatter plots
```
plot1 <- FeatureScatter(Lgr5Cre_Setdb1KO, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Lgr5Cre_Setdb1KO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("feauture_scatter_Setdb1KO", width=15, height=8)
plot1 + plot2
dev.off()
```
## Normalize the data
```
Lgr5Cre_Setdb1KO = NormalizeData(Lgr5Cre_Setdb1KO, normalization.method = "LogNormalize", scale.factor = 10000)
```
## Find Variable Features
```
Lgr5Cre_Setdb1KO = FindVariableFeatures(Lgr5Cre_Setdb1KO, selection.method = "vst", nfeatures = 2000)
```
## Identify the 10 most highly variable genes
```
top10 = head(VariableFeatures(Lgr5Cre_Setdb1KO), 10)
```
## plot variable features with and without labels
```
pdf("variable_features_Setdb1KO", width=15, height=8)
plot1 <- VariableFeaturePlot(Lgr5Cre_Setdb1KO)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
dev.off()
```
## Scaling the data
```
all.genes <- rownames(Lgr5Cre_Setdb1KO)
Lgr5Cre_Setdb1KO <- ScaleData(Lgr5Cre_Setdb1KO, features = all.genes)
```
## Perform linear dimensional reduction
```
Lgr5Cre_Setdb1KO <- RunPCA(Lgr5Cre_Setdb1KO, features = VariableFeatures(object = Lgr5Cre_Setdb1KO))
```
## Visualize PCA results
```
print(Lgr5Cre_Setdb1KO[["pca"]], dims = 1:5, nfeatures = 5)
pdf("Viz_DimPlot_PCA", width=15, height=8)
VizDimLoadings(Lgr5Cre_Setdb1KO, dims = 1:2, reduction = "pca")
dev.off()
```
##
```
pdf("DimPlot_PCA", width=15, height=8)
DimPlot(Lgr5Cre_Setdb1KO, reduction = "pca")
dev.off()
```
## DimHeatmap plot
```
pdf("DimHeatmap", width=15, height=8)
DimHeatmap(Lgr5Cre_Setdb1KO, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()
```
## Determine the dimensionality of the dataset
```
Lgr5Cre_Setdb1KO <- JackStraw(Lgr5Cre_Setdb1KO, num.replicate = 100)
Lgr5Cre_Setdb1KO <- ScoreJackStraw(Lgr5Cre_Setdb1KO, dims = 1:20)
```
## JackstrawPlot
```
pdf("JackStrawPlot", width=15, height=8)
JackStrawPlot(Lgr5Cre_Setdb1KO, dims = 1:15)
dev.off()
```
## ElbowPlot
```
pdf("elbow_plot", width=15, height=8)
ElbowPlot(Lgr5Cre_Setdb1KO)
dev.off()
```
## Cluster the cells
```
Lgr5Cre_Setdb1KO <- FindNeighbors(Lgr5Cre_Setdb1KO, dims = 1:10)
Lgr5Cre_Setdb1KO <- FindClusters(Lgr5Cre_Setdb1KO, resolution = 0.5)
```
## Run non-linear dimensional reduction
```
pdf("umapPlot", width=15, height=8)
Lgr5Cre_Setdb1KO <- RunUMAP(Lgr5Cre_Setdb1KO, dims = 1:10)
DimPlot(Lgr5Cre_Setdb1KO, reduction = "umap")
dev.off()
```
## SplitObjects
```
SplitObject(Lgr5Cre_Setdb1KO, split.by = "ident")
n_cells=(FetchData(Lgr5Cre_Setdb1KO, var=c("ident", "orig.ident")) %>% dplyr::count(ident, orig.ident) %>% tidyr::spread(ident, n))
pdf("umapPlot_replicates", width=15, height=8)
dev.off()
```