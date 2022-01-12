library(dplyr)
library(Seurat)
library(patchwork)
library(DoubletFinder)
library(data.table)
library(magrittr)
library(ggplot2)
library(cowplot)
library(stringr)

########################################################################################################################################################################################################################

# CREATE OBJECT Lgr5Cre_SETDB1KO

# Load dataset A (Object A)
Lgr5Cre_Setdb1KO_A.data=Read10X(data.dir = "/media/dimbo/10T/data/talianidis_data/scRNA_Seq/Talianidis_GFP_2_fixed/analysis/Cellranger/run_count_21L004866/outs/filtered_feature_bc_matrix/")
Lgr5Cre_Setdb1KO_A = CreateSeuratObject(counts = Lgr5Cre_Setdb1KO_A.data, project = "Setdb1KO A", min.cells = 3, min.features = 200)
# Load dataset B (Object B)
Lgr5Cre_Setdb1KO_B.data <- Read10X(data.dir = "/media/dimbo/10T/data/talianidis_data/scRNA_Seq/Talianidis_GFP_2_fixed/analysis/Cellranger/run_count_21L004870/outs/filtered_feature_bc_matrix/")
Lgr5Cre_Setdb1KO_B = CreateSeuratObject(counts = Lgr5Cre_Setdb1KO_B.data, project = "Setdb1KO B", min.cells = 3, min.features = 200)

# Merge objects
Lgr5Cre_Setdb1KO=merge(Lgr5Cre_Setdb1KO_A, y=Lgr5Cre_Setdb1KO_B, add.cell.ids=c("Rep_A","Rep_B"), project="Lgr5Cre_Setdb1KO")

##################################################################################################################################################################################################################

# CREATE OBJECT LGR5Cre_WT

# Load dataset A (Object A)
Lgr5Cre_WT_A.data=Read10X(data.dir = "/media/dimbo/10T/data/talianidis_data/scRNA_Seq/Talianidis_GFP_2_fixed/analysis/Cellranger/run_count_21L004858/outs/filtered_feature_bc_matrix/")
Lgr5Cre_WT_A = CreateSeuratObject(counts = Lgr5Cre_WT_A.data, project = "Lgr5Cre A", min.cells = 3, min.features = 200)
# Load dataset B (Object B)
Lgr5Cre_WT_B.data <- Read10X(data.dir = "/media/dimbo/10T/data/talianidis_data/scRNA_Seq/Talianidis_GFP_2_fixed/analysis/Cellranger/run_count_21L004862/outs/filtered_feature_bc_matrix/")
Lgr5Cre_WT_B = CreateSeuratObject(counts = Lgr5Cre_WT_B.data, project = "Lgr5Cre B", min.cells = 3, min.features = 200)

# Merge objects
Lgr5Cre_WT=merge(Lgr5Cre_WT_A, y=Lgr5Cre_WT_B, add.cell.ids=c("Rep_A","Rep_B"), project="Lgr5Cre_WT")

##################################################################################################################################################################################################################

#FIND DOUBLETS IN EACH CONDITION SEPARATELY

Lgr5Cre_Setdb1KO <- NormalizeData(Lgr5Cre_Setdb1KO, normalization.method = "LogNormalize", scale.factor = 10000)

Lgr5Cre_Setdb1KO <- FindVariableFeatures(Lgr5Cre_Setdb1KO, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(Lgr5Cre_Setdb1KO)
Lgr5Cre_Setdb1KO <- ScaleData(Lgr5Cre_Setdb1KO, features = all.genes)

Lgr5Cre_Setdb1KO <- RunPCA(Lgr5Cre_Setdb1KO, features = VariableFeatures(object = Lgr5Cre_Setdb1KO))

Lgr5Cre_Setdb1KO <- RunUMAP(Lgr5Cre_Setdb1KO, dims = 1:10)

nExp <- round(ncol(Lgr5Cre_Setdb1KO) * 0.069)
Lgr5Cre_Setdb1KO.filt <- doubletFinder_v3(Lgr5Cre_Setdb1KO, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

DF.name = colnames(Lgr5Cre_Setdb1KO.filt@meta.data)[grepl("DF.classification", colnames(Lgr5Cre_Setdb1KO.filt@meta.data))]

pdf("CowPlot_Setdb1KO.pdf", width=15, height=8)
cowplot::plot_grid(ncol = 2, DimPlot(Lgr5Cre_Setdb1KO.filt, group.by = DF.name) + NoAxes())
dev.off()

pdf("VlnPlot_Setdb1KO.pdf", width=15, height=8)
VlnPlot(Lgr5Cre_Setdb1KO.filt, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
dev.off()

Lgr5Cre_Setdb1KO.filt = Lgr5Cre_Setdb1KO.filt[, Lgr5Cre_Setdb1KO.filt@meta.data[, DF.name] == "Singlet"]

##################################################################################################################################################################################################################

Lgr5Cre_WT <- NormalizeData(Lgr5Cre_WT, normalization.method = "LogNormalize", scale.factor = 10000)

Lgr5Cre_WT <- FindVariableFeatures(Lgr5Cre_WT, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(Lgr5Cre_WT)
Lgr5Cre_WT <- ScaleData(Lgr5Cre_WT, features = all.genes)

Lgr5Cre_WT <- RunPCA(Lgr5Cre_WT, features = VariableFeatures(object = Lgr5Cre_Setdb1KO))

Lgr5Cre_WT <- RunUMAP(Lgr5Cre_WT, dims = 1:10)

nExp <- round(ncol(Lgr5Cre_WT) * 0.057)
Lgr5Cre_WT.filt <- doubletFinder_v3(Lgr5Cre_WT, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

DF.name = colnames(Lgr5Cre_WT.filt@meta.data)[grepl("DF.classification", colnames(Lgr5Cre_WT.filt@meta.data))]

pdf("CowPlot_WT.pdf", width=15, height=8)
cowplot::plot_grid(ncol = 2, DimPlot(Lgr5Cre_WT.filt, group.by = DF.name) + NoAxes())
dev.off()

pdf("VlnPlot_WT.pdf", width=15, height=8)
VlnPlot(Lgr5Cre_WT.filt, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
dev.off()

##################################################################################################################################################################################################################

# Merge Setdb1KO with WT to new object Lgr5Cre_MERGED
Lgr5Cre_MERGED=merge(Lgr5Cre_Setdb1KO.filt, y=Lgr5Cre_WT.filt, add.cell.ids=c("Setdb1KO","Lgr5Cre"), project="Lgr5Cre_MERGED")

# Change replicates names within object. e.g(Setdb1KO_repA & Setdb1KO_repB -> Setdb1KO)
Lgr5Cre_MERGED$orig.ident=plyr::mapvalues(x=Lgr5Cre_MERGED$orig.ident, from = c("Setdb1KO A", "Setdb1KO B", "Lgr5Cre A", "Lgr5Cre B"), to = c("Setdb1KO", "Setdb1KO", "Lgr5Cre", "Lgr5Cre"))

##################################################################################################################################################################################################################

# QC and selecting cells
Lgr5Cre_MERGED[["percent.mt"]] <- PercentageFeatureSet(Lgr5Cre_MERGED, pattern = "^mt-")

# Violin Plot
pdf("violin_plot_Lgr5Cre.pdf", width=18, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Violin Plot with merged replicates
pdf("violin_plot_Lgr5Cre_merged.pdf", width=18, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident", ncol = 3)
dev.off()

# Feature Scatter Plot
plot1 <- FeatureScatter(Lgr5Cre_MERGED, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Lgr5Cre_MERGED, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("feauture_scatter.pdf", width=15, height=8)
plot1 + plot2
dev.off()

# Feature Scatter Plot
plot1 <- FeatureScatter(Lgr5Cre_MERGED, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident")
plot2 <- FeatureScatter(Lgr5Cre_MERGED, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")
pdf("feauture_scatter_merged.pdf", width=15, height=8)
plot1 + plot2
dev.off()

# Filtering cells 
Lgr5Cre_MERGED <- subset(Lgr5Cre_MERGED, subset = nFeature_RNA > 200 & nFeature_RNA < 3800 & percent.mt < 5)

# Normalize the data
Lgr5Cre_MERGED <- NormalizeData(Lgr5Cre_MERGED, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features
Lgr5Cre_MERGED <- FindVariableFeatures(Lgr5Cre_MERGED, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Lgr5Cre_MERGED), 10)

# Plot variable features
pdf("variable_features_Lgr5Cre.pdf", width=15, height=8)
plot1 <- VariableFeaturePlot(Lgr5Cre_MERGED)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()

# Scaling the data
all.genes <- rownames(Lgr5Cre_MERGED)
Lgr5Cre_MERGED <- ScaleData(Lgr5Cre_MERGED, features = all.genes)

# Performing linear dimensional reduction
Lgr5Cre_MERGED <- RunPCA(Lgr5Cre_MERGED, features = VariableFeatures(object = Lgr5Cre_MERGED))
print(Lgr5Cre_MERGED[["pca"]], dims = 1:5, nfeatures = 5)

# VizDim Plot
pdf("Viz_DimPlot_PCA_Lgr5Cre.pdf", width=15, height=8)
VizDimLoadings(Lgr5Cre_MERGED, dims = 1:2, reduction = "pca")
dev.off()

# DimPlot
pdf("DimPlot_PCA_Lgr5Cre.pdf", width=15, height=8)
VizDimLoadings(Lgr5Cre_MERGED, dims = 1:2, reduction = "pca")
dev.off()

# Dim Heatmap
pdf("DimHeatmap_Lgr5Cre.pdf", width=15, height=8)
DimHeatmap(Lgr5Cre_MERGED, dims = 1, cells = 500, balanced = TRUE)
dev.off()

# Determine the dimensionality of the dataset
Lgr5Cre_MERGED <- JackStraw(Lgr5Cre_MERGED, num.replicate = 100)
Lgr5Cre_MERGED <- ScoreJackStraw(Lgr5Cre_MERGED, dims = 1:20)

pdf("JackstrawPlot_Lgr5Cre.pdf", width=15, height=8)
JackStrawPlot(Lgr5Cre_MERGED, dims = 1:15)
dev.off()

pdf("ElbowPlot_Lgr5Cre.pdf", width=15, height=8)
ElbowPlot(Lgr5Cre_MERGED)
dev.off()

# Cluster the Cells
Lgr5Cre_MERGED <- FindNeighbors(Lgr5Cre_MERGED, dims = 1:10)
Lgr5Cre_MERGED <- FindClusters(Lgr5Cre_MERGED, resolution = 0.5)

# Run non linear dimensional reduction
reticulate::py_install(packages = 'umap-learn')
Lgr5Cre_MERGED <- RunUMAP(Lgr5Cre_MERGED, dims = 1:10)

# Umap plot
pdf("umapPlot_Lgr5Cre.pdf", width=15, height=8)
DimPlot(Lgr5Cre_MERGED, reduction = "umap")
dev.off()

# Umap plot merged
pdf("umapPlot_Lgr5Cre_merged.pdf", width=15, height=8)
DimPlot(Lgr5Cre_MERGED, reduction = "umap", group.by = "orig.ident")
dev.off()

# Split objects 
SplitObject(Lgr5Cre_MERGED, split.by = "ident")
n_cells=(FetchData(Lgr5Cre_MERGED, var=c("ident", "orig.ident")) %>% dplyr::count(ident, orig.ident) %>% tidyr::spread(ident, n))

# Umap plot merged replicates per condition
pdf("umapPlot_Lgr5Cre_conditions.pdf", width=15, height=8)
DimPlot(Lgr5Cre_MERGED,label=TRUE, split.by="orig.ident") + NoLegend()
dev.off()

# Feauture plot top 10 genes
pdf("feature_plot_top_10_MostVariantGenes.pdf", width=15, height=8)
FeaturePlot(Lgr5Cre_MERGED, features = top10)
dev.off()

# Dot plot per top 10 genes
pdf("DotPlot_Lgr5cre.pdf", width=15, height=8)
DotPlot(Lgr5Cre_MERGED, cols = c("blue", "orange"), group.by = "orig.ident", features = top10)
dev.off()

# Dot plot per selected Markers
pdf("Top10_Markers.pdf", width=15, height=8)
DotPlot(Lgr5Cre_MERGED, features = c("Reg3g", "Gsdmc4", "Prss32", "Krt8", "Elf3", "Sis", "Fabp1", "Hnf4a", "Tmigd1", "Fabp6"), split.by = "orig.ident")
dev.off()

##################################################################################################################################################################################################################

# Extract metadata
md = Lgr5Cre_MERGED@meta.data %>% as.data.table()
md[, .N, by = c("orig.ident", "seurat_clusters")]
md[, .N, by = c("orig.ident", "seurat_clusters")] %>% dcast(., orig.ident ~ seurat_clusters, value.var = "N")

# Extract and save all genes
write.table(all.genes, file = "all_genes.tsv", sep = "\t", row.names = TRUE)

# Find all markers distinguishing cluster 0 from the rest of the clusters
cluster0.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 0, ident.2 = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14), min.pct = 0.25, logfc.threshold = 0.37)
cluster1.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 1, ident.2 = c(0,2,3,4,5,6,7,8,9,10,11,12,13,14), min.pct = 0.25, logfc.threshold = 0.37)
cluster2.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 2, ident.2 = c(0,1,3,4,5,6,7,8,9,10,11,12,13,14), min.pct = 0.25, logfc.threshold = 0.37)
cluster3.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 3, ident.2 = c(0,1,2,4,5,6,7,8,9,10,11,12,13,14), min.pct = 0.25, logfc.threshold = 0.37)
cluster4.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 4, ident.2 = c(0,1,2,3,5,6,7,8,9,10,11,12,13,14), min.pct = 0.25, logfc.threshold = 0.37)
cluster5.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 5, ident.2 = c(0,1,2,3,4,6,7,8,9,10,11,12,13,14), min.pct = 0.25, logfc.threshold = 0.37)
cluster6.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 6, ident.2 = c(0,1,2,3,4,5,7,8,9,10,11,12,13,14), min.pct = 0.25, logfc.threshold = 0.37)
cluster7.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 7, ident.2 = c(0,1,2,3,4,5,6,8,9,10,11,12,13,14), min.pct = 0.25, logfc.threshold = 0.37)
cluster8.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 8, ident.2 = c(0,1,2,3,4,5,6,7,9,10,11,12,13,14), min.pct = 0.25, logfc.threshold = 0.37)
cluster9.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 9, ident.2 = c(0,1,2,3,4,5,6,7,8,10,11,12,13,14), min.pct = 0.25, logfc.threshold = 0.37)
cluster10.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 10, ident.2 = c(0,1,2,3,4,5,6,7,8,9,11,12,13,14), min.pct = 0.25, logfc.threshold = 0.37)
cluster11.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 11, ident.2 = c(0,1,2,3,4,5,6,7,8,9,10,12,13,14), min.pct = 0.25, logfc.threshold = 0.37)
cluster12.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 12, ident.2 = c(0,1,2,3,4,5,6,7,8,9,10,11,13,14), min.pct = 0.25, logfc.threshold = 0.37)
cluster13.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 13, ident.2 = c(0,1,2,3,4,5,6,7,8,9,10,11,12,14), min.pct = 0.25, logfc.threshold = 0.37)
cluster14.markers <- FindMarkers(Lgr5Cre_MERGED,ident.1 = 14, ident.2 = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13), min.pct = 0.25, logfc.threshold = 0.37)

# Save files
write.table(cluster0.markers, file = "cluster0.markers.tsv", sep = "\t", row.names = TRUE)
write.table(cluster1.markers, file = "cluster1.markers.tsv", sep = "\t", row.names = TRUE)
write.table(cluster2.markers, file = "cluster2.markers.tsv", sep = "\t", row.names = TRUE)
write.table(cluster3.markers, file = "cluster3.markers.tsv", sep = "\t", row.names = TRUE)
write.table(cluster4.markers, file = "cluster4.markers.tsv", sep = "\t", row.names = TRUE)
write.table(cluster5.markers, file = "cluster5.markers.tsv", sep = "\t", row.names = TRUE)
write.table(cluster6.markers, file = "cluster6.markers.tsv", sep = "\t", row.names = TRUE)
write.table(cluster7.markers, file = "cluster7.markers.tsv", sep = "\t", row.names = TRUE)
write.table(cluster8.markers, file = "cluster8.markers.tsv", sep = "\t", row.names = TRUE)
write.table(cluster9.markers, file = "cluster9.markers.tsv", sep = "\t", row.names = TRUE)
write.table(cluster10.markers, file = "cluster10.markers.tsv", sep = "\t", row.names = TRUE)
write.table(cluster11.markers, file = "cluster11.markers.tsv", sep = "\t", row.names = TRUE)
write.table(cluster12.markers, file = "cluster12.markers.tsv", sep = "\t", row.names = TRUE)
write.table(cluster13.markers, file = "cluster13.markers.tsv", sep = "\t", row.names = TRUE)
write.table(cluster14.markers, file = "cluster14.markers.tsv", sep = "\t", row.names = TRUE)

# Mat code
#Lgr5Cre_MERGED.markers = FindAllMarkers(Lgr5Cre_MERGED, only.pos = F, min.pct = 0.25, logfc.threshold = 0.37)

##################################################################################################################################################################################################################

# Feature Plots based on most variable features per cluster

pdf("Feature plot cluster_0 Markers 1-5.pdf", width=20, height=30)
FeaturePlot(Lgr5Cre_MERGED, features = c("Olfm4", "Gkn3", "Ifitm3", "Jaml", "2410006H16Rik"), split.by = "orig.ident", cols = c("grey" ,"blue"))
dev.off()


##################################################################################################################################################################################################################

# DotPlots per Cluster

pdf("Enterocyte_Mature Markers.pdf", width=15, height=8)
p = DotPlot(Lgr5Cre_MERGED,features = c("Elf3","Sis","Fabp1","Hnf4aos","Hnf4a", "Hnf4g", "Tmigd1","Fabp6","Slc51b","Slc51a","Mep1a","Fam151a","Naaladl1","Slc34a2","Plb1","Nudt4","Dpep1","Pmp22","Xpnpep2","Muc3","Neu1","Clec2h","Phgr1","Prss30","Aldob","Alpi","Apoa1","Apoa4","Lct"),cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Enterocyte Mature Markers")
dev.off()

##################################################################################################################################################################################################################

# Dot plots per condition

pdf("Enterocyte Immature Markers Grouped.pdf", width=15, height=8)
p =  DotPlot(Lgr5Cre_MERGED,features = c("Reg3g", "Gsdmc4", "Prss32", "Krt8"), group.by = "orig.ident" ,cols = c("blue","orange"))
p + ggtitle("Enterocyte Immature Markers") + RotatedAxis()
dev.off()

##################################################################################################################################################################################################################

# Umap renamed cluster
pdf("Umap_renamed_Clusterd Markers.pdf", width=15, height=8)
Lgr5Cre_MERGED_Renamed=RenameIdents(Lgr5Cre_MERGED,  `0` = "Stem", `1` = "Enterocyte Progenitor I", `2` = "Goblet I", `3` = "Stem - TA", `4` = "?", `5` = "Enterocyte Progenitor II", `6` = "Goblet II", `7` = "Enterocyte Mature I", `8` = "Paneth", `9` = "Enterocyte Mature II", `10` = "Enderoendocrine I", `11` = "Tuft", `12` = "Enterocyte immature", `13` = "Enteroendocrine II", `14` = "Necroptosis")
DimPlot(Lgr5Cre_MERGED_Renamed, label = TRUE)
dev.off()

##################################################################################################################################################################################################################

# Split violin plots per markers
pdf("Violin Plot_Birc5.pdf", width=15, height=10)
plots <- VlnPlot(Lgr5Cre_MERGED, features = c(""), split.by = "orig.ident", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()

# Selected Markers

pdf("Violin Plot_Sis.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("Sis"), split.by = "orig.ident")
dev.off()

# For cluster 0

pdf("Violin Plot_Olfm4.pdf", width=15, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("Olfm4"), split.by = "orig.ident")
dev.off()

##################################################################################################################################################################################################################

# Split Ridge plots per markers

pdf("RidgePlot_Olfm4.pdf", width=15, height=10)
RidgePlot(Lgr5Cre_MERGED_Renamed, features = c("Olfm4"), ncol = 2)
dev.off()

##################################################################################################################################################################################################################

# Cell cycle

# Read in the expression matrix The first row is a header row, the first column is rownames
exp.mat <- read.table(file = "/media/dimbo/10T/data/talianidis_data/scRNA_Seq/Talianidis_GFP_2_fixed/analysis/Seurat/cell_cycle/", header = TRUE, as.is = TRUE, row.names = 1)

s.genes <- cc.genes$s.genes
s.genes <- tolower(s.genes)
s.genes = str_to_title(s.genes)

g2m.genes <- cc.genes$g2m.genes
g2m.genes <- tolower(g2m.genes)
g2m.genes = str_to_title(g2m.genes)

# Assign Cell-Cycle scores
Lgr5Cre_MERGED <- CellCycleScoring(Lgr5Cre_MERGED, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

pdf("PCA on Cell Cycle Genes.pdf", width=15, height=8)
Lgr5Cre_MERGED <- RunPCA(Lgr5Cre_MERGED, features = c(s.genes, g2m.genes))
DimPlot(Lgr5Cre_MERGED)
dev.off()

# Regress out cell cycle scores during data scaling
Lgr5Cre_MERGED <- ScaleData(Lgr5Cre_MERGED, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Lgr5Cre_MERGED))

Lgr5Cre_MERGED <- RunPCA(Lgr5Cre_MERGED, features = VariableFeatures(Lgr5Cre_MERGED), nfeatures.print = 10)

Lgr5Cre_MERGED <- RunPCA(Lgr5Cre_MERGED, features = c(s.genes, g2m.genes))
DimPlot(Lgr5Cre_MERGED)

##################################################################################################################################################################################################################

# Create dot plots - on specific clusters and selected markers

pdf("Setdb1KO_Enterocyte Mature top_markers.pdf", width=15, height=8)
p =  DotPlot(Setdb1KO,features = c("Smim24","St3gal4","2200002D01Rik","Krt19","Dmbt1","Gna11","Prap1","Gpx1","Serpinb6a","Reg1","Rbp2","Fabp1","Apoa1","Gsta1","Fabp2","Arg2","Gstm3","Adh6a","Sis","Apoc3","S100g"), idents = c(0,3),cols = c("blue","orange"))
p + ggtitle("Setdb1KO") + RotatedAxis()
dev.off()


##################################################################################################################################################################################################################

# Create metadata from main object. e.g export WT and Setdb1KO from Merged object and create two new objects

split_obj = SplitObject(Lgr5Cre_MERGED, split.by = "orig.ident")

# See how the objects are named
View(split_obj)

Setdb1KO = split_obj[["Setdb1KO"]]

WT = split_obj[["Lgr5Cre"]]

##################################################################################################################################################################################################################

# Calculate averaged expression values for each identity class
avg_exp_Setdb1KO = (AverageExpression(object = Setdb1KO))
avg_exp_WT = (AverageExpression(object = WT))

write.table(avg_exp_Setdb1KO, file = "avg_exp_Setdb1KO.txt", sep = "\t",
            row.names = TRUE, col.names = NA)
write.table(avg_exp_WT, file = "avg_exp_WT.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

##################################################################################################################################################################################################################

# Create violin plots for setdb1KO and WT metadata (merged normalized-processed and then split) for nfeatures 
# ncount per cluster.

pdf("Setdb1KO-ViolinPlot_per_Cluster.pdf", width=18, height=10)
VlnPlot(Setdb1KO, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
dev.off()

pdf("WT-ViolinPlot_per_Cluster.pdf", width=18, height=10)
VlnPlot(WT, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
dev.off()

pdf("Merged_ViolinPlot_per_Cluster.pdf", width=18, height=10)
VlnPlot(Lgr5Cre_MERGED, features = c("nFeature_RNA", "nCount_RNA"), split.by = "orig.ident", ncol = 2)
dev.off()

# Keep cells with only 2000 reads
#Lgr5Cre_MERGED <- subset(Lgr5Cre_MERGED, subset = nCount_RNA > 2000)
