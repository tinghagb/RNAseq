#clean environment
ls()
remove(list = ls())
#Can also clear global envi Ctrl+Shift+F10

#loadling Seurat and other required packages
library(dplyr)
library(Seurat)
library(patchwork)
#Inspect the package version
packageVersion('Seurat')

# Load the tissue ILC dataset
ILC.data_wt1 <- Read10X(data.dir = "~/Workspace/GSE102299_RAW/wt_PBS_1")
ILC.data_wt2 <- Read10X(data.dir = "~/Workspace/GSE102299_RAW/wt_PBS_2")
ILC.data_hdm1 <- Read10X(data.dir = "~/Workspace/GSE102299_RAW/wt_HDM_1")
ILC.data_hdm2 <- Read10X(data.dir = "~/Workspace/GSE102299_RAW/wt_HDM_2")


# Initialize the Seurat object with the raw (non-normalized data).
ILC_wt1 <- CreateSeuratObject(counts = ILC.data_wt1, project = "PBS1", min.cells = 3, min.features = 200)
ILC_wt2 <- CreateSeuratObject(counts = ILC.data_wt2, project = "PBS2", min.cells = 3, min.features = 200)
ILC_hdm1 <- CreateSeuratObject(counts = ILC.data_hdm1, project = "HDM1", min.cells = 3, min.features = 200)
ILC_hdm2 <- CreateSeuratObject(counts = ILC.data_hdm2, project = "HDM2", min.cells = 3, min.features = 200)


#Merge differenct tissue ILC to one dataset
ILC_total<- merge(ILC_wt1, y = c(ILC_wt2, ILC_hdm1, ILC_hdm2), cell.ids = c("PBS1", "PBS2", "HDM1", "HDM2"), project = "lung ILC")


# What does the matrix look like?
# Lets examine a few genes in the first hundred cells
wt1 <- ILC.data_wt1[c("Il5", "Il13", "Chat"), 1:100]
wt2 <- ILC.data_wt2[c("Il5", "Il13", "Chat"), 1:100]
hdm1 <- ILC.data_hdm1[c("Il5", "Il13", "Chat"), 1:100]
hdm2 <- ILC.data_hdm2[c("Il5", "Il13", "Chat"), 1:100]


#[Optional]Export the raw counts. Can take a look by Excel
write.csv(wt1, "~/Workspace/GSE102299_RAW/raw_wt1.csv", row.names = TRUE)
write.csv(wt2, "~/Workspace/GSE102299_RAW/raw_wt2.csv", row.names = TRUE)
write.csv(hdm1, "~/Workspace/GSE102299_RAW/raw_hdm1.csv", row.names = TRUE)
write.csv(hdm2, "~/Workspace/GSE102299_RAW/raw_hdm2.csv", row.names = TRUE)
#Know the total number of cells in each file. For example, hdm1 have 3344 cells 

#[Optional]Examine the total raw counts of some gene of interest
x <- ILC.data_hdm1[c("Il5", "Chat"), 1:3344]
#Inspect the type of the data x
class(x)
#Change to regular matrix to be eligible for calculation
as.matrix(x)
class(x)
x
rowSums(x, na.rm = FALSE)

#Inspect the data matrix size
dense.size_lung <- object.size(as.matrix(ILC.data_wt1))
dense.size_lung

sparse.size_lung<- object.size(ILC.data_wt1)
sparse.size_lung

#QC and selecting 
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
ILC_total[["percent.mt"]] <- PercentageFeatureSet(ILC_total, pattern = "^MT-")

#Inspect QC metric
#The number of unique genes and total molecules are automatically calculated during CreateSeuratObject
#You can find them stored in the object meta data
# Show QC metrics for the first 5 cells
head(ILC_total@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(ILC_total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(ILC_total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ILC_total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#select the high quality one
ILC_total <- subset(ILC_total, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 5)
# Numbers can be changed between different batch of expeirment. Try to cut the high of feature

#NORMALIZE
#By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.
ILC_total_norm <- NormalizeData(ILC_total, normalization.method = "LogNormalize", scale.factor = 10000)


#Identification of highly variable features (feature selection)
ILC_total_norm <- FindVariableFeatures(ILC_total_norm, selection.method = "vst", nfeatures = 2000)


# Identify the 20 most highly variable genes
top20_norm <- head(VariableFeatures(ILC_total_norm), 20)
top20_norm
# plot variable features with and without labels
plot3 <- VariableFeaturePlot(ILC_total_norm)
plot3
plot4 <- LabelPoints(plot = plot3, points = top20_norm, repel = TRUE)
plot4

#Scaling the data; apply linear transformation that is a preprocessing step prior to PCA
#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in ILC_total[["RNA"]]@scale.data
all.genes <- rownames(ILC_total_norm)
ILC_total_norm <- ScaleData(ILC_total_norm, features = all.genes)

#Perform linear dimensional reduction
ILC_total_norm <- RunPCA(ILC_total_norm, features = VariableFeatures(object = ILC_total_norm))

#Seurat provides several useful ways of visualizing both cells and features that define the PCA, including VizDimReduction, DimPlot, and DimHeatmap
#I prefer #2
#1:VizDimLoadings
VizDimLoadings(ILC_total_norm, dims = 1:2, reduction = "pca")
#2
DimPlot(ILC_total_norm, reduction = "pca")
#3
DimHeatmap(ILC_total_norm, dims = 1:3, cells = 500, balanced = TRUE)

#Cluster the cells
#As in PhenoGraph, we first construct a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity). This step is performed using the FindNeighbors function, and takes as input the previously defined dimensionality of the dataset (first 10 PCs)
ILC_total_norm <- FindNeighbors(ILC_total_norm, dims = 1:10)
# Check at cluster IDs of the first 5 cells
head(Idents(ILC_total_norm), 5)

#non-linear dimensional reduction (UMAP/tSNE)
# UMAP is better with both local and global reduction
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages ='umap-learn')
ILC_total_norm <- RunUMAP(ILC_total_norm, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(ILC_total_norm, reduction = "umap", cols = c("red", "orange", "blue", "purple"))
# GOT UMAP!!!

#plot your gene of interest
VlnPlot(ILC_total_norm, features = c("Il5", "Il13", "Areg", "Chat"))
VlnPlot(ILC_total_norm, features = c("Gata3", "Klrg1", "Il1rl1", "Il17rb", "Tcf12", "Tcf3", "Tcf4", "Id1", "Id2", "Id3", "Id4"))

FeaturePlot(ILC_total_norm, features = c("Gata3", "Il7r", "Il1rl1", "Il17rb", "Klrg1", "Alox5","Ltb4r1", "Cysltr1", "Cysltr2", "Lta4h", "Ltc4s"))
FeaturePlot(ILC_total_norm, features = c("Il5", "Il13", "Areg", "Arg1", "Ptgs1", "Ptgs2","Hpgds", "Ptgdr1", "Ptgdr2"))
FeaturePlot(ILC_total_norm, features = c("Chat", "Chrm3", "Chrna1", "Chrna1os", "Chrna7", "Chrna9", "Chrnb1","Chrnb2", "Chrne"))
FeaturePlot(ILC_total_norm, features = c("Gata3","Il1rl1", "Il17rb", "Klrg1", "Tcf12", "Tcf3", "Tcf4", "Id1", "Id2", "Id3", "Id4"))
FeaturePlot(ILC_total_norm, features = c("Klrg1","Arg1", "Ccl5", "Il13", "Amica1", "Cd7", "Ms4a4b", "Klrk1", "Klf2", "Cd3d", "Cd3e", "Stmn1", "Il18r1", "Cd2", "Tcf7", "Ccl5", "Nkg7", "Ms4a6b"))

#Differential analysis
ILC.markers <- FindAllMarkers(ILC_total_norm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ILC.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_logFC)

#Save the list
write.csv(ILC.markers, "~/Workspace/GSE102299_RAW/ILC_All_marker.csv", row.names = TRUE)

#[Optional] Let the machine differentiate different cell types
#To cluster the cells, we next apply modularity optimization techniques such as the Louvain algorithm (default) or SLM [SLM, Blondel et al., Journal of Statistical Mechanics], to iteratively group cells together, with the goal of optimizing the standard modularity function. The FindClusters function implements this procedure
#We find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells
ILC_total_norm <- FindClusters(ILC_total_norm, resolution = 0.5)
# resolution = 0.4 have less groups, 1.2 have more groups
# Look at cluster IDs of the first 5 cells
head(Idents(ILC_total_norm), 5)

#non-linear dimensional reduction (UMAP/tSNE)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages ='umap-learn')
ILC_total_norm <- RunUMAP(ILC_total_norm, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(ILC_total_norm, reduction = "umap")


#Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 1
#it identifes positive and negative markers of a single cluster (specified in ident.1
#The min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells, and the thresh.test argument requires a feature to be differentially expressed (on average) by some amount between the two groups.
cluster1.markers <- FindMarkers(ILC_total_norm, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
#same code can be applied to cluster 2, 3, ..etc

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(ILC_total_norm, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)


#Diiferential analysis of each cluster
ILC.markers_2 <- FindAllMarkers(ILC_total_norm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ILC.markers_2 %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_logFC)
#Save the list
write.csv(ILC.markers_2, "~/Workspace/GSE102299_RAW/ILC_All_marker_cluster.csv", row.names = TRUE)

#Method 1 : draw heatmap on selected DEGs
##Draw the heatmap of Klrg1+- DEGs!!!!!!
features <- c("Klrg1","Arg1", "Ccl5", "Il13", "Amica1", "Cd7", "Ms4a4b", "Klrk1", "Klf2", "Cd3d", "Cd3e", "Stmn1", "Il18r1", "Cd2", "Tcf7", "Ccl5", "Nkg7", "Ms4a6b")
#Load the library
library(ggplot2)
DoHeatmap(ILC_total_norm, features = features) + NoLegend()

#Method 2 : draw heatmap of all DEGs
features <- rownames(ILC.markers)
#To designate colors in gradient
DoHeatmap(ILC_total_norm , features = features) + scale_fill_gradientn(colors = c("purple", "black", "yellow"))


