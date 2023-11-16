# Code: Integrated analysis of multiple samples
#Compiled: November 16, 2023
#The provider：Xueyan Zhang
#This code describes the analysis process of single cell sequencing results in the article Yingjie Fu is submitting.
### Setup the Seurat objects
#Single-cell sequencing data (GSE182135) from the GEO database were downloaded locally.  We first read in the ten samples and set up the Seurat objects. 
#Single-cell RNA profiling of the pharyngeal endoderm of the mouse embryo. 
#On embryonic days 9.5, 10.5, 11.5, and 12.5, pharyngeal endoderm was isolated by a combination of dissection (mechanical 
#isolation of the whole pharynx) and cell sorting (purification of viable, Pax9+ epithelia). Data include two or three samples per timepoint (biological replicates).
#Ten samples were integrated. 
library("scde")
library("DEsingle")
library("Seurat")
library("dplyr")
library("magrittr")
library("SingleCellExperiment")
library("patchwork")
library("cowplot")
library("stringr")
library("tidyverse")
library("export")
# step1:Quality control and filtration of 10 samples.
#Set up object
#1
E9_5_rep2 <- Read10X( data.dir="E9_5_rep2_2018AUG07/")
E9_5_rep2 <- CreateSeuratObject(counts = E9_5_rep2, project = "E9_5_rep2", min.cells = 3, min.features = 200)
#The proportion of mitochondrial RNA was calculated
E9_5_rep2[["percent.mt"]] <- PercentageFeatureSet(E9_5_rep2, pattern = "^MT-") 
VlnPlot(E9_5_rep2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
## Data filtering
E9_5_rep2 <- subset(E9_5_rep2, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 1)
dim(E9_5_rep2)
#2
E9_5_rep1 <- Read10X( data.dir="E9_5_rep1_2018JULY18/")
E9_5_rep1 <- CreateSeuratObject(counts = E9_5_rep1, project = "E9_5_rep1", min.cells = 3, min.features = 200)
E9_5_rep1[["percent.mt"]] <- PercentageFeatureSet(E9_5_rep1, pattern = "^MT-")
VlnPlot(E9_5_rep1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E9_5_rep1 <- subset(E9_5_rep1, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 1)
dim(E9_5_rep1)
#3
E10_5_rep7 <- Read10X( data.dir="E10_5_rep7_2018AUG16/")
E10_5_rep7 <- CreateSeuratObject(counts = E10_5_rep7, project = "E10_5_rep7", min.cells = 3, min.features = 200)
E10_5_rep7[["percent.mt"]] <- PercentageFeatureSet(E10_5_rep7, pattern = "^MT-")
VlnPlot(E10_5_rep7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E10_5_rep7 <- subset(E10_5_rep7, subset = nFeature_RNA > 1000 & nFeature_RNA < 6500 & percent.mt < 1)
dim(E10_5_rep7)
#4
E10_5_rep6 <- Read10X( data.dir="E10_5_rep6_2018SEP22/")
E10_5_rep6 <- CreateSeuratObject(counts = E10_5_rep6, project = "E10_5_rep6", min.cells = 3, min.features = 200)
E10_5_rep6[["percent.mt"]] <- PercentageFeatureSet(E10_5_rep6, pattern = "^MT-")
VlnPlot(E10_5_rep6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E10_5_rep6 <- subset(E10_5_rep6, subset = nFeature_RNA > 1000 & nFeature_RNA < 6500 & percent.mt < 1)
dim(E10_5_rep6)
#5
E10_5_rep5 <- Read10X( data.dir="E10_5_rep5_2018AUG16/")
E10_5_rep5 <- CreateSeuratObject(counts = E10_5_rep5, project = "E10_5_rep5", min.cells = 3, min.features = 200)
E10_5_rep5[["percent.mt"]] <- PercentageFeatureSet(E10_5_rep5, pattern = "^MT-")
VlnPlot(E10_5_rep5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E10_5_rep5 <- subset(E10_5_rep5, subset = nFeature_RNA > 1000 & nFeature_RNA < 6500 & percent.mt < 1)
dim(E10_5_rep5)
#6
E11_5_rep3 <- Read10X( data.dir="E11_5_rep3_2018SEP22/")
E11_5_rep3 <- CreateSeuratObject(counts = E11_5_rep3, project = "E11_5_rep3", min.cells = 3, min.features = 200)
E11_5_rep3[["percent.mt"]] <- PercentageFeatureSet(E11_5_rep3, pattern = "^MT-")
VlnPlot(E11_5_rep3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E11_5_rep3 <- subset(E11_5_rep3, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 1)
dim(E11_5_rep3)
#7
E11_5_rep2B <- Read10X( data.dir="E11_5_rep2B_2018JULY27/")
E11_5_rep2B <- CreateSeuratObject(counts = E11_5_rep2B, project = "E11_5_rep2B", min.cells = 3, min.features = 200)
E11_5_rep2B[["percent.mt"]] <- PercentageFeatureSet(E11_5_rep2B, pattern = "^MT-")
VlnPlot(E11_5_rep2B, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E11_5_rep2B <- subset(E11_5_rep2B, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 1)
dim(E11_5_rep2B)
#8
E11_5_rep2 <- Read10X( data.dir="E11_5_rep2_2018JUNE26/")
E11_5_rep2 <- CreateSeuratObject(counts = E11_5_rep2, project = "E11_5_rep2", min.cells = 3, min.features = 200)
E11_5_rep2[["percent.mt"]] <- PercentageFeatureSet(E11_5_rep2, pattern = "^MT-")
VlnPlot(E11_5_rep2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E11_5_rep2 <- subset(E11_5_rep2, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 1)
dim(E11_5_rep2)
#9
E12_5_rep2 <- Read10X( data.dir="E12_5_rep2_2018JUNE28/")
E12_5_rep2 <- CreateSeuratObject(counts = E12_5_rep2, project = "E12_5_rep2", min.cells = 3, min.features = 200)
E12_5_rep2[["percent.mt"]] <- PercentageFeatureSet(E12_5_rep2, pattern = "^MT-")
VlnPlot(E12_5_rep2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E12_5_rep2 <- subset(E12_5_rep2, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 1)
dim(E12_5_rep2)
#10
E12_5_rep1 <- Read10X( data.dir="E12_5_rep1_2018JUNE28/")
E12_5_rep1 <- CreateSeuratObject(counts = E12_5_rep1, project = "E12_5_rep1", min.cells = 3, min.features = 200)
E12_5_rep1[["percent.mt"]] <- PercentageFeatureSet(E12_5_rep1, pattern = "^MT-")
VlnPlot(E12_5_rep1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E12_5_rep1 <- subset(E12_5_rep1, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 1)
dim(E12_5_rep1)
##########################################################
######Merge 10 Seurat objects 
test.seu <- merge(E9_5_rep2, y = c(E9_5_rep1,E10_5_rep7,E10_5_rep6,E10_5_rep5,E11_5_rep3,
                                   E11_5_rep2B,E11_5_rep2,E12_5_rep2,E12_5_rep1), project = "ALL")
test.seu
head(colnames(test.seu))
unique(sapply(X = strsplit(colnames(test.seu), split = "_"), FUN = "[", 1))
table(test.seu$orig.ident)

#Now we can run a single integrated analysis on all cells!
##Standard procedure.
## Run the standard workflow for visualization and clustering.
test.seu <- NormalizeData(test.seu, normalization.method = "LogNormalize", scale.factor = 10000)
test.seu <- FindVariableFeatures(test.seu, selection.method = "vst", nfeatures = 2000)
test.seu <- ScaleData(test.seu, features = rownames(test.seu))
test.seu <- RunPCA(test.seu, features = VariableFeatures(test.seu),npcs = 100)
ElbowPlot(test.seu)
## t-SNE and Clustering
test.seu <- FindNeighbors(test.seu, dims = 1:40)
test.seu <- FindClusters(test.seu, resolution = 1.5)
test.seu <- RunUMAP(test.seu, dims = 1:40)
test.seu <- RunTSNE(test.seu, dims = 1:40)

## Visualization of umap
p3 <- DimPlot(test.seu, reduction = "umap", group.by = "orig.ident", pt.size=0.5) 
p4 <- DimPlot(test.seu, reduction = "umap", group.by = "ident", pt.size=0.5, label = TRUE,repel = TRUE) 
fig_umap <- plot_grid(p3, p4, labels = c('orig.ident','ident'),align = "v",ncol = 2) 
fig_umap
ggsave(filename = "all_umap.pdf", plot = fig_umap, device = 'pdf', width = 30, height = 15, units = 'cm')

#Save the file
saveRDS(test.seu,file = "2-1012umap.rds")
test.seu <- readRDS(file = "2-1012umap.rds")

#Visualization
celltype_marker=c("Pax9", "Neurod1","Hoxa3","Pdgfra",
                  "Epcam","Rxrg","Sox10","Fstl1")
VlnPlot(test.seu,features = celltype_marker,pt.size = 0,ncol = 2)
ggsave(filename = "all-marker.png",device = "png",width = 44,height = 33,units = "cm")

#Cell populations with high pax9 and EPCAM gene expression were extracted
cd4_sce1 = test.seu[,!(test.seu@meta.data$seurat_clusters %in% c(12,23,27,29,31,32,36,38,39,41,42,44))]
#Run the standard workflow for visualization and clustering.
cd4_sce1 <- NormalizeData(cd4_sce1, normalization.method = "LogNormalize", scale.factor = 10000)
cd4_sce1 <- FindVariableFeatures(cd4_sce1, selection.method = "vst", nfeatures = 2000)
cd4_sce1 <- ScaleData(cd4_sce1, features = rownames(cd4_sce1))
cd4_sce1 <- RunPCA(cd4_sce1, features = VariableFeatures(cd4_sce1),npcs = 100)
DimHeatmap(cd4_sce1, dims =1:30, cells = 500, balanced = TRUE)
ElbowPlot(cd4_sce1, ndims =100)
cd4_sce1 <- FindNeighbors(cd4_sce1, dims = 1:60)
cd4_sce1 <- FindClusters(cd4_sce1, resolution = 0.8)
cd4_sce1 <- RunUMAP(cd4_sce1, dims = 1:60)
cd4_sce1 <- RunTSNE(cd4_sce1, dims = 1:60)
#Save the file
saveRDS(cd4_sce1,file = "1017-60pctest.rds")
## Visualization
cd4_sce1@meta.data$orig.ident <- factor(cd4_sce1@meta.data$orig.ident, ordered = TRUE, 
                                        levels = c("E9_5_rep1","E9_5_rep2","E10_5_rep5","E10_5_rep6","E10_5_rep7","E11_5_rep2",
                                                    "E11_5_rep2B","E11_5_rep3","E12_5_rep1","E12_5_rep2"))
p1 <- DimPlot(cd4_sce1, reduction = "umap", group.by = "orig.ident", pt.size=0.5) 
p1
p2 <- DimPlot(cd4_sce1, reduction = "umap", group.by = "ident", pt.size=0.5, label = TRUE,repel = TRUE) 
p2
fig_umap <- plot_grid(p1, p2, labels = c('orig.ident','ident'),align = "v",ncol = 2) 
fig_umap
ggsave(filename = "umap.pdf", plot = fig_umap, device = 'pdf', width = 40, height = 15, units = 'cm')

### Cell Type Annotation
#We can explore these marker genes for each cluster and use them to annotate our clusters as specific cell types.
###Identify conserved cell type markers.
PPE_marker <- c("Eya1","Pax9","Hoxa3","Foxn1","Bmp4","Six1","Pax1","Tbx1","Gcm2","Epcam")
#Violin diagram
VlnPlot(cd4_sce1,features = PPE_marker,adjust=1,pt.size = 0.1,ncol = 2)
ggsave(filename = "Annotation-Vln.pdf", device = 'pdf', width = 40, height = 25, units = 'cm')
#manual annotation
table(cd4_sce1@meta.data$seurat_clusters)
Third_pharyngeal_pouch <- c(8,13,14)
Unknown <- c(0:7,9:12,15:26)
current.cluster.ids <- c(Third_pharyngeal_pouch,
                         Unknown)
new.cluster.ids <- c(rep("Third pharyngeal pouch",length(Third_pharyngeal_pouch)),
                     rep("Unknown",length(Unknown)))

cd4_sce1@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(cd4_sce1@meta.data$seurat_clusters)), 
                                               from = current.cluster.ids, to = new.cluster.ids)
table(cd4_sce1@meta.data$celltype)#2366 3PPE 
#umap
colors <-  c("#F2C11D","#BEBEBE")
p2 <- DimPlot(cd4_sce1, reduction = "umap", group.by = "celltype", pt.size=0.5,cols = colors,repel = TRUE)
ggsave(filename = "3PPE_umap.pdf", plot = p2, device = 'pdf', width = 18, height = 12, units = 'cm')

#We extracted the subset of Third pharyngeal pouch
PPE3 <- cd4_sce1[,(cd4_sce1@meta.data$seurat_clusters %in% c(8,13,14))]
table(PPE3@meta.data$orig.ident)
Idents(PPE3) <- "orig.ident"  
#Select subgroups for mapping
#Figure1_F
vln_df = PPE3[,PPE3@meta.data$orig.ident %in% c("E9_5_rep1",'E10_5_rep5',"E11_5_rep2","E12_5_rep1")]
VlnPlot(vln_df,features = c("Hoxa3","Eya1","Six1","Pax9","Epcam"),pt.size = 2,adjust = 1,ncol = 2)
ggsave(filename = "3PPE-vln.pdf", plot = p2, device = 'pdf', width = 35, height = 25, units = 'cm')

#和Figure1_D
F1_AD <- c("Epcam","Pax9","Bmp4","Hoxa3","Pax1","Ptprc","Six1","Tbx1") 
for(i in F1_AD){
  print(paste0("the current column is ",i))
  FeaturePlot(cd4_sce1,features = i,reduction = "umap", min.cutoff = "q10", max.cutoff = "q90", pt.size = 1, repel = F,label = F,
              order = T,cols = c("#DCDCDC", "#DC143C"))+
    scale_x_continuous("")+scale_y_continuous("")+
    theme_bw()+ 
    theme( 
      panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
      axis.ticks = element_blank(),axis.text = element_blank(), 
      plot.title = element_text(hjust = 0.5,size=14)
    )
  ggsave(filename = paste0(i,"_umap.pdf"), width = 10, height = 8, 
         device = 'pdf', units = 'cm')
}

#Figure1_E
F1_E <- c("Aire","Cldn4","Cldn3","Foxn1","Ly75","Psmb11","Krt5","Krt8") 
for(i in F1_E){
  print(paste0("the current column is ",i))
  FeaturePlot(cd4_sce1,features = i,reduction = "umap", min.cutoff = "q10", max.cutoff = "q90", pt.size = 1, repel = F,label = F,
              order = T,cols = c('#D5F0E2','#4682B4'))+
    scale_x_continuous("")+scale_y_continuous("")+
    theme_bw()+ 
    theme( 
      panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
      axis.ticks = element_blank(),axis.text = element_blank(), 
      plot.title = element_text(hjust = 0.5,size=14)
    )
  ggsave(filename = paste0(i,"_umap.pdf"), width = 10, height = 8, 
         device = 'pdf', units = 'cm')
}

