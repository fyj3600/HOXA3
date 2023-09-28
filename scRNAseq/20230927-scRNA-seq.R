# Code: Integrated analysis of multiple samples
#Compiled: September 27, 2023
#The provider：Xueyan Zhang
#This code describes the analysis process of single cell sequencing results in the article Yingjie Fu is submitting.
### Setup the Seurat objects
#Single-cell sequencing data (GSE182135) from the GEO database were downloaded locally.  We first read in the ten samples and set up the Seurat objects. 
#Single-cell RNA profiling of the pharyngeal endoderm of the mouse embryo. On embryonic days 9.5, 10.5, 11.5, and 12.5, pharyngeal endoderm was isolated by a combination of dissection (mechanical isolation of the whole pharynx) and cell sorting (purification of viable, Pax9+ epithelia). Data include two or three samples per timepoint (biological replicates).
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
# The hypervariable genes screened in each dataset were extracted.
install.packages("venn") 
library(venn) 
alldata.list<- c(E9_5_rep2,E9_5_rep1,E10_5_rep7,E10_5_rep6,E10_5_rep5,E11_5_rep3,
                 E11_5_rep2B,E11_5_rep2,E12_5_rep2,E12_5_rep1)
hvgs_per_dataset <- lapply(alldata.list, function(x) {
  x@assays$RNA@var.features
})
head(hvgs_per_dataset)
venn::venn(hvgs_per_dataset, opacity = 0.4, zcolor = (scales::hue_pal())(3), cexsn = 1, 
           cexil = 1, lwd = 1, col = "white", frame = F, borders = NA)


### Perform integration
#We then identify anchors using the function, which takes a list of Seurat objects as input, and use these anchors to integrate the two datasets together with .`FindIntegrationAnchors``IntegrateData`

# Integration
#We then identify anchors using the function, which takes a list of Seurat objects as input, and use these anchors to integrate the two datasets together with
testAB.anchors <- FindIntegrationAnchors(object.list = list(E9_5_rep2,
                                                            E9_5_rep1,E12_5_rep2,E12_5_rep1,E11_5_rep3,E11_5_rep2B,
                                                            E11_5_rep2,E10_5_rep7,E10_5_rep6,E10_5_rep5), dims = 1:30)
testAB.integrated <- IntegrateData(anchorset = testAB.anchors, dims = 1:30,
                                   new.assay.name = "CCA")

### Perform an integrated analysis
#Now we can run a single integrated analysis on all cells!
##Standard procedure.
## Run the standard workflow for visualization and clustering.
testAB.integrated <- ScaleData(testAB.integrated, features = rownames(testAB.integrated))
testAB.integrated <- RunPCA(testAB.integrated, npcs = 100, verbose = FALSE)
ElbowPlot(testAB.integrated)
DimHeatmap(testAB.integrated, dims = 1:100, cells = 500, balanced = TRUE)
library(data.table)
class(testAB.integrated)
testAB.integrated <- FindNeighbors(testAB.integrated, dims = 1:50)
testAB.integrated <- FindClusters(testAB.integrated, resolution = 2)
## t-SNE and Clustering
testAB.integrated <- RunUMAP(testAB.integrated, dims = 1:50)
testAB.integrated <- RunTSNE(testAB.integrated, dims = 1:50)
testAB.integrated$patient=str_replace(testAB.integrated$orig.ident,"_.*$","")
#Save the file
saveRDS(testAB.integrated,file = "seurat-all.rds")
testAB.integrated <- readRDS(file = "seurat-all.rds")

## Visualization
p1 <- DimPlot(testAB.integrated, reduction = "tsne", group.by = "orig.ident", 
              pt.size=0.5) 
p2 <- DimPlot(testAB.integrated, reduction = "tsne", group.by = "ident",
              pt.size=0.5, label = TRUE,repel = TRUE) 
fig_umap <- plot_grid(p1, p2, labels = c('orig.ident','ident'),align = "v",ncol = 2) 
ggsave(filename = "seurat_all_tsne.pdf", plot = fig_umap, device = 'pdf', 
       width = 30, height = 15, units = 'cm')

p3 <- DimPlot(testAB.integrated, reduction = "umap", group.by = "orig.ident",
              pt.size=0.5) 
p4 <- DimPlot(testAB.integrated, reduction = "umap", group.by = "ident", 
              pt.size=0.5, label = TRUE,repel = TRUE) 
fig_umap <- plot_grid(p3, p4, labels = c('patient','ident'),align = "v",ncol = 2) 
ggsave(filename = "seurat_all_umap.pdf", plot = fig_umap, device = 'pdf', 
       width = 30, height = 15, units = 'cm')

celltype_marker=c("Pax9", "Neurod1","Hoxa3","Pdgfra"," Dlx5","Esam",
                  "Epcam","Rxrg","Sox10","Fstl1")
#Visualization
VlnPlot(testAB.integrated,features = celltype_marker,pt.size = 0,ncol = 2)
table(testAB.integrated@meta.data$seurat_clusters)
#Cell populations with low pax9 and EPCAM gene expression were extracted
cd4_sce1 = testAB.integrated[,!(testAB.integrated@meta.data$seurat_clusters %in% c(6,26,39,40,43,45))]
#Run the standard workflow for visualization and clustering.
cd4_sce1 <- ScaleData(cd4_sce1, features = rownames(cd4_sce1))
cd4_sce1 <- RunPCA(cd4_sce1, features = VariableFeatures(cd4_sce1),npcs = 100)
DimHeatmap(cd4_sce1, dims =1:30, cells = 500, balanced = TRUE)
ElbowPlot(cd4_sce1, ndims =100)
cd4_sce1 <- FindNeighbors(cd4_sce1, dims = 1:60)
cd4_sce1 <- FindClusters(cd4_sce1, resolution = 0.5)
cd4_sce1 <- RunUMAP(cd4_sce1, dims = 1:60)
cd4_sce1 <- RunTSNE(cd4_sce1, dims = 1:60)
## Visualization
p1 <- DimPlot(cd4_sce1, reduction = "tsne", group.by = "orig.ident", pt.size=0.5)
p2 <- DimPlot(cd4_sce1, reduction = "tsne", group.by = "ident",   pt.size=0.5,
              label = TRUE,repel = TRUE)
fig_tsne <- plot_grid(p1, p2, labels = c('orig.ident','ident'),align = "v",ncol = 2)
fig_tsne

p3 <- DimPlot(cd4_sce1, reduction = "umap", group.by = "orig.ident", pt.size=0.5) 
p4 <- DimPlot(cd4_sce1, reduction = "umap", group.by = "ident", pt.size=0.5,
              label = TRUE,repel = TRUE) 
fig_umap <- plot_grid(p3, p4, labels = c('orig.ident','ident'),align = "v",ncol = 2) 
fig_umap


### Cell Type Annotation

#We can explore these marker genes for each cluster and use them to annotate our clusters as specific cell types.

###Identify conserved cell type markers.
TPP_marker <- c("Eya1","Pax9","Hoxa3","Foxn1","Bmp4","Six1","Pax1","Tbx1")
#Violin diagram
VlnPlot(cd4_sce1,features = TPP_marker,pt.size = 0,ncol = 2)
table(cd4_sce1@meta.data$seurat_clusters)
Third_pharyngeal_pouch <- c(0,1,5,9,16,20)
Unknown <- c(2,3,4,6,7,8,10,11,12,13,14,15,17,18,19,21,22,23,24)
current.cluster.ids <- c(Third_pharyngeal_pouch,
                         Unknown)
new.cluster.ids <- c(rep("Third pharyngeal pouch",length(Third_pharyngeal_pouch)),
                     rep("Unknown",length(Unknown)))

cd4_sce1@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character
                                                              (cd4_sce1@meta.data$seurat_clusters)), 
                                               from = current.cluster.ids, to = new.cluster.ids)
head(cd4_sce1@meta.data)
#tsne
plotCB=as.data.frame(cd4_sce1@meta.data%>%filter(seurat_clusters!="13" &
                                                   seurat_clusters!="15"))[,"CB"]
colors <-  c("#FFD700","#DCDCDC")
DimPlot(cd4_sce1, reduction = "tsne", group.by = "celltype", pt.size=0.5,cols = colors,label = TRUE,repel = TRUE)

p1 <- DimPlot(cd4_sce1, reduction = "tsne", group.by = "orig.ident", pt.size=0.5)
p2 <- DimPlot(cd4_sce1, reduction = "tsne", group.by = "celltype",
              pt.size=0.5,cols = colors,repel = TRUE)
fig_tsne <- plot_grid(p1, p2, labels = c('orig.ident','celltype'),align = "v",ncol = 2)
fig_tsne

#We extracted a subset of Third pharyngeal pouch
TPP <- cd4_sce1[,(cd4_sce1@meta.data$seurat_clusters %in% c(0,1,5,9,16,20))]
table(TPP@meta.data$seurat_clusters)
Idents(TPP) <- "orig.ident" 

vln_df = TPP[,TPP@meta.data$orig.ident %in% c("E9_5_rep1",'E10_5_rep7',"E11_5_rep2","E12_5_rep2")]
table(vln_df@meta.data$orig.ident)
p1 <- VlnPlot(subset(vln_df,Hoxa3 > 0.2 ),"Hoxa3",adjust = .8,pt.size=0)+ theme_bw()+
  geom_boxplot(width=.2,col="black",fill="white")+  NoLegend()
p2 <- VlnPlot(subset(vln_df,Eya1 > 0.2 ),"Eya1",adjust = .8,pt.size=0)+ theme_bw()+
  geom_boxplot(width=.2,col="black",fill="white")+  NoLegend()
p3 <- VlnPlot(subset(vln_df,Six1 > 0.2 ),"Six1",adjust = .8,pt.size=0)+ theme_bw()+
  geom_boxplot(width=.2,col="black",fill="white")+  NoLegend()
p4 <- VlnPlot(subset(vln_df,Epcam > 0.2 ),"Epcam",adjust = .8,pt.size=0)+ theme_bw()+
  geom_boxplot(width=.2,col="black",fill="white") +NoLegend()
p5 <- VlnPlot(subset(vln_df,Pax9 > 0.2 ),"Pax9",adjust = .8,pt.size=0)+ theme_bw()+
  geom_boxplot(width=.2,col="black",fill="white") +NoLegend()

library(patchwork)
(p1|p2)/(p3|p4)/(p5|plot_spacer())
ggsave(filename = "TPP_Violin .pdf", width = 22, height = 20, 
       device = 'pdf', units = 'cm')

#Figure1_A和Figure1_D
gene <- c("Epcam","Pax9","Bmp4","Hoxa3","Pax1","Ptprc","Six1","Tbx1") 
for(i in gene){
  print(paste0("the current column is ",i))
  FeaturePlot(cd4_sce1,features = i,reduction = "tsne", min.cutoff = "q10", max.cutoff = "q90", pt.size = 1, repel = F,label = F,
              order = T,cols = c("#DCDCDC", "#DC143C"))+
    scale_x_continuous("")+scale_y_continuous("")+
    theme_bw()+ 
    theme( 
      panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
      axis.ticks = element_blank(),axis.text = element_blank(), 
      plot.title = element_text(hjust = 0.5,size=14)
    )
  ggsave(filename = paste0(i,"_tsne.pdf"), width = 10, height = 8, 
         device = 'pdf', units = 'cm')
}

#Figure1_E
F1_E <- c("Aire","Cldn4","Cldn3","Foxn1","Ly75","Psmb11","Krt5","Krt8") 
for(i in F1_E){
  print(paste0("the current column is ",i))
  FeaturePlot(cd4_sce1,features = i,reduction = "tsne", min.cutoff = "q10", max.cutoff = "q90", pt.size = 1, repel = F,label = F,
              order = T,cols = c('#D5F0E2','#4682B4'))+
    scale_x_continuous("")+scale_y_continuous("")+
    theme_bw()+ 
    theme( 
      panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
      axis.ticks = element_blank(),axis.text = element_blank(), 
      plot.title = element_text(hjust = 0.5,size=14)
    )
  ggsave(filename = paste0(i,"_tsne.pdf"), width = 10, height = 8, 
         device = 'pdf', units = 'cm')
}

