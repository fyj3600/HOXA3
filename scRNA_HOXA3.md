```R
# yingjie fu
#GEO:GSE182135, 10个样本合并
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

# step1:导入10个样本数据,并过滤
setwd("E:/scRNA-seq/20221010/原始数据/1")
rm(list=ls())#清空环境变量
#1
E9_5_rep2 <- Read10X( data.dir="E9_5_rep2_2018AUG07/")
#将矩阵转换为Seurat对象
E9_5_rep2 <- CreateSeuratObject(counts = E9_5_rep2, project = "E9_5_rep2", min.cells = 3, min.features = 200)
E9_5_rep2[["percent.mt"]] <- PercentageFeatureSet(E9_5_rep2, pattern = "^MT-") #计算线粒体RNA比例
VlnPlot(E9_5_rep2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) #数值分布的小提琴图
## 过滤
E9_5_rep2 <- subset(E9_5_rep2, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 1)
dim(E9_5_rep2)
#2
E9_5_rep1 <- Read10X( data.dir="E9_5_rep1_2018JULY18/")
E9_5_rep1 <- CreateSeuratObject(counts = E9_5_rep1, project = "E9_5_rep1", min.cells = 3, min.features = 200)
E9_5_rep1[["percent.mt"]] <- PercentageFeatureSet(E9_5_rep1, pattern = "^MT-")
VlnPlot(E9_5_rep1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E9_5_rep1 <- subset(E9_5_rep1, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 1)
dim(E9_5_rep1)
#8
E10_5_rep7 <- Read10X( data.dir="E10_5_rep7_2018AUG16/")
E10_5_rep7 <- CreateSeuratObject(counts = E10_5_rep7, project = "E10_5_rep7", min.cells = 3, min.features = 200)
E10_5_rep7[["percent.mt"]] <- PercentageFeatureSet(E10_5_rep7, pattern = "^MT-")
VlnPlot(E10_5_rep7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E10_5_rep7 <- subset(E10_5_rep7, subset = nFeature_RNA > 1000 & nFeature_RNA < 6500 & percent.mt < 1)
dim(E10_5_rep7)
#9
E10_5_rep6 <- Read10X( data.dir="E10_5_rep6_2018SEP22/")
E10_5_rep6 <- CreateSeuratObject(counts = E10_5_rep6, project = "E10_5_rep6", min.cells = 3, min.features = 200)
E10_5_rep6[["percent.mt"]] <- PercentageFeatureSet(E10_5_rep6, pattern = "^MT-")
VlnPlot(E10_5_rep6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E10_5_rep6 <- subset(E10_5_rep6, subset = nFeature_RNA > 1000 & nFeature_RNA < 6500 & percent.mt < 1)
dim(E10_5_rep6)
#10
E10_5_rep5 <- Read10X( data.dir="E10_5_rep5_2018AUG16/")
E10_5_rep5 <- CreateSeuratObject(counts = E10_5_rep5, project = "E10_5_rep5", min.cells = 3, min.features = 200)
E10_5_rep5[["percent.mt"]] <- PercentageFeatureSet(E10_5_rep5, pattern = "^MT-")
VlnPlot(E10_5_rep5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E10_5_rep5 <- subset(E10_5_rep5, subset = nFeature_RNA > 1000 & nFeature_RNA < 6500 & percent.mt < 1)
dim(E10_5_rep5)
#5
E11_5_rep3 <- Read10X( data.dir="E11_5_rep3_2018SEP22/")
E11_5_rep3 <- CreateSeuratObject(counts = E11_5_rep3, project = "E11_5_rep3", min.cells = 3, min.features = 200)
E11_5_rep3[["percent.mt"]] <- PercentageFeatureSet(E11_5_rep3, pattern = "^MT-")
VlnPlot(E11_5_rep3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E11_5_rep3 <- subset(E11_5_rep3, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 1)
dim(E11_5_rep3)
#6
E11_5_rep2B <- Read10X( data.dir="E11_5_rep2B_2018JULY27/")
E11_5_rep2B <- CreateSeuratObject(counts = E11_5_rep2B, project = "E11_5_rep2B", min.cells = 3, min.features = 200)
E11_5_rep2B[["percent.mt"]] <- PercentageFeatureSet(E11_5_rep2B, pattern = "^MT-")
VlnPlot(E11_5_rep2B, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E11_5_rep2B <- subset(E11_5_rep2B, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 1)
dim(E11_5_rep2B)
#7
E11_5_rep2 <- Read10X( data.dir="E11_5_rep2_2018JUNE26/")
E11_5_rep2 <- CreateSeuratObject(counts = E11_5_rep2, project = "E11_5_rep2", min.cells = 3, min.features = 200)
E11_5_rep2[["percent.mt"]] <- PercentageFeatureSet(E11_5_rep2, pattern = "^MT-")
VlnPlot(E11_5_rep2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E11_5_rep2 <- subset(E11_5_rep2, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 1)
dim(E11_5_rep2)
#3
E12_5_rep2 <- Read10X( data.dir="E12_5_rep2_2018JUNE28/")
E12_5_rep2 <- CreateSeuratObject(counts = E12_5_rep2, project = "E12_5_rep2", min.cells = 3, min.features = 200)
E12_5_rep2[["percent.mt"]] <- PercentageFeatureSet(E12_5_rep2, pattern = "^MT-")
VlnPlot(E12_5_rep2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E12_5_rep2 <- subset(E12_5_rep2, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 1)
dim(E12_5_rep2)
#4
E12_5_rep1 <- Read10X( data.dir="E12_5_rep1_2018JUNE28/")
E12_5_rep1 <- CreateSeuratObject(counts = E12_5_rep1, project = "E12_5_rep1", min.cells = 3, min.features = 200)
E12_5_rep1[["percent.mt"]] <- PercentageFeatureSet(E12_5_rep1, pattern = "^MT-")
VlnPlot(E12_5_rep1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E12_5_rep1 <- subset(E12_5_rep1, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 1)
dim(E12_5_rep1)
##########################################################
# 提取每个数据集中筛选到的高变基因
alldata.list<- c(E9_5_rep2,E9_5_rep1,E10_5_rep7,E10_5_rep6,E10_5_rep5,E11_5_rep3,E11_5_rep2B,E11_5_rep2,E12_5_rep2,E12_5_rep1)
hvgs_per_dataset <- lapply(alldata.list, function(x) {
  x@assays$RNA@var.features
})
head(hvgs_per_dataset)
#install.packages("venn") #安装
library(venn) #导
# 绘制venn图查看不同样本中高变基因的重叠情况
venn::venn(hvgs_per_dataset, opacity = 0.4, zcolor = (scales::hue_pal())(3), cexsn = 1, 
           cexil = 1, lwd = 1, col = "white", frame = F, borders = NA)

### Integration ----
testAB.anchors <- FindIntegrationAnchors(object.list = list(E9_5_rep2,E9_5_rep1,E12_5_rep2,E12_5_rep1,E11_5_rep3,E11_5_rep2B,E11_5_rep2,E10_5_rep7,E10_5_rep6,E10_5_rep5), dims = 1:30)
testAB.integrated <- IntegrateData(anchorset = testAB.anchors, dims = 1:30, new.assay.name = "CCA")
#这一步之后就多了一个整合后的assay（原先有一个RNA的assay），整合前后的数据分别存储在这两个assay中
testAB.integrated ##57717samples
names(testAB.integrated@assays)
## [1] "RNA" "CCA"
# by default, Seurat now sets the integrated assay as the default assay, so any
# operation you now perform will be on the ingegrated data.
testAB.integrated@active.assay
## [1] "CCA"

dim(testAB.integrated[["RNA"]]@counts)
dim(testAB.integrated[["RNA"]]@data)
dim(testAB.integrated[["CCA"]]@counts) #因为是从RNA这个assay的data矩阵开始整合的，所以这个矩阵为空
dim(testAB.integrated[["CCA"]]@data)
##后续仍然是标准流程，基于上面得到的整合data矩阵
# Run the standard workflow for visualization and clustering
testAB.integrated <- ScaleData(testAB.integrated, features = rownames(testAB.integrated))
testAB.integrated <- RunPCA(testAB.integrated, npcs = 100, verbose = FALSE)#npcs : 要计算和存储的PC总数（默认为50
ElbowPlot(testAB.integrated)
DimHeatmap(testAB.integrated, dims = 1:100, cells = 500, balanced = TRUE)
library(data.table)
saveRDS(testAB.integrated,file = "1011seurat-all.rds")
##或者直接读入之前保存的数据
testAB.integrated <- readRDS(file = "1011seurat-all.rds")
class(testAB.integrated)


testAB.integrated <- FindNeighbors(testAB.integrated, dims = 1:50)
testAB.integrated <- FindClusters(testAB.integrated, resolution = 2)
testAB.integrated <- RunUMAP(testAB.integrated, dims = 1:50)
testAB.integrated <- RunTSNE(testAB.integrated, dims = 1:50)

##
#library(cowplot)
#library(stringr)
#library(tidyverse)
testAB.integrated$patient=str_replace(testAB.integrated$orig.ident,"_.*$","")

#后面两个参数用来添加文本标签
p1 <- DimPlot(testAB.integrated, reduction = "tsne", group.by = "orig.ident", pt.size=0.5) 
p1
p2 <- DimPlot(testAB.integrated, reduction = "tsne", group.by = "ident", pt.size=0.5, label = TRUE,repel = TRUE) 
p2
fig_umap <- plot_grid(p1, p2, labels = c('orig.ident','ident'),align = "v",ncol = 2) 
fig_umap
ggsave(filename = "0916整合tsne.pdf", plot = fig_umap, device = 'pdf', width = 30, height = 15, units = 'cm')


p3 <- DimPlot(testAB.integrated, reduction = "umap", group.by = "orig.ident", pt.size=0.5) 
p3
p4 <- DimPlot(testAB.integrated, reduction = "umap", group.by = "ident", pt.size=0.5, label = TRUE,repel = TRUE) 
p4
fig_umap <- plot_grid(p3, p4, labels = c('patient','ident'),align = "v",ncol = 2) 
fig_umap
ggsave(filename = "0916整合umap.pdf", plot = fig_umap, device = 'pdf', width = 30, height = 15, units = 'cm')


saveRDS(testAB.integrated,file = "0916umap.rds")
testAB.integrated <- readRDS(file = "0916umap.rds")

###区分一下常见marker,文章中的marker
celltype_marker=c("Pax9", "Neurod1","Hoxa3","Pdgfra"," Dlx5","Esam",
                  "Epcam","Rxrg","Sox10","Fstl1")
##做小提琴图
VlnPlot(testAB.integrated,features = celltype_marker,pt.size = 0,ncol = 2)

table(testAB.integrated@meta.data$seurat_clusters)
#####################################################################
#提取指定单细胞亚群重新分析
pax9_sce1 = testAB.integrated[,!(testAB.integrated@meta.data$seurat_clusters %in% c(6,26,39,40,43,45))]
pax9_sce1##53693
saveRDS(pax9_sce1,file = "0916删减umap.rds")
pax9_sce1 <- readRDS(file = "0916删减umap.rds")

cd4_sce1=pax9_sce1

cd4_sce1 <- ScaleData(cd4_sce1, features = rownames(cd4_sce1))
cd4_sce1 <- RunPCA(cd4_sce1, features = VariableFeatures(cd4_sce1),npcs = 100)
DimHeatmap(cd4_sce1, dims =1:30, cells = 500, balanced = TRUE)
ElbowPlot(cd4_sce1, ndims =100)##碎石图
cd4_sce1 <- FindNeighbors(cd4_sce1, dims = 1:60)
cd4_sce1 <- FindClusters(cd4_sce1, resolution = 0.5)
cd4_sce1 <- RunUMAP(cd4_sce1, dims = 1:60)
cd4_sce1 <- RunTSNE(cd4_sce1, dims = 1:60)
library(stringr)
library(tidyverse)
library(export)
library(cowplot)
p1 <- DimPlot(cd4_sce1, reduction = "tsne", group.by = "orig.ident", pt.size=0.5)
p2 <- DimPlot(cd4_sce1, reduction = "tsne", group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE)

fig_tsne <- plot_grid(p1, p2, labels = c('orig.ident','ident'),align = "v",ncol = 2)
fig_tsne
ggsave(filename = "0916删减_tsne.pdf", plot = fig_tsne, device = 'pdf', width = 27, height = 12, units = 'cm')

p3 <- DimPlot(cd4_sce1, reduction = "umap", group.by = "orig.ident", pt.size=0.5) 
p3
p4 <- DimPlot(cd4_sce1, reduction = "umap", group.by = "ident", pt.size=0.5, label = TRUE,repel = TRUE) 
p4
fig_umap <- plot_grid(p3, p4, labels = c('orig.ident','ident'),align = "v",ncol = 2) 
fig_umap
ggsave(filename = "0916删减_umap.pdf", plot = fig_umap, device = 'pdf', width = 30, height = 15, units = 'cm')

saveRDS(cd4_sce1,file = "0916删减-60pctest.rds")
cd4_sce1 <- readRDS(file = "0916删减-60pctest.rds")

###########################################################

####细胞类型注释
###区分一下常见marker,文章中的marker
PPE_marker <- c("Eya1","Pax9","Hoxa3","Foxn1","Bmp4","Six1","Pax1","Tbx1")
##做小提琴图
VlnPlot(cd4_sce1,features = PPE_marker,pt.size = 0,ncol = 2)
##鉴定完大类，把注释信息添加上去
table(cd4_sce1@meta.data$seurat_clusters)
Third_pharyngeal_pouch <- c(0,1,5,9,16,20)
Unknown <- c(2,3,4,6,7,8,10,11,12,13,14,15,17,18,19,21,22,23,24)
current.cluster.ids <- c(Third_pharyngeal_pouch,
                         Unknown)
new.cluster.ids <- c(rep("Third pharyngeal pouch",length(Third_pharyngeal_pouch)),
                     rep("Unknown",length(Unknown)))

cd4_sce1@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(cd4_sce1@meta.data$seurat_clusters)), 
                                          from = current.cluster.ids, to = new.cluster.ids)

##新的cd4_sce1@meta.data是这样的
##统计一下各种细胞的数目

head(cd4_sce1@meta.data)
##，画最后的tsne图

plotCB=as.data.frame(cd4_sce1@meta.data%>%filter(seurat_clusters!="13" &
                                                   seurat_clusters!="15"))[,"CB"]
colors <-  c("#FFD700","#DCDCDC")#
DimPlot(cd4_sce1, reduction = "tsne", group.by = "celltype", pt.size=0.5,cols = colors,label = TRUE,repel = TRUE)

p1 <- DimPlot(cd4_sce1, reduction = "tsne", group.by = "orig.ident", pt.size=0.5)
p2 <- DimPlot(cd4_sce1, reduction = "tsne", group.by = "celltype", pt.size=0.5,cols = colors,repel = TRUE)

fig_tsne <- plot_grid(p1, p2, labels = c('orig.ident','celltype'),align = "v",ncol = 2)
fig_tsne
ggsave(filename = "注释_tsne.pdf", plot = fig_tsne, device = 'pdf', width = 30, height = 12, units = 'cm')


saveRDS(cd4_sce1,file = "注释.rds") #保存test.seu对象，下次可以直接调用，不用重复以上步骤,cells = plotCB
cd4_sce1 <- readRDS(file ="注释.rds")

########################################################
##提取3PPE细胞群
PPE <- cd4_sce1[,(cd4_sce1@meta.data$seurat_clusters %in% c(0,1,5,9,16,20))]
table(PPE@meta.data$seurat_clusters)
Idents(PPE) <- "orig.ident"  ##细胞身份，属于哪个亚群

##挑选亚群作图
vln_df = PPE[,PPE@meta.data$orig.ident %in% c("E9_5_rep1",'E10_5_rep7',"E11_5_rep2","E12_5_rep2")]

table(vln_df@meta.data$orig.ident)
###"Eya1","Six1","Pax9","Hoxa3","Epcam"
marker=c("Eya1","Six1","Pax9","Epcam","Hoxa3")
VlnPlot(vln_df,features = marker,pt.size = 0,adjust = 1,ncol = 2)
#我们知道一个关键的参数scale = "width"导致了这种局面，其他没有出现小提琴的应该是零值比例太多
#数据过滤。
VlnPlot(subset(vln_df,Hoxa3 > 0.2 ), "Hoxa3",adjust = .5,pt.size=0)+ theme_bw()+NoLegend()
#改腰围adjust = .5

p1 <- VlnPlot(subset(vln_df,Hoxa3 > 0.2 ),"Hoxa3",adjust = .8,pt.size=0)+ theme_bw()+geom_boxplot(width=.2,col="black",fill="white")+  NoLegend()
p2 <- VlnPlot(subset(vln_df,Eya1 > 0.2 ),"Eya1",adjust = .8,pt.size=0)+ theme_bw()+geom_boxplot(width=.2,col="black",fill="white")+  NoLegend()
p3 <- VlnPlot(subset(vln_df,Six1 > 0.2 ),"Six1",adjust = .8,pt.size=0)+ theme_bw()+geom_boxplot(width=.2,col="black",fill="white")+  NoLegend()
p4 <- VlnPlot(subset(vln_df,Epcam > 0.2 ),"Epcam",adjust = .8,pt.size=0)+ theme_bw()+geom_boxplot(width=.2,col="black",fill="white") +NoLegend()
p5 <- VlnPlot(subset(vln_df,Pax9 > 0.2 ),"Pax9",adjust = .8,pt.size=0)+ theme_bw()+geom_boxplot(width=.2,col="black",fill="white") +NoLegend()

library(patchwork)##拼图
#p1+p2+p3+p4+p5+plot_layout(nrow=5)
(p1|p2)/(p3|p4)/(p5|plot_spacer())
ggsave(filename = "3PPE小提琴.pdf", width = 22, height = 20, 
       device = 'pdf', units = 'cm')





#批量画图
gene <- c("Epcam","Pax9","Bmp4","Hoxa3","Pax1","Ptprc","Six1","Tbx1") #Figure1_A和Figure1_D
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

#Figure1_E批量画图
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




```

