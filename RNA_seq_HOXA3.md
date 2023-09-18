

```R
#HOXA3抑制剂处理hESCs来源的第三咽囊内胚层细胞前后RNA-seq数据分析

#######################################################################
#step1. 计算FPKM
rm(list = ls())
#1. 获取基因长度
library(dplyr)
install.packages("BiocManager")
BiocManager::install("GenomicFeatures")
library(GenomicFeatures)
#NCBI数据库中下载人类注释文件Homo_sapiens.GRCh38.98.chr.gtf到本地。
#打开注释文件并整理
txdb <- makeTxDbFromGFF("Homo_sapiens.GRCh38.98.chr.gtf", format="gtf")
exons_gene <- exonsBy(txdb,by="gene")
exons_gene_lens <- lapply(exons_gene, function(x){sum(width(reduce(x)))})
#查看
class(exons_gene_lens)
length(exons_gene_lens)
#转换成data frame
exons_gene_lens1 <- as.data.frame(exons_gene_lens)
class(exons_gene_lens1)
dim(exons_gene_lens1)
#转置
exons_gene_lens1 <- t(exons_gene_lens1)
dim(exons_gene_lens1)
#将文件保存
write.csv(exons_gene_lens1, file="gene_Length.csv")
#重新读入
exons_gene_lens2 <- read.csv("gene_Length.csv",header=T)
colnames(exons_gene_lens2) <- c("Geneid","Length")
##载入reads counts的原始数据。
readscounts <- read.table("counts_matrix.txt", header=T)
###将exons_gene_lens1和readscounts文件合并
mycounts <- merge(exons_gene_lens2, readscounts, by="Geneid", all=F)
dim(mycounts)
mycounts[1:3,1:3]#查看
##整理数据。
rownames(mycounts)<-mycounts[,1]
mycounts<-mycounts[,-1]
mycounts[1:3,1:3]

#2.计算TPM,FPKM
#计算TPM
kb <- mycounts$Length / 1000
countdata <- mycounts[,2:7]#第一列Length不需要了
rpk <- countdata / kb##计算每个值的rpk
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
head(tpm)
write.csv(tpm, file="mRNA_tpm.csv")
##计算FPKM
fpkm <- t(t(rpk)/colSums(countdata) * 10^6) 
head(fpkm)
write.csv(fpkm, file="mRNA_fpkm.csv")
#####################################################################

#step2. 基因ID转换

if(!requireNamespace("BiocManager",quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")
library(biomaRt)
#选择需要的数据库
mymart<-useMart("ENSEMBL_MART_ENSEMBL")
#选定人类基因组Human genes (GRCh38.p13)
txdb <- useDataset("hsapiens_gene_ensembl",mart=mymart)
#打开自己的数据，确定数据类型。
mycounts<-read.table("counts_matrix.txt",header=TRUE)
dim(mycounts)
#取出基因名
my_names<-mycounts[,1]
#提取想要的注释类型
my_symbols<-getBM(attributes=c("external_gene_name","ensembl_gene_id",
                               "hgnc_symbol"),#我们要获取的注释类型
                  filters="ensembl_gene_id",#我们已知的注释类型
                  values=my_names,#已知的注释类型的数据
                  mart=txdb)#选定的数据库的基因组
head(my_symbols)
#与我们的数据合并
colnames(my_symbols)<-c('Gene_name','Geneid','symbol')
head(mycounts)
mycounts <- merge(my_symbols, mycounts, by="Geneid", all=F)
mycounts<-mycounts[,-c(1,2)]
###保存数据
rows <- rownames(unique(mycounts['symbol']))#去除重复
mycounts <- mycounts[rows,]
write.table(mycounts,"counts_matrix_symbol.txt",row.names=FALSE)

#step3. DESeq2做差异分析

rm(list = ls()) #清除环境内存
library(tidyverse)
library(DESeq2)

mycounts<-read.table("counts_matrix_symbol.txt", header=TRUE)
head(mycounts)
#更改行名
rownames(mycounts)<-mycounts[,1]
mycounts <- mycounts[,-1]
head(mycounts)

library(stringr)
#创建二分组数据框
condition = rep(c("control","treatment"),times = c(3,3))
sample_name = colnames(mycounts)
colData <- data.frame(sample_name,condition)
head(colData)

###构建dds矩阵
dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
##对dds进行标准化
dds <- DESeq(dds)
dds
sizeFactors(dds)
##提出结果"treatment"vs"control"
res <- results(dds, contrast=c("condition", "treatment","control"))
res <- res[order(res$pvalue),]
head(res)
summary(res)
class(res)
write.csv(res,file="all_DESeq2_results.csv")
table(res$pvalue<0.05)#统计
#筛选想要的列
diff_gene_deseq2 <-subset(res, pvalue < 0.05 & abs(log2FoldChange) > 1)
dim(diff_gene_deseq2)  
head(diff_gene_deseq2)
diff_gene_deseq2 <- na.omit(diff_gene_deseq2)#返回不含NA的值
write.csv(diff_gene_deseq2,file= "DESeq2_diffExpression.csv")

#step4. 差异基因热图
rm(list = ls()) #清除环境内存
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(export)
# 加载原始表达矩阵
dat <-read.table("counts_matrix_symbol.txt",header=TRUE)
##修改行名
rownames(dat)<-dat[,1]
dat <- dat[,-1]
head(dat)

#提取出Figure3E热图所需基因
F3_E <-dat[c('BMP4','TBX1','PBX1','HOXA3','ETV5'),]
#创建分组数据
condition = rep(c("control","treatment"),times = c(3,3))
sample_name = colnames(F3_E)
colData <- data.frame(sample_name,condition)
head(colData)

dat1 <- t(scale(t(F3_E)))

##正常情况下，需要修改这个名字。
group <- data.frame(colData$condition)
rownames(group)=colnames(dat1)
##绘制热图
g2 <- pheatmap(dat1,show_colnames =T, show_rownames = T, cluster_cols = F,
               border=FALSE,
               cellwidth = 40,
               cellheight = 40,
               color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
g2		 
graph2ppt(file="Figure3E.ppt", width=10, aspectr=1)

#提取出Figure3O热图所需基因
F3_O <-dat[c('PTPRT','NRTN','TLX3','INSL3','CTHRC1','ANGPT1'),]
dat2 <- t(scale(t(F3_O)))
##正常情况下，需要修改这个名字。
group <- data.frame(colData$condition)
rownames(group)=colnames(dat2)
##绘制热图
g2 <- pheatmap(dat2,show_colnames =T, show_rownames = T, cluster_cols = F,
               border=FALSE,
               cellwidth = 40,
               cellheight = 40,
               color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
g2		 
graph2ppt(file="Figure3O.ppt", width=10, aspectr=1)

#提取出Figure2B热图所需基因
F2_B <-dat[c('WNT7B','WNT4','WNT9A','WNT7A','WNT8B','LRP5','LRP5','DVL3','FRAT1','GSK3B','AXIN1',
             'APC','CTNNB1','TCF4','CCND1'),]
dat3 <- t(scale(t(F2_B)))
##正常情况下，需要修改这个名字。
group <- data.frame(colData$condition)
rownames(group)=colnames(dat3)
##绘制热图
g2 <- pheatmap(dat3,show_colnames =T, show_rownames = T, cluster_cols = F,
               border_color = "white",
               cellwidth = 40,
               cellheight = 20,
               color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
g2		 
graph2ppt(file="Figure2B.ppt", width=10, aspectr=1)




```