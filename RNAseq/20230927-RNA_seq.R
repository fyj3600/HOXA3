#RNA-seq data analysis of hESCs-derived third pharyngeal capsule endodermal cells #before and after HOXA3 inhibitor treatment
#####################################################################
#Compiled: September 27, 2023
#The provider：Xueyan Zhang

#step1. calculate FPKM
#Get gene length
library(dplyr)
install.packages("BiocManager")
BiocManager::install("GenomicFeatures")
library(GenomicFeatures)
#Human annotation files (Homo_sapiens.GRCh38.98.chr.gtf) were downloaded from the NCBI database.
#Open the Human annotation file and process it
txdb <- makeTxDbFromGFF("Homo_sapiens.GRCh38.98.chr.gtf", format="gtf")
exons_gene <- exonsBy(txdb,by="gene")
exons_gene_lens <- lapply(exons_gene, function(x){sum(width(reduce(x)))})
#convert an object into a data frame
exons_gene_lens1 <- as.data.frame(exons_gene_lens)
class(exons_gene_lens1)
dim(exons_gene_lens1)
#Swap of rows and columns
exons_gene_lens1 <- t(exons_gene_lens1)
dim(exons_gene_lens1)
#save file
write.csv(exons_gene_lens1, file="gene_Length.csv")
#Re-read the file
exons_gene_lens2 <- read.csv("gene_Length.csv",header=T)
colnames(exons_gene_lens2) <- c("Geneid","Length")
##Read our sequencing data. 
readscounts <- read.table("counts_matrix.txt", header=T)
###Merge gene length files and our sequencing data.
mycounts <- merge(exons_gene_lens2, readscounts, by="Geneid", all=F)
head(mycounts)
##Change a line name.
rownames(mycounts)<-mycounts[,1]
mycounts<-mycounts[,-1]
#calculate TPM and FPKM
#TPM
kb <- mycounts$Length / 1000
countdata <- mycounts[,2:7]
rpk <- countdata / kb
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
head(tpm)
write.csv(tpm, file="mRNA_tpm.csv")
##FPKM
fpkm <- t(t(rpk)/colSums(countdata) * 10^6) 
head(fpkm)
write.csv(fpkm, file="mRNA_fpkm.csv")

###############################################################
#step2. Gene ID is converted to gene name
if(!requireNamespace("BiocManager",quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")
library(biomaRt)
#Select a database
mymart<-useMart("ENSEMBL_MART_ENSEMBL")
#The selected Human Genome (GRCh38.p13)
txdb <- useDataset("hsapiens_gene_ensembl",mart=mymart)
#Open our own data and determine the data type.
mycounts<-read.table("counts_matrix.txt",header=TRUE)
dim(mycounts)
#Take out the gene name
my_names<-mycounts[,1]
#Extract the desired annotation type
my_symbols<-getBM(attributes=c("external_gene_name","ensembl_gene_id",
                               "hgnc_symbol"),
                  filters="ensembl_gene_id",
                  values=my_names,
                  mart=txdb)
head(my_symbols)
colnames(my_symbols)<-c('Gene_name','Geneid','symbol')
head(mycounts)
mycounts <- merge(my_symbols, mycounts, by="Geneid", all=F)
mycounts<-mycounts[,-c(1,2)]
###save file
rows <- rownames(unique(mycounts['symbol']))
mycounts <- mycounts[rows,]
write.table(mycounts,"counts_matrix_symbol.txt",row.names=FALSE)

##################################################################
#step3. analysis of difference (DESeq2).
library(tidyverse)
library(DESeq2)
library(stringr)
#Read our sequencing data. 
mycounts<-read.table("counts_matrix_symbol.txt", header=TRUE)
head(mycounts)
rownames(mycounts)<-mycounts[,1]
mycounts <- mycounts[,-1]
head(mycounts)
#Create a group
condition = rep(c("control","treatment"),times = c(3,3))
sample_name = colnames(mycounts)
colData <- data.frame(sample_name,condition)
head(colData)
###Construct matrix
dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
dds <- DESeq(dds)
dds
sizeFactors(dds)
res <- results(dds, contrast=c("condition", "treatment","control"))
res <- res[order(res$pvalue),]
head(res)
summary(res)
class(res)
#The differential gene files were saved
write.csv(res,file="all_DESeq2_results.csv")
diff_gene_deseq2 <-subset(res, pvalue < 0.05 & abs(log2FoldChange) > 1)
dim(diff_gene_deseq2)  
head(diff_gene_deseq2)
diff_gene_deseq2 <- na.omit(diff_gene_deseq2)
write.csv(diff_gene_deseq2,file= "DESeq2_diffExpression.csv")

####################################################################
#step4. Heat map of differentially expressed genes.
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(export)
#Read our sequencing data. 
dat <-read.table("counts_matrix_symbol.txt",header=TRUE)
rownames(dat)<-dat[,1]
dat <- dat[,-1]
head(dat)
#Figure3E
F3_E <-dat[c('BMP4','TBX1','PBX1','HOXA3','ETV5'),]
#Create a group
condition = rep(c("control","treatment"),times = c(3,3))
sample_name = colnames(F3_E)
colData <- data.frame(sample_name,condition)
head(colData)
dat1 <- t(scale(t(F3_E)))
group <- data.frame(colData$condition)
rownames(group)=colnames(dat1)
g1 <- pheatmap(dat1,show_colnames =T, show_rownames = T, cluster_cols = F,
               border=FALSE,
               cellwidth = 40,
               cellheight = 40,
               color = colorRampPalette(c("navy", "white", "firebrick3"))(50))	 
graph2ppt(file="Figure3E.ppt", width=10, aspectr=1)
#Figure3O
F3_O <-dat[c('PTPRT','NRTN','TLX3','INSL3','CTHRC1','ANGPT1'),]
dat2 <- t(scale(t(F3_O)))
group <- data.frame(colData$condition)
rownames(group)=colnames(dat2)
g2 <- pheatmap(dat2,show_colnames =T, show_rownames = T, cluster_cols = F,
               border=FALSE,
               cellwidth = 40,
               cellheight = 40,
               color = colorRampPalette(c("navy", "white", "firebrick3"))(50))	 
graph2ppt(file="Figure3O.ppt", width=10, aspectr=1)
#Figure2B
F2_B <-dat[c('WNT7B','WNT4','WNT9A','WNT7A','WNT8B','LRP5',
             'LRP5','DVL3','FRAT1','GSK3B','AXIN1',
             'APC','CTNNB1','TCF4','CCND1'),]
dat3 <- t(scale(t(F2_B)))
group <- data.frame(colData$condition)
rownames(group)=colnames(dat3)
g3 <- pheatmap(dat3,show_colnames =T, show_rownames = T, cluster_cols = F,
               border_color = "white",
               cellwidth = 40,
               cellheight = 20,
               color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
graph2ppt(file="Figure2B.ppt", width=10, aspectr=1)

