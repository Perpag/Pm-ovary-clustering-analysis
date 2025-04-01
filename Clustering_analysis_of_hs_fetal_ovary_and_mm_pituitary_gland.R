library(pheatmap)
library(tidyverse)
library(cowplot)
library(dplyr)
library(Seurat)

#human fetal ovary dataaset
#hs_7
hs_7.data <- Read10X(data.dir = "./hs_7")
hs_7_barcodes <- read.delim(file = "./hs_7/barcodes.tsv", stringsAsFactors = F, header = F)
hs_7 <- CreateSeuratObject(hs_7.data, project = "hs_7", min.cells = 3, min.features = , meta.data = rownames(hs_7_barcodes$v1)) #check values for cells and features
VlnPlot(hs_7, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
hs_7<- subset(hs_7, subset = nFeature_RNA > 2500 & nFeature_RNA < 5000)

#hs_9
hs_9.data <- Read10X(data.dir = "./hs_9")
hs_9_barcodes <- read.delim(file = "./hs_9/barcodes.tsv", stringsAsFactors = F, header = F)
hs_9 <- CreateSeuratObject(hs_9.data, project = "hs_9", min.cells = 3, min.features = , meta.data = rownames(hs_9_barcodes$v1)) #check values for cells and features
VlnPlot(hs_9, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
hs_9<- subset(hs_9, subset = nFeature_RNA > 1500 & nFeature_RNA < 5000)

#hs_10
hs_10.data <- Read10X(data.dir = "./hs_10")
hs_10_barcodes <- read.delim(file = "./hs_10/barcodes.tsv", stringsAsFactors = F, header = F)
hs_10 <- CreateSeuratObject(hs_10.data, project = "hs_10", min.cells = 3, min.features = , meta.data = rownames(hs_10_barcodes$v1)) #check values for cells and features
VlnPlot(hs_10, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
hs_10<- subset(hs_10, subset = nFeature_RNA > 2000 & nFeature_RNA < 6000)

#hs_13
hs_13.data <- Read10X(data.dir = "./hs_13")
hs_13_barcodes <- read.delim(file = "./hs_13/barcodes.tsv", stringsAsFactors = F, header = F)
hs_13 <- CreateSeuratObject(hs_13.data, project = "hs_13", min.cells = 3, min.features = , meta.data = rownames(hs_13_barcodes$v1)) #check values for cells and features
VlnPlot(hs_13, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
hs_13<- subset(hs_13, subset = nFeature_RNA > 2000 & nFeature_RNA < 6000)

#hs_16
hs_16.data <- Read10X(data.dir = "./hs_16")
hs_16_barcodes <- read.delim(file = "./hs_16/barcodes.tsv", stringsAsFactors = F, header = F)
hs_16 <- CreateSeuratObject(hs_16.data, project = "hs_16", min.cells = 3, min.features = , meta.data = rownames(hs_16_barcodes$v1)) #check values for cells and features
VlnPlot(hs_16, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
hs_16<- subset(hs_16, subset = nFeature_RNA > 2000 & nFeature_RNA < 5000)

hs_7 <- RenameCells(hs_7, add.cell.id = "hs_7")
hs_9 <- RenameCells(hs_9, add.cell.id = "hs_9")
hs_10 <- RenameCells(hs_10, add.cell.id = "hs_10")
hs_13 <- RenameCells(hs_13, add.cell.id = "hs_13")
hs_16 <- RenameCells(hs_16, add.cell.id = "hs_16")


hs_ovary<-merge(hs_7, c(hs_9, hs_10, hs_13, hs_16))


hs_ovary[["percent.mt"]] <- PercentageFeatureSet(hs_ovary, pattern = "^MT-")
hs_ovary <- subset(hs_ovary, subset = percent.mt < 5)

hs_ovary <- NormalizeData(hs_ovary)
hs_ovary <- FindVariableFeatures(hs_ovary, selection.method = "vst", nfeatures = 2000)
hs_ovary<- ScaleData(object = hs_ovary)
hs_ovary<- RunPCA(object = hs_ovary, features = VariableFeatures(object = hs_ovary))
hs_ovary_jackstraw<-JackStraw(hs_ovary, reduction = "pca", assay = NULL, dims = 50,
                  num.replicate = 100, prop.freq = 0.01, verbose = TRUE,
                  maxit = 1000)
ScoreJackStraw(hs_ovary_jackstraw, reduction = "pca", dims = 1:50,
               score.thresh = 1e-05, do.plot = TRUE)
hs_ovary<- FindNeighbors(object = hs_ovary, dims = 1:50)
hs_ovary<- FindClusters(object = hs_ovary, resolution = 1)
hs_ovary <- RunUMAP(object = hs_ovary, dims = 1:50)
hs_ovary[["RNA"]] <- JoinLayers(hs_ovarys[["RNA"]])

hs_ovary<-subset(hs_ovary, downsample=2000)


hs_ovary<- RenameIdents(object = hs_ovary, `33` = "34 ")
hs_ovary<- RenameIdents(object = hs_ovary, `32` = "33 ")
hs_ovary<- RenameIdents(object = hs_ovary, `31` = "32 ")
hs_ovary<- RenameIdents(object = hs_ovary, `30` = "31 ")
hs_ovary<- RenameIdents(object = hs_ovary, `29` = "30 ")
hs_ovary<- RenameIdents(object = hs_ovary, `28` = "29 ")
hs_ovary<- RenameIdents(object = hs_ovary, `27` = "28 ")
hs_ovary<- RenameIdents(object = hs_ovary, `26` = "27 ")
hs_ovary<- RenameIdents(object = hs_ovary, `25` = "26 ")
hs_ovary<- RenameIdents(object = hs_ovary, `24` = "25 ")
hs_ovary<- RenameIdents(object = hs_ovary, `23` = "24 ")
hs_ovary<- RenameIdents(object = hs_ovary, `22` = "23 ")
hs_ovary<- RenameIdents(object = hs_ovary, `21` = "22 ")
hs_ovary<- RenameIdents(object = hs_ovary, `20` = "21 ")
hs_ovary<- RenameIdents(object = hs_ovary, `19` = "20 ")
hs_ovary<- RenameIdents(object = hs_ovary, `18` = "19 ")
hs_ovary<- RenameIdents(object = hs_ovary, `17` = "18 ")
hs_ovary<- RenameIdents(object = hs_ovary, `16` = "17 ")
hs_ovary<- RenameIdents(object = hs_ovary, `15` = "16 ")
hs_ovary<- RenameIdents(object = hs_ovary, `14` = "15 ")
hs_ovary<- RenameIdents(object = hs_ovary, `13` = "14 ")
hs_ovary<- RenameIdents(object = hs_ovary, `12` = "13 ")
hs_ovary<- RenameIdents(object = hs_ovary, `11` = "12 ")
hs_ovary<- RenameIdents(object = hs_ovary, `10` = "11 ")
hs_ovary<- RenameIdents(object = hs_ovary, `9` = "10 ")
hs_ovary<- RenameIdents(object = hs_ovary, `8` = "9 ")
hs_ovary<- RenameIdents(object = hs_ovary, `7` = "8 ")
hs_ovary<- RenameIdents(object = hs_ovary, `6` = "7 ")
hs_ovary<- RenameIdents(object = hs_ovary, `5` = "6 ")
hs_ovary<- RenameIdents(object = hs_ovary, `4` = "5 ")
hs_ovary<- RenameIdents(object = hs_ovary, `3` = "4 ")
hs_ovary<- RenameIdents(object = hs_ovary, `2` = "3 ")
hs_ovary<- RenameIdents(object = hs_ovary, `1` = "2 ")
hs_ovary<- RenameIdents(object = hs_ovary, `0` = "1 ")

#mouse pituitary gland dataset
#mm_1
mm_1.data <- Read10X(data.dir = "./mm_1")
mm_1_barcodes <- read.delim(file = "./mm_1/barcodes.tsv", stringsAsFactors = F, header = F)
mm_1 <- CreateSeuratObject(mm_1.data, project = "mm_1", min.cells = 3, min.features = , meta.data = rownames(mm_1_barcodes$v1)) #check values for cells and features
VlnPlot(mm_1, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
mm_1<- subset(mm_1, subset = nFeature_RNA > 1500 & nFeature_RNA < 5000)

#mm_2
mm_2.data <- Read10X(data.dir = "./mm_2")
mm_2_barcodes <- read.delim(file = "./mm_2/barcodes.tsv", stringsAsFactors = F, header = F)
mm_2 <- CreateSeuratObject(mm_2.data, project = "mm_2", min.cells = 3, min.features = , meta.data = rownames(mm_2_barcodes$v1)) #check values for cells and features
VlnPlot(mm_2, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
mm_2<- subset(mm_2, subset = nFeature_RNA > 1500 & nFeature_RNA < 5000)
mm_1[["percent.mt"]] <- PercentageFeatureSet(mm_1, pattern = "^MT-")
mm_1 <- subset(mm_1, subset = percent.mt < 5)

mm_1<- RenameCells(mm_1, add.cell.id = "mm_1")
mm_2<- RenameCells(mm_2, add.cell.id = "mm_2")

pg <- list(mm_1, mm_2)

for (i in 1:length(pd)) {
  
  pd[[i]] <- NormalizeData(pg[[i]], 
                           verbose = FALSE)
  pg[[i]] <- FindVariableFeatures(pg[[i]], 
                                  selection.method = "vst",
                                  nfeatures = 2000,
                                  verbose = FALSE)
  
}
pg_anchors <- FindIntegrationAnchors(object.list = pg,
                                     dims = 1:50)

pg_integrated <- IntegrateData(anchorset = pg_anchors, dims = 1:50)

#Clustering                                                                                         
pg_integrated <- ScaleData(object = pg_integrated)
pg_integrated <- RunPCA(object = pg_integrated, features = VariableFeatures(object = pg_integrated))
pg_jackstraw<-JackStraw(pg_integrated, reduction = "pca", assay = NULL, dims = 50,
                 num.replicate = 100, prop.freq = 0.01, verbose = TRUE,
                 maxit = 1000)
ScoreJackStraw(pg_jackstraw, reduction = "pca", dims = 1:50,
               score.thresh = 1e-05, do.plot = TRUE)

pg_integrated <- FindNeighbors(object = pg_integrated, dims = 1:50)
pg_integrated<- FindClusters(object = pg_integrated, resolution = 1)

#UMAP
pg_integrated <- RunUMAP(object = pg_integrated, dims = 1:50)
DimPlot(pg_integrated, group.by="orig.ident")
DimPlot(object = pg_integrated, reduction = "umap", label = TRUE)
DimPlot(pg_integrated, group.by="orig.ident", split.by = "orig.ident")
pg_integrated[["RNA"]] <- JoinLayers(pg_integrated[["RNA"]])
