library(pheatmap)
library(tidyverse)
library(cowplot)
library(dplyr)
library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(patchwork)
library(EnhancedVolcano)
set.seed(312)

#loading the data
#1
scf1.data <- Read10X(data.dir = "./scf1")
scf1_barcodes <- read.delim(file = "./scf1/barcodes.tsv", stringsAsFactors = F, header = F)
scf1 <- CreateSeuratObject(scf1.data, project = "scf1", min.cells = 3, min.features = , meta.data = rownames(scf1_barcodes$v1)) #check values for cells and features
VlnPlot(scf1, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
scf1<- subset(scf1, subset = nFeature_RNA > 400 & nFeature_RNA < 5000)

#2
scf5.data <- Read10X(data.dir = "./scf5")
scf5_barcodes <- read.delim(file = "./scf1/barcodes.tsv", stringsAsFactors = F, header = F)
scf5 <- CreateSeuratObject(scf5.data, project = "scf5", min.cells = 3, min.features = , meta.data = rownames(scf5_barcodes$v1)) #check values for cells and features
VlnPlot(scf5, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
scf5<- subset(scf5, subset = nFeature_RNA > 400 & nFeature_RNA < 5000)

#3
scf9.data <- Read10X(data.dir = "./scf9")
scf9_barcodes <- read.delim(file = "./scf9/barcodes.tsv", stringsAsFactors = F, header = F)
scf9 <- CreateSeuratObject(scf9.data, project = "scf9", min.cells = 3, min.features = , meta.data = rownames(scf9_barcodes$v1)) #check values for cells and features
VlnPlot(scf9, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
scf9<- subset(scf9, subset = nFeature_RNA > 400 & nFeature_RNA < 5000)

#Renamimg
scf1 <- RenameCells(scf1, add.cell.id = "scf1")
scf5 <- RenameCells(scf5, add.cell.id = "scf5")
scf9 <- RenameCells(scf9, add.cell.id = "scf9")

pm_ovary_list <- list(scf1, scf5, scf9)

for (i in 1:length(pm_ovary_list)) {
  
  pm_ovary_list[[i]] <- NormalizeData(pm_ovary_list[[i]], 
                                      verbose = FALSE)
  pm_ovary_list[[i]] <- FindVariableFeatures(pm_ovary_list[[i]], 
                                             selection.method = "vst",
                                             nfeatures = 2000,
                                             verbose = FALSE)
  
}
pm_ovary_anchors <- FindIntegrationAnchors(object.list = pm_ovary_list,
                                           dims = 1:50)

pm_ovary <- IntegrateData(anchorset = pm_ovary_anchors, dims = 1:50)

#Clustering
pm_ovary_integrated <- ScaleData(object = pm_ovary)
pm_ovary_integrated <- RunPCA(object = pm_ovary_integrated, features = VariableFeatures(object = pm_ovary_integrated))
pm_ovary_jackstraw<-JackStraw(pm_ovary_integrated, reduction = "pca", assay = NULL, dims = 50,
               num.replicate = 100, prop.freq = 0.01, verbose = TRUE,
               maxit = 1000)
ScoreJackStraw(pm_ovary_jackstraw, reduction = "pca", dims = 1:50,
               score.thresh = 1e-05, do.plot = TRUE)
pm_ovary_integrated <- FindNeighbors(object = pm_ovary_integrated, dims = 1:50)
pm_ovary_integrated<- FindClusters(object = pm_ovary_integrated, resolution = 1)

#UMAP
pm_ovary_integrated_umap <- RunUMAP(object = pm_ovary_integrated, dims = 1:50)
DimPlot(pm_ovary_integrated_umap, group.by="orig.ident")
DimPlot(object = pm_ovary_integrated_umap, reduction = "umap", label = TRUE)
DimPlot(pm_ovary_integrated_umap, group.by="orig.ident", split.by = "orig.ident")

pm_ovary_integrated_umap<- RenameIdents(object = pm_ovary_integrated_umap, `21` = "27")
pm_ovary_integrated_umap<- RenameIdents(object = pm_ovary_integrated_umap, `24` = "26")
pm_ovary_integrated_umap<- RenameIdents(object = pm_ovary_integrated_umap, `16` = "25")
pm_ovary_integrated_umap<- RenameIdents(object = pm_ovary_integrated_umap, `2` = "24")
pm_ovary_integrated_umap<- RenameIdents(object = pm_ovary_integrated_umap, `10` = "23")
pm_ovary_integrated_umap<- RenameIdents(object = pm_ovary_integrated_umap, `20` = "22")
pm_ovary_integrated_umap<- RenameIdents(object = pm_ovary_integrated_umap, `18` = "21")
pm_ovary_integrated_umap<- RenameIdents(object = pm_ovary_integrated_umap, `22` = "20")
pm_ovary_integrated_umap<- RenameIdents(object = pm_ovary_integrated_umap, `2` = "19")
pm_ovary_integrated_umap<- RenameIdents(object = pm_ovary_integrated_umap, `23` = "18")
pm_ovary_integrated_umap<- RenameIdents(object = pm_ovary_integrated_umap, `7` = "17")
pm_ovary_integrated_umap<- RenameIdents(object = pm_ovary_integrated_umap, `9` = "16")
pm_ovary_integrated_umap<- RenameIdents(object = pm_ovary_integrated_umap, `15` = "15")
pm_ovary_integrated_umap<- RenameIdents(object = pm_ovary_integrated_umap, `25` = "14")
pm_ovary_integrated_umap<- RenameIdents(object = pm_ovary_integrated_umap, `14` = "13")
pm_ovary_integrated_umap<- RenameIdents(object = pm_ovary_integrated_umap, `11` = "12")
pm_ovary_integrated_umap<- RenameIdents(object = pm_ovary_integrated_umap, `17` = "11")
pm_ovary_integrated_umap<- RenameIdents(object = pm_ovary_integrated_umap, `3` = "10")
pm_ovary_integrated_umap<- RenameIdents(object = pm_ovary_integrated_umap, `0` = "9")
pm_ovary_integrated_umap<- RenameIdents(object = pm_ovary_integrated_umap, `13` = "8")
pm_ovary_integrated_umap<- RenameIdents(object = pm_ovary_integrated_umap, `1` = "7")
pm_ovary_integrated_umap<- RenameIdents(object = pm_ovary_integrated_umap, `5` = "6")
pm_ovary_integrated_umap<- RenameIdents(object = pm_ovary_integrated_umap, `12` = "5")
pm_ovary_integrated_umap<- RenameIdents(object = pm_ovary_integrated_umap, `19` = "4")
pm_ovary_integrated_umap<- RenameIdents(object = pm_ovary_integrated_umap, `26` = "3")
pm_ovary_integrated_umap<- RenameIdents(object = pm_ovary_integrated_umap, `8` = "2")
pm_ovary_integrated_umap<- RenameIdents(object = pm_ovary_integrated_umap, `6` = "1")

markers <- FindAllMarkers(object = pm_ovary_integrated_umap, assay = "RNA", only.pos = TRUE, min.pct = 0.01 )

#cell type groups
epithelium<- WhichCells(pm_ovary_integrated_umap, idents = c("1", "2", "3", "4"))
follicle1<- WhichCells(pm_ovary_integrated_umap, idents = c("5", "6", "7", "8"))
follicle2<- WhichCells(pm_ovary_integrated_umap, idents = c("9", "10", "11", "12"))
germ<- WhichCells(pm_ovary_integrated_umap, idents = c("13", "14"))
oocytes<- WhichCells(pm_ovary_integrated_umap, idents = c("15", "16", "17", "18", "19"))
neurons<- WhichCells(pm_ovary_integrated_umap, idents = "20")
muscles<- WhichCells(pm_ovary_integrated_umap, idents = c("21", "22"))
follicle3<- WhichCells(pm_ovary_integrated_umap, idents = c("23", "24", "25"))
immune<- WhichCells(pm_ovary_integrated_umap, idents = c("26", "27"))

DimPlot(pm_ovary_integrated_umap, label=TRUE,reduction = "umap", cells.highlight= list(epithelium, neurons,germ,oocytes,muscles, follicle3, follicle1, follicle2, immune), cols.highlight = c("orange", "salmon","pink","#be80ff6d","#00bf7661","green", "yellow", "#3da1ff78"), cols= "lightblue", sizes.highlight = 0.1)     


# trajectory analysis with monocle3
cds<-subset(pm_ovary_integrated_umap, idents=c("13", "14","15", "16", "17", "18", "19"))
cds <- as.cell_data_set(cds)
cds <- cluster_cells(cds, resolution=1e-3)
p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)
cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)
cds <- order_cells(cds)
plot_cells(cds,
color_cells_by = "pseudotime",label_cell_groups=FALSE,label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)

EnhancedVolcano(Supplementary_file_4, lab = Supplementary_file_4$Description, x = "log2FC", y = "pvalue", xlab = bquote(~log[2]~"average fold chage"), ylab = bquote(~Log[10]~italic(pval-adj)), pCutoff = 5e-2, FCcutoff = 1, legendLabels = c("NS", expression(Log[2]~FC, "pval-adj", expression(pval-adj~and~log[2]~FC))), legendPosition = "bottom")


