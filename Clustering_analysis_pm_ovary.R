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


#monocle3
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


#gene_ontology_enrichment
# Load GAF file (replace with actual file path)
gaf_data <- read.delim("GCF_015706575.1_ASM1570657v1_gene_ontology.gaf", 
                       header = FALSE, sep = "\t", comment.char = "!")


# Extract relevant columns (GeneID and GO Term)
gene2go <- gaf_data[, c(3, 5)]  # Column 2 = Gene ID, Column 5 = GO term
colnames(gene2go) <- c("GeneID", "GO_Term")

# Convert to list format required by topGO
gene2go_list <- split(gene2go$GO_Term, gene2go$GeneID)

# Separate upregulated and downregulated genes
upregulated_genes <- Supplementary_file_4[Supplementary_file_4$log2FC > 0, ]
downregulated_genes <- Supplementary_file_4[Supplementary_file_4$log2FC < 0, ]

# Define Upregulated (log2FC > 0) and Downregulated (log2FC < 0)
Supplementary_file_4$Regulation <- ifelse(deg_data$log2FC > 0, "Upregulated", "Downregulated")

# Create a gene factor list for topGO (1 = DEG, 0 = non-DEG)
geneList <- factor(as.integer(gene2go$GeneID %in% Supplementary_file_4$GeneID))
names(geneList) <- gene2go$GeneID

# Check data
head(Supplementary_file_4)

# Create topGO data object
GOdata <- new("topGOdata", 
              ontology = "BP",  # Biological Process; can also use "MF" or "CC"
              allGenes = geneList,
              annot = annFUN.gene2GO, 
              gene2GO = gene2go_list)

# Run Fisher's test for enrichment
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

# Get top enriched GO terms
go_results <- GenTable(GOdata, classicFisher = resultFisher, topNodes = 10)

# View results
print(go_results)

# Extract significant GO terms
sig_go_terms <- go_results$GO.ID

# Extract genes associated with these GO terms
significant_genes <- unique(unlist(lapply(sig_go_terms, function(go) {
  genesInTerm(GOdata, go)
})))

gene2go <- gene2go[, c(2, 1)]  # Ensure GO_Term is first, GeneID second
colnames(gene2go) <- c("GO_Term", "GeneID")  # Rename columns

# Filter only genes in significant GO terms
upregulated_filtered <- upregulated_genes[upregulated_genes$GeneID %in% significant_genes, ]

# Perform enrichment for upregulated genes
ego_up <- enricher(
  gene = as.character(upregulated_filtered$GeneID), 
  TERM2GENE = gene2go,  
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05
)

# Plot enrichment for upregulated genes
dotplot(ego_up, showCategory = 20) + labs(title = "GO Enrichment for Upregulated Genes")

# Filter only genes in significant GO terms
downregulated_filtered <- downregulated_genes[downregulated_genes$GeneID %in% significant_genes, ]

# Perform enrichment for downregulated genes
ego_down <- enricher(
  gene = as.character(downregulated_filtered$GeneID), 
  TERM2GENE = gene2go,  
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05
)

# Plot enrichment for downregulated genes
dotplot(ego_down, showCategory = 20) + labs(title = "GO Enrichment for Downregulated Genes")


#Cross species comparison-SAMap

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


library(ggalluvial)
library(networkD3)
library(dplyr)
library(ggplot2)
library(reshape2)
library(Seurat)
library(SeuratData)
library(SeuratDisk)

pm<-pm_ovary_integrated_umap
sp<-sp_ovary
spl<-sp_3dpf_integrated_umap
dm<-dm_ovary
dr<-dr_ovary
hs<-hs_ovary
mm<-mm_ovary
mm_pg<-pg_integrated

pm[["RNA3"]] <- as(object = pm[["RNA"]], Class = "Assay")
DefaultAssay(pm) <- "RNA3"
pm[["RNA"]] <- NULL
pm <- RenameAssays(object = pm, RNA3 = 'RNA')
SaveH5Seurat(pm, filename = "pm_ovary.h5Seurat")
Convert("pm_ovary.h5Seurat", dest = "h5ad")

sp[["RNA3"]] <- as(object = sp[["RNA"]], Class = "Assay")
DefaultAssay(sp) <- "RNA3"
sp[["RNA"]] <- NULL
sp <- RenameAssays(object = sp, RNA3 = 'RNA')
SaveH5Seurat(sp, filename = "sp_adult_ovary.h5Seurat")
Convert("sp_adult_ovary.h5Seurat", dest = "h5ad")

spl[["RNA3"]] <- as(object = spl[["RNA"]], Class = "Assay")
DefaultAssay(spl) <- "RNA3"
spl[["RNA"]] <- NULL
spl <- RenameAssays(object = spl, RNA3 = 'RNA')
SaveH5Seurat(spl, filename = "sp_3df_larva.h5Seurat")
Convert("sp_3df_larva.h5Seurat", dest = "h5ad")

mm[["RNA3"]] <- as(object = mm[["RNA"]], Class = "Assay")
DefaultAssay(mm) <- "RNA3"
mm[["RNA"]] <- NULL
mm <- RenameAssays(object = mm, RNA3 = 'RNA')
SaveH5Seurat(mm, filename = "mm_adult_ovary.h5Seurat")
Convert("mm_adult_ovary.h5Seurat", dest = "h5ad")

dr[["RNA3"]] <- as(object = dr[["RNA"]], Class = "Assay")
DefaultAssay(dr) <- "RNA3"
dr[["RNA"]] <- NULL
dr <- RenameAssays(object = dr, RNA3 = 'RNA')
SaveH5Seurat(dr, filename = "dr_adult_ovary.h5Seurat")
Convert("dr_adult_ovary.h5Seurat", dest = "h5ad")

dm[["RNA3"]] <- as(object = dm[["RNA"]], Class = "Assay")
DefaultAssay(dm) <- "RNA3"
dm[["RNA"]] <- NULL
dm <- RenameAssays(object = dm, RNA3 = 'RNA')
SaveH5Seurat(dm, filename = "dm_adult_ovary.h5Seurat")
Convert("dm_adult_ovary.h5Seurat", dest = "h5ad")

hs[["RNA3"]] <- as(object = hs[["RNA"]], Class = "Assay")
DefaultAssay(hs) <- "RNA3"
hs[["RNA"]] <- NULL
hs <- RenameAssays(object = hs, RNA3 = 'RNA')
SaveH5Seurat(hs, filename = "hs_fetal_ovary.h5Seurat")
Convert("hs_fetal_ovary.h5Seurat", dest = "h5ad")

mm[["RNA3"]] <- as(object = mm[["RNA"]], Class = "Assay")
DefaultAssay(mm) <- "RNA3"
mm[["RNA"]] <- NULL
hs <- RenameAssays(object = mm, RNA3 = 'RNA')
SaveH5Seurat(mm, filename = "mm_pituitary_gland.h5Seurat")
Convert("mm_pituitary_gland.h5Seurat", dest = "h5ad")


#In python

python3

from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
from samalg import SAM
import pandas as pd

% bash map_genes.sh --tr1 /directory of P.miniata proteome fasta file --t1 prot --n1 pm --tr2 //directory of species tested proteome fasta file --t2 prot --n2 xx #xx= spl or dm or dr or hs or mm or sp or mm_pg depending on the comparison

fn1 = "directory of pm_ovary.h5ad"
fn2 = "directory of species compared to"

filenames = {'pm':fn1,'xx':fn2} #xx= spl or dm or dr or hs or mm or sp or mm_pg depending on the comparison

sm = SAMAP(
  filenames,
  f_maps = 'data/maps/',
  save_processed=True #if False, do not save the processed results to `*_pr.h5ad`
)

sm.run(pairwise=True)
samap = sm.samap # SAM object with 2 species stitched together
keys = {'pm':'seurat_clusters','xx':'seurat_clusters'} #xx= spl or dm or dr or hs or mm or sp or mm_pg depending on the comparison
D,MappingTable = get_mapping_scores(sm,keys,n_top = 0)
D.head()
df = MappingTable

print(df)

df.to_csv('mapping_table_sea_star_vs_x_species.csv')

#In R Studio

df <- load mapping scores
data1 <- melt(df)
data1$X <- as.factor(data1$X)

threshold <- 0.1  # threshold value
data1 <- filter(data1, data1$value >= threshold)

# Transform the data for ggalluvial
alluvial_data <- data1 %>%
  mutate(species1 = factor(data1$X, levels = unique(data1$X)),
         species2 = factor(data1$variable, levels = unique(data1$variable)))

ggplot(alluvial_data,
       aes(axis1 = species1, axis2 = species2, y = data1$value)) +
  geom_alluvium(aes(fill = data1$value), width = 0.25) +
  geom_stratum(width = 0.25, fill = "white", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Species 1", "Species 2"), expand = c(0.15, 0.05)) +
  scale_fill_gradient(low = "lightblue", high = "blue", name = "Alignment Score")+
  theme_minimal()+theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  )+theme(
    axis.text.x = element_blank(),  # Remove x-axis labels
    axis.ticks.x = element_blank(), # Remove x-axis ticks
    axis.text.y = element_blank(),  # Remove y-axis labels
    axis.ticks.y = element_blank()  # Remove y-axis ticks
  )

ggsave("Sankey_plot_pm_ovary_vs_x species.tiff", units="in", width=12, height=13
       , dpi=300, compression = 'lzw')

