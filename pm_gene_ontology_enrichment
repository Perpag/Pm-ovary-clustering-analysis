library(topGO)
library(clusterProfiler)
library(ggplot2)
library(dplyr)
library(tidyr)

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
