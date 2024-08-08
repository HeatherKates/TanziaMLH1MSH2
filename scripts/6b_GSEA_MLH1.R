# Load the libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)

#load data
MLH1_res <- readRDS("/blue/zhangw/hkates/Tanzia_RNAseq/results/deseq2/MLH1_res.Rds")
MLH1_dds <- readRDS("/blue/zhangw/hkates/Tanzia_RNAseq/results/deseq2/MLH1_dds.Rds")


# Convert DESeqMLH1_results to a data frame
MLH1_res_df <- as.data.frame(MLH1_res)
MLH1_res_df$external_gene_name <- rownames(MLH1_res_df)

#read in gene mapping df
gene_mapping <- readRDS("/blue/zhangw/hkates/Tanzia_RNAseq/results/deseq2/MLH1_gene_mapping.Rds")

# Merge with the gene mapping to retain ENSEMBL IDs and gene names
MLH1_res_df <- merge(MLH1_res_df, gene_mapping, by = "external_gene_name")

# Create a named vector of log2FoldChange values with ENSEMBL IDs
geneList <- MLH1_res_df$log2FoldChange
names(geneList) <- MLH1_res_df$ensembl_gene_id

# Sort the gene list in decreasing order
geneList <- sort(geneList, decreasing = TRUE)

gseaGO <- gseGO(geneList, 
                OrgDb = org.Hs.eg.db, 
                keyType = "ENSEMBL", 
                ont = "ALL", 
                pvalueCutoff = 0.05, 
                verbose = TRUE)

# View the MLH1_results
head(gseaGO)

# Visualize the MLH1_results
dotplot(gseaGO)

