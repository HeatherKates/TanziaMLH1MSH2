# Load necessary libraries
library(tximport)
library(DESeq2)
library(biomaRt)

# Define the directory containing the quantification results
dir <- "/blue/zhangw/hkates/Tanzia_RNAseq/results/salmon/"

# Get the list of quantification directories starting with "MLH1"
dirs <- list.files(dir, pattern="^MSH2.*_quant$", full.names=TRUE)

# Create a named vector of files for tximport
files <- file.path(dirs, "quant.sf")
names(files) <- basename(dirs)

# Define the conditions based on filenames
sample_names <- names(files)
conditions <- ifelse(grepl("KO", sample_names), "KO", "HA")

# Create colData dataframe
colData <- data.frame(
  row.names = sample_names,
  condition = factor(conditions, levels = c("HA", "KO"))
)

# Use biomaRt to get the transcript-to-gene mapping
options(timeout = 600)  # 10 minutes
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
tx2gene <- getBM(attributes=c("ensembl_transcript_id_version", "ensembl_gene_id"),
                 mart=mart)

# Import quantification data using tximport and summarize to gene level
txi <- tximport(files, type="salmon", tx2gene=tx2gene, txOut=FALSE)

# Create DESeq2 dataset
dds <- DESeqDataSetFromTximport(txi, colData, design=~condition)

# replace with gene names
# Get the gene names for the Ensembl gene IDs
gene_ids <- rownames(dds)
genes <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"), 
               filters="ensembl_gene_id", 
               values=gene_ids, 
               mart=mart)

# Ensure that all gene IDs are mapped
genes <- genes[match(gene_ids, genes$ensembl_gene_id),]

# Store the mapping in a separate data frame
gene_mapping <- data.frame(ensembl_gene_id = gene_ids, 
                           external_gene_name = genes$external_gene_name)


# Replace row names in DESeq2 dataset
rownames(dds) <- genes$external_gene_name

# Verify the change
rownames(dds)[1:10]

# Differential expression analysis
dds <- DESeq(dds)
res <- results(dds)

# Save results
MSH2_dds <- dds
MSH2_res <- res
MSH2_gene_mapping <- gene_mapping

write.csv(as.data.frame(MSH2_res), file="/blue/zhangw/hkates/Tanzia_RNAseq/results/deseq2/MSH2_DESeq2_gene_results.csv")
save.image(file="/blue/zhangw/hkates/Tanzia_RNAseq/results/deseq2/MSH2_DESeq2_result.RDATA")

saveRDS(MSH2_dds,"/blue/zhangw/hkates/Tanzia_RNAseq/results/deseq2/MSH2_dds.Rds")
saveRDS(MSH2_res,"/blue/zhangw/hkates/Tanzia_RNAseq/results/deseq2/MSH2_res.Rds")
saveRDS(MSH2_gene_mapping,"/blue/zhangw/hkates/Tanzia_RNAseq/results/deseq2/MSH2_gene_mapping.Rds")
