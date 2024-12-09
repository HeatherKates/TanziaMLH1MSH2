library(edgeR)
library(BiocParallel)
register(SerialParam())
library(dplyr)
# Set up parallelization to use 4 cores
bp_param <- MulticoreParam(workers = 6)
library(DiffBind)

samples <- read.csv("diffbind_samplesheet.csv")
samples$SampleID <- gsub("../","",samples$SampleID)
samples$Condition <- gsub("../","",samples$Condition)
samples <- samples %>% filter(Factor=="HA")
samples <- samples %>% select(!Tissue.1)
# Initialize DiffBind object
dbaObj <- NULL
dbaObj <- dba(sampleSheet=samples)
# Count reads in consensus peaks using 6 cores
dbaObj <- dba.count(dbaObj, summits=250, bParallel = TRUE,minOverlap = 1)
# Define contrasts (R4 vs. KO)
dbaObj <- dba.contrast(dbaObj, minMembers=2, reorderMeta=list(Condition="KO"))
# Perform differential binding analysis
dbaObj <- dba.analyze(dbaObj)
saveRDS(dbaObj,file="../results/01_diffbind/MLH1_dbaObj.RDATA")
# Generate and inspect the report of differentially bound regions
diffPeaks <- dba.report(dbaObj,th=1)
write.csv(diffPeaks, file = "../results/01_diffbind/MLH1_differential_binding_results.csv")
diffPeaksIn <- read.csv("../results/01_diffbind/MLH1_differential_binding_results.csv")

# Subtract KO peaks from R4 peaks:
library(GenomicRanges)
library(rtracklayer)

# Load called peaks for R4 and KO
KO_peaks_1 <- import(samples$Peaks[[1]])
KO_peaks_2 <- import(samples$Peaks[[2]])
R4_peaks_1 <- import(samples$Peaks[[3]])
R4_peaks_2 <- import(samples$Peaks[[4]])

# Combine KO peaks into a single set (non-redundant union)
combined_KO_peaks <- GenomicRanges::reduce(c(KO_peaks_1, KO_peaks_2))

# Combine R4 peaks into a single set (non-redundant union)
combined_R4_peaks <- reduce(c(R4_peaks_1, R4_peaks_2))

# Subtract KO peaks from R4 peaks to get unique R4 peaks
unique_r4_peaks <- setdiff(combined_R4_peaks, combined_KO_peaks)

# Export the unique R4 peaks to a new BED file
export(unique_r4_peaks, "../results/01_diffbind/MLH1_unique_peaks_results.bed")

# Function to write peaks data frame to a BED file
write_peaks_to_bed <- function(peaks_df, output_file) {
  
  # Ensure required columns are present
  if (!all(c("seqnames", "start", "end", "X", "Fold") %in% colnames(peaks_df))) {
    stop("The data frame does not contain the necessary columns: 'seqnames', 'start', 'end', 'X', 'Fold'")
  }
  
  # Create the BED format data frame
  bed <- data.frame(
    chr = peaks_df$seqnames,      # Chromosome column
    start = peaks_df$start,       # Start position
    end = peaks_df$end,           # End position
    name = peaks_df$X,            # Peak ID or use `.` if no name is required
    score = peaks_df$Fold         # Log2 fold-change as score
  )
  
  # Write the BED file
  write.table(bed, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  message("BED file written to: ", output_file)
}

# Example usage with your dataframe 'diffPeaksIn'
write_peaks_to_bed(diffPeaksIn, "../results/01_diffbind/MLH1R4_filtered_differential_binding_results.bed")

MLH1R4_filtered_peaks <- diffPeaksIn 

# Write the results
write.csv(MLH1R4_filtered_peaks,"../results/01_diffbind/MLH1R4_filtered_differential_binding_results.csv")
#write_peaks_to_bed(MLH1R4_filtered_peaks, "../results/01_diffbind/MLH1R4_filtered_differential_binding_results.bed")

#From this point, you can view filtered peaks and per-sample reads in IGV by loading:
#"../7_diffbind/filtered_differential_binding_results.bed"
#"../6_bigwig/MLH1KO-1.bw","../6_bigwig/MLH1KO-2.bw","../6_bigwig/MLH1KO-3.bw","../6_bigwig/MLH1R4-1.bw","../6_bigwig/MLH1R4-2.bw","../6_bigwig/MLH1R4-3.bw")
#Write all the results to an excel

##Annotate peaks with chipseeker
# Load necessary libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # You can change the genome version accordingly
library(clusterProfiler)

# Load the peaks (ensure peaks are in BED format)
peaks <- readPeakFile("../results/01_diffbind/MLH1R4_filtered_differential_binding_results.bed")
peaks <- unique_r4_peaks
library(GenomicRanges)

# Ensure peaks are in the GRanges format, then add the "chr" prefix
seqlevelsStyle(peaks) <- "UCSC"  # This will add "chr" to the chromosome names

# Annotate the peaks to genes
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Now, run the annotation
peak_annotation <- annotatePeak(peaks, tssRegion=c(-3000, 3000), TxDb=txdb)


# Add gene symbols

# Load the necessary library for gene symbol mapping
library(org.Hs.eg.db)  # For human genes, change if using another species

# Extract the Entrez IDs from the annotation result
entrez_ids <- as.data.frame(peak_annotation)$geneId

# Map Entrez IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db, keys = entrez_ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")

# Save the annotated results
peak_annotation_df <- as.data.frame(peak_annotation)

# Add the gene symbols to the peak_annotation result
peak_annotation_df$geneSymbol <- gene_symbols

#Add peak annotation to diffbind results
# Add a "X" column
peak_annotation_df$X=rownames(peak_annotation_df)

# Add DE results

#Read in DE results
DE_results <- read.csv("../../../RNAseq/results/4-6_deseq2/MLH1_DESeq2_gene_results.csv")

# Specify the direction of the log2FoldChange
DE_results <- DE_results %>% dplyr::rename(log2FoldChange_KOvsR4=log2FoldChange)
DE_results$log2FoldChange_R4vsKO <- -DE_results$log2FoldChange_KOvsR4
DE_results <- DE_results %>% dplyr::select(-log2FoldChange_KOvsR4)
# Rename ensemble ID col
DE_results <- DE_results %>% dplyr::rename(ensembl_gene_id=X)

# Add entrez ID for merge with peak annotation
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Extract the Ensembl Gene IDs from DE_results dataframe
ensembl_ids <- DE_results$ensembl_gene_id

# Use biomaRt to retrieve the corresponding Entrez Gene IDs
gene_info <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                   filters = "ensembl_gene_id", values = ensembl_ids, mart = ensembl)

# Merge the Entrez Gene IDs with your DE_results dataframe based on the Ensembl Gene IDs
DE_results <- merge(DE_results,gene_info,by="ensembl_gene_id",all.x=TRUE)

# Add DE_result to all colnames except the last 2 col which is geneSymbol and entrez ID
colnames(DE_results)[1:7] <- paste0("DE_result.",colnames(DE_results)[1:7])

# Merge DE result with peaks
peaks_with_anno <- merge(peak_annotation_df,DE_results,by.y="entrezgene_id",by.x="geneId",all.x=TRUE)

# Rename "X" to be more informative
peaks_with_anno <- peaks_with_anno %>% dplyr::rename(peak_ID=X)

# Remove suffix x from merges
colnames(peaks_with_anno) <- gsub(".x$","",colnames(peaks_with_anno))

# Remove chipseeker cols redundant with macs2 cols
peaks_with_anno <- peaks_with_anno %>% dplyr::select(c(-geneChr,-geneStart,-geneEnd))
## Add README

# Create a data frame with the column names, descriptions, and the source of the data
README_df <- data.frame(
  Column = colnames(peaks_with_anno),
  Description = c(
    "Gene ID from annotation",
    "Peak ID from macs2 peak calling",
    "Start position of the peak (macs2)",
    "End position of the peak (macs2)",
    "Width of the peak (macs2)",
    "Strand of the peak (macs2)",
    "Peak annotation (e.g., promoter, exon, intron) from ChIPseeker",
    "Length of the annotated gene (ChIPseeker)",
    "Strand of the annotated gene (ChIPseeker)",
    "Transcript ID from ChIPseeker annotation",
    "Distance from the peak to the nearest TSS (ChIPseeker)",
    "Gene symbol from ChIPseeker annotation (ChIPseeker)",
    "Peak ID (macs2)",
    "Ensembl gene ID from DESeq2 results (DESeq2)",
    "Base mean expression level (DESeq2)",
    "Standard error of the log2 fold change (DESeq2)",
    "Statistical test statistic (DESeq2)",
    "P-value for differential expression (DESeq2)",
    "FDR-adjusted p-value for differential expression (DESeq2)",
    "Gene symbol from DESeq2 results",
    "Log2 fold change from DE analysis"
  ),
  `Analysis that produced field` = c(
    "ChIPseeker",
    "macs2",
    "macs2",
    "macs2",
    "macs2",
    "macs2",
    "ChIPseeker",
    "ChIPseeker",
    "ChIPseeker",
    "ChIPseeker",
    "ChIPseeker",
    "ChIPseeker",
    "macs2",
    "DESeq2",
    "DESeq2",
    "DESeq2",
    "DESeq2",
    "DESeq2",
    "DESeq2",
    "DESeq2",
    "DESeq2"
  ),
  stringsAsFactors = FALSE
)


library(writexl)

dfs <- list(
  README = README_df,
  peaks_with_anno = peaks_with_anno
)

#Add the gene symbol to the annotation
# Write the list of data frames to an Excel file using writexl
write_xlsx(dfs, "../results/01_diffbind/MLH1R4_peaks.xlsx")

