library(DiffBind)
library(edgeR)
library(BiocParallel)
library(dplyr)
# Set up parallelization to use 4 cores
bp_param <- MulticoreParam(workers = 6)

samples <- data.frame(
  SampleID = c("MSH2KO-1", "MSH2KO-2", "MSH2KO-3", "MSH2R4-1", "MSH2R4-2", "MSH2R4-3"),
  Tissue = c(rep("MSH2KO", 3), rep("MSH2R4", 3)),  #  'Tissue' can represent conditions or leave this column out.
  Factor = c(rep("MSH2KO", 3), rep("MSH2R4", 3)),
  Condition = c(rep("MSH2KO", 3), rep("MSH2R4", 3)),  # The experimental conditions
  Treatment = rep("None", 6),  # Since there's no treatment in this setup
  Replicate = c(1, 2, 3, 1, 2, 3),  # Indicating replicates
  bamReads = c("../6_bigwig/MSH2KO-1_sort_n.bam", "../6_bigwig/MSH2KO-2_sort_n.bam", "../6_bigwig/MSH2KO-3_sort_n.bam",
               "../6_bigwig/MSH2R4-1_sort_n.bam", "../6_bigwig/MSH2R4-2_sort_n.bam", "../6_bigwig/MSH2R4-3_sort_n.bam"),
  Peaks = c("../5_genrich/MSH2KO-1.narrowPeak", "../5_genrich/MSH2KO-2.narrowPeak", "../5_genrich/MSH2KO-3.narrowPeak",
            "../5_genrich/MSH2R4-1.narrowPeak", "../5_genrich/MSH2R4-2.narrowPeak", "../5_genrich/MSH2R4-3.narrowPeak"),
  PeakCaller = rep("narrow", 6),  # For Genrich, 'narrow' for narrowPeak format
  stringsAsFactors = FALSE
)


# Initialize DiffBind object
dbaObj <- NULL
dbaObj <- dba(sampleSheet=samples)
# Count reads in peaks

# Count reads in consensus peaks using 6 cores
dbaObj <- dba.count(dbaObj, summits=250, bParallel = TRUE,minOverlap = 2)

# Define contrasts (KO vs HLA)
dbaObj <- dba.contrast(dbaObj, categories = DBA_CONDITION)

# Perform differential binding analysis
dbaObj <- dba.analyze(dbaObj)
saveRDS(dbaObj,file="../7_diffbind/MSH2_dbaObj.RDATA")
# Generate and inspect the report of differentially bound regions
diffPeaks <- dba.report(dbaObj)
write.csv(diffPeaks, file = "../7_diffbind/MSH2_differential_binding_results.csv")
diffPeaksIn <- read.csv("../7_diffbind/MSH2_differential_binding_results.csv")


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
write_peaks_to_bed(diffPeaksIn, "../7_diffbind/MSH2_differential_binding_results.bed")

#Make a README
# Define column names and their descriptions
column_descriptions <- data.frame(
  Heading = c("X", "seqnames", "start", "end", "width", "strand", "Conc", 
              "Conc_MSH2R4", "Conc_MSH2KO", "Fold", "p.value", "FDR"),
  Description = c(
    "Row index or identifier (not relevant for analysis)",
    "Chromosome or contig name where the peak is located",
    "Start position of the peak on the chromosome",
    "End position of the peak on the chromosome",
    "Width of the peak (end - start)",
    "Strand information (+/-) where the peak is located, if applicable",
    "Average read concentration (normalized counts) across all samples",
    "Average read concentration in MSH2 samples (wild-type or tagged)",
    "Average read concentration in KO (knockout) samples",
    "Log2 fold change of read concentration between MSH2KO and MSH2R4 samples",
    "Raw p-value for differential binding between MSH2KO and MSH2R4 samples",
    "False discovery rate (adjusted p-value) for differential binding"
  ),
  stringsAsFactors = FALSE
)

#MSH2 peaks

# Set thresholds
min_Conc_MSH2R4 <- 2  # Adjust this value based on data and expectations
max_Conc_MSH2KO <- 0.2  # Adjust based on definition of "low" in KO

# Filter peaks for positive fold change, low concentration in KO, and reasonably high concentration in MSH2
MSH2R4_filtered_peaks <- diffPeaksIn %>%
  filter(Fold > 0 & Conc_MSH2KO < max_Conc_MSH2KO & Conc_MSH2R4 > min_Conc_MSH2R4)

# Check the filtered results
write.csv(MSH2R4_filtered_peaks,"../7_diffbind/MSH2R4_filtered_differential_binding_results.csv")
write_peaks_to_bed(MSH2R4_filtered_peaks, "../7_diffbind/MSH2R4_filtered_differential_binding_results.bed")

#MSH2KO peaks

# Set thresholds
min_Conc_MSH2KO <- 2  # Adjust this value based on data and expectations
max_Conc_MSH2R4 <- 0.2  # Adjust based on definition of "low" in MSH2

# Filter peaks for positive fold change, low concentration in KO, and reasonably high concentration in MSH2
MSH2KO_filtered_peaks <- diffPeaksIn %>%
  filter(Fold < 0 & Conc_MSH2R4 < max_Conc_MSH2R4 & Conc_MSH2KO > min_Conc_MSH2KO)

# Check the filtered results
write.csv(MSH2KO_filtered_peaks,"../7_diffbind/MSH2KO_filtered_differential_binding_results.csv")
write_peaks_to_bed(MSH2KO_filtered_peaks, "../7_diffbind/MSH2KO_filtered_differential_binding_results.bed")

#From this point, you can view filtered peaks and per-sample reads in IGV by loading:
#"../7_diffbind/filtered_differential_binding_results.bed"
#"../6_bigwig/MSH2KO-1.bw","../6_bigwig/MSH2KO-2.bw","../6_bigwig/MSH2KO-3.bw","../6_bigwig/MSH2R4-1.bw","../6_bigwig/MSH2R4-2.bw","../6_bigwig/MSH2R4-3.bw")
#Write all the results to an excel

##Annotate peaks with chipseeker
# Load necessary libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # You can change the genome version accordingly
library(clusterProfiler)

# Load the peaks (ensure peaks are in BED format)
peaks <- readPeakFile("../7_diffbind/MSH2_differential_binding_results.bed")

library(GenomicRanges)

# Ensure peaks are in the GRanges format, then add the "chr" prefix
seqlevelsStyle(peaks) <- "UCSC"  # This will add "chr" to the chromosome names

# Annotate the peaks to genes
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Now, run the annotation
peak_annotation <- annotatePeak(peaks, tssRegion=c(-3000, 3000), TxDb=txdb)

# Save the annotated results
peak_annotation_df <- as.data.frame(peak_annotation)

# Add column descriptions for the annotation results
column_descriptions_anno <- data.frame(
  Column = c(
    "seqnames", "start", "end", "width", "strand", 
    "V4", "V5", "annotation", "geneChr", "geneStart", 
    "geneEnd", "geneLength", "geneStrand", "geneId", 
    "transcriptId", "distanceToTSS"
  ),
  Description = c(
    "Chromosome name where the peak is located",
    "Start position of the peak",
    "End position of the peak",
    "Width of the peak (end - start)",
    "Strand of the peak (+ or -)",
    "Additional metadata (e.g., score or rank of peak)", 
    "Additional metadata (e.g., signal value or score)",
    "Genomic annotation of the peak (e.g., promoter, exon, intron, intergenic)",
    "Chromosome where the nearest gene is located",
    "Start position of the nearest gene",
    "End position of the nearest gene",
    "Length of the nearest gene",
    "Strand of the nearest gene (+ or -)",
    "Ensembl or RefSeq ID of the nearest gene",
    "Ensembl or RefSeq ID of the nearest transcript",
    "Distance from the peak to the transcription start site (TSS) of the nearest gene"
  ),
  stringsAsFactors = FALSE
)


library(openxlsx)
# Create a new workbook
wb <- createWorkbook()

# List of data frames and their corresponding names
dfs <- list(README = column_descriptions,
            differential_peaks = diffPeaksIn,
            MLH1R4_filtered_peaks = MSH2R4_filtered_peaks,
            MLH1KO_filtered_peaks = MSH2KO_filtered_peaks)

# Add each data frame to the workbook as a separate sheet
for (df_name in names(dfs)) {
  addWorksheet(wb, df_name)  # Create a sheet named after the data frame
  writeData(wb, df_name, dfs[[df_name]])  # Write the data to the sheet
}

# Save the workbook to an Excel file
saveWorkbook(wb, "../7_diffbind/MSH2_Differential_binding_results.xlsx", overwrite = TRUE)

