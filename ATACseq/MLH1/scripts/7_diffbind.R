library(DiffBind)
library(edgeR)
library(BiocParallel)
library(dplyr)

# Set up parallelization to use 4 cores
bp_param <- MulticoreParam(workers = 6)

samples <- data.frame(
  SampleID = c("MLH1KO-1", "MLH1KO-2", "MLH1KO-3", "MLH1R4-1", "MLH1R4-2", "MLH1R4-3"),
  Tissue = c(rep("MLH1KO", 3), rep("MLH1R4", 3)),  #  'Tissue' can represent conditions or leave this column out.
  Factor = c(rep("MLH1KO", 3), rep("MLH1R4", 3)),
  Condition = c(rep("MLH1KO", 3), rep("MLH1R4", 3)),  # The experimental conditions
  Treatment = rep("None", 6),  # Since there's no treatment in this setup
  Replicate = c(1, 2, 3, 1, 2, 3),  # Indicating replicates
  bamReads = c("../6_bigwig/MLH1KO-1_sort_n.bam", "../6_bigwig/MLH1KO-2_sort_n.bam", "../6_bigwig/MLH1KO-3_sort_n.bam",
               "../6_bigwig/MLH1R4-1_sort_n.bam", "../6_bigwig/MLH1R4-2_sort_n.bam", "../6_bigwig/MLH1R4-3_sort_n.bam"),
  Peaks = c("../5_genrich/MLH1KO-1.narrowPeak", "../5_genrich/MLH1KO-2.narrowPeak", "../5_genrich/MLH1KO-3.narrowPeak",
            "../5_genrich/MLH1R4-1.narrowPeak", "../5_genrich/MLH1R4-2.narrowPeak", "../5_genrich/MLH1R4-3.narrowPeak"),
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
saveRDS(dbaObj,file="../7_diffbind/MLH1_dbaObj.RDATA")

# Generate and inspect the report of differentially bound regions
diffPeaks <- dba.report(dbaObj)
head(diffPeaks)
write.csv(diffPeaks, file = "../7_diffbind/MLH1_differential_binding_results.csv")
diffPeaksIn <- read.csv("../7_diffbind/MLH1_differential_binding_results.csv")


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
write_peaks_to_bed(diffPeaksIn, "../7_diffbind/MLH1_differential_binding_results.bed")

#Make a README
# Define column names and their descriptions
column_descriptions <- data.frame(
  Heading = c("X", "seqnames", "start", "end", "width", "strand", "Conc", 
              "Conc_MLH1", "Conc_KO", "Fold", "p.value", "FDR"),
  Description = c(
    "Row index or identifier (not relevant for analysis)",
    "Chromosome or contig name where the peak is located",
    "Start position of the peak on the chromosome",
    "End position of the peak on the chromosome",
    "Width of the peak (end - start)",
    "Strand information (+/-) where the peak is located, if applicable",
    "Average read concentration (normalized counts) across all samples",
    "Average read concentration in MLH1 samples (wild-type or tagged)",
    "Average read concentration in KO (knockout) samples",
    "Log2 fold change of read concentration between KO and MLH1 samples",
    "Raw p-value for differential binding between KO and MLH1 samples",
    "False discovery rate (adjusted p-value) for differential binding"
  ),
  stringsAsFactors = FALSE
)

#Filter peaks only in MLH1R4 (vs MLH1KO)

# Set thresholds
min_Conc_MLH1 <- 2  # Adjust this value based on data and expectations
max_Conc_MLH1KO <- 0.2  # Adjust based on definition of "low" in KO

# Filter peaks for positive fold change, low concentration in KO, and reasonably high concentration in MLH1
MLH1R4_filtered_peaks <- diffPeaksIn %>%
  filter(Fold > 0 & Conc_MLH1KO < max_Conc_MLH1KO & Conc_MLH1 > min_Conc_MLH1)

# Check the filtered results
write.csv(MLH1R4_filtered_peaks,"../7_diffbind/MLH1R4_filtered_differential_binding_results.csv")
write_peaks_to_bed(MLH1R4_filtered_peaks, "../7_diffbind/MLH1R4_filtered_differential_binding_results.bed")

#From this point, you can view filtered peaks and per-sample reads in IGV by loading:
#"../7_diffbind/filtered_differential_binding_results.bed"
#"../6_bigwig/MLH1KO-1.bw","../6_bigwig/MLH1KO-2.bw","../6_bigwig/MLH1KO-3.bw","../6_bigwig/MLH1R4-1.bw","../6_bigwig/MLH1R4-2.bw","../6_bigwig/MLH1R4-3.bw")

#Filter peaks only in MLH1KO

# Set thresholds
min_Conc_MLH1KO <- 2  # Adjust this value based on data and expectations
max_Conc_MLH1 <- 0.2  # Adjust based on definition of "low" in KO

# Filter peaks for positive fold change, low concentration in KO, and reasonably high concentration in MLH1
MLH1KO_filtered_peaks <- diffPeaksIn %>%
  filter(Fold < 0 & Conc_MLH1 < max_Conc_MLH1 & Conc_MLH1KO > min_Conc_MLH1KO)

# Check the filtered results
write.csv(MLH1KO_filtered_peaks,"../7_diffbind/MLH1KO_filtered_differential_binding_results.csv")
write_peaks_to_bed(MLH1KO_filtered_peaks, "../7_diffbind/MLH1KO_filtered_differential_binding_results.bed")

#From this point, you can view filtered peaks and per-sample reads in IGV by loading:
#"../7_diffbind/filtered_differential_binding_results.bed"
#"../6_bigwig/MLH1KO-1.bw","../6_bigwig/MLH1KO-2.bw","../6_bigwig/MLH1KO-3.bw","../6_bigwig/MLH1R4-1.bw","../6_bigwig/MLH1R4-2.bw","../6_bigwig/MLH1R4-3.bw")

#Write all the results to an excel
library(openxlsx)
# Create a new workbook
wb <- createWorkbook()

# List of data frames and their corresponding names
dfs <- list(README = column_descriptions,
            differential_peaks = diffPeaksIn,
            MLH1R4_filtered_peaks = MLH1R4_filtered_peaks,
            MLH1KO_filtered_peaks = MLH1KO_filtered_peaks)

# Add each data frame to the workbook as a separate sheet
for (df_name in names(dfs)) {
  addWorksheet(wb, df_name)  # Create a sheet named after the data frame
  writeData(wb, df_name, dfs[[df_name]])  # Write the data to the sheet
}

# Save the workbook to an Excel file
saveWorkbook(wb, "../7_diffbind/MLH1_Differential_binding_results.xlsx", overwrite = TRUE)
