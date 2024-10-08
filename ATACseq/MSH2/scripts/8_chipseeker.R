# Load necessary libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # You can change the genome version accordingly
library(clusterProfiler)

# Load the peaks (ensure peaks are in BED format)
peaks_group1 <- readPeakFile("../7_diffbind/MSH2KO_filtered_differential_binding_results.bed")
peaks_group2 <- readPeakFile("../7_diffbind/MSH2R4_filtered_differential_binding_results.bed")

library(GenomicRanges)

# Ensure peaks are in the GRanges format, then add the "chr" prefix
seqlevelsStyle(peaks_group1) <- "UCSC"  # This will add "chr" to the chromosome names
seqlevelsStyle(peaks_group2) <- "UCSC"

# Annotate the peaks to genes
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Now, run the annotation
peak_annotation_group1 <- annotatePeak(peaks_group1, tssRegion=c(-3000, 3000), TxDb=txdb)
peak_annotation_group2 <- annotatePeak(peaks_group2, tssRegion=c(-3000, 3000), TxDb=txdb)

# Visualize the genomic distribution of peaks
plotAnnoBar(peak_annotation_group1)
plotAnnoBar(peak_annotation_group2)

# Save the annotated results
write.table(as.data.frame(peak_annotation_group1), file="../8_chipseek/MSH2KO_annotated_peaks.csv", sep=",", row.names=FALSE)
write.table(as.data.frame(peak_annotation_group2), file="../8_chipseek/MSH2R4_annotated_peaks.csv", sep=",", row.names=FALSE)

