#!/bin/bash
module load homer

# Add "chr" prefix to the chromosome column in the BED file
sed 's/^/chr/' ../7_diffbind/MSH2R4_filtered_differential_binding_results.bed > ../9_HOMER/MSH2R4_filtered_differential_binding_results_chr.bed
# Add "chr" prefix to the chromosome column in the BED file
sed 's/^/chr/' ../7_diffbind/MSH2KO_filtered_differential_binding_results.bed > ../9_HOMER/MSH2KO_filtered_differential_binding_results_chr.bed

findMotifsGenome.pl ../9_HOMER/MSH2KO_filtered_differential_binding_results_chr.bed hg38 ../9_HOMER/MSH2KO_peaks_motifs -size 200
findMotifsGenome.pl ../9_HOMER/MSH2R4_filtered_differential_binding_results_chr.bed hg38 ../9_HOMER/MSH2R4_peaks_motifs -size 200

