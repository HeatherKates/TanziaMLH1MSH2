#!/bin/bash

# Load Salmon module
module load salmon/1.10.1

# Create the directory for the Salmon index
mkdir -p salmon_index
path="salmon_index"

# Download the transcriptome and move it to the designated directory
wget -P $path ftp://ftp.ensembl.org/pub/release-104/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz
gunzip $path/Mus_musculus.GRCm39.cdna.all.fa.gz

# Download the GTF file and move it to the designated directory
wget -P $path ftp://ftp.ensembl.org/pub/release-104/gtf/mus_musculus/Mus_musculus.GRCm39.104.gtf.gz
gunzip $path/Mus_musculus.GRCm39.104.gtf.gz

# Create the Salmon index in the designated directory
salmon index -t $path/Mus_musculus.GRCm39.cdna.all.fa -i $path/salmon_index_mouse
