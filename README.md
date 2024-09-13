# RNAseq, ATACseq, and CUT&RUN Analysis for MSH2 and MLH1 in Tumor Metastasis

This repository contains scripts and files related to the analysis of RNAseq, ATACseq, and CUT&RUN data to investigate the contrasting roles of MSH2 and MLH1 in tumor metastasis.

## Overview

The code provided here is fully reproducible. Except for the original `fastq.gz` files, all other files and results can be recreated using the scripts included in this repository. Due to size constraints and data ownership, large files and non public files are omitted. For access to these files, please see the `.gitignore` and email [hkates@ufl.edu](mailto:hkates@ufl.edu).

## Directory Structure

- **RNAseq**
  * **scripts/**: Contains all the scripts necessary to reproduce the analysis.  
    * script files are named beginning with a sequential order that they can be run in for full reproducibility. (A comes before B).
  * **results/**: Directory for storing analysis results.

- **ATACseq**
  * **scripts/**
  * **results/**
  
## Reproducibility

To ensure full reproducibility, follow the steps below:

1. **Download Original Data**: Obtain the original `fastq.gz` files or run scripts on hipergator where they will access data in /orange (persmission is for users in zhangw only)
2. **Run Scripts**: Use the scripts provided in the `scripts` directory to process the data and generate results.

## Contact

For access to large files or any other inquiries, please contact:

Heather Kates  
Email: [hkates@ufl.edu](mailto:hkates@ufl.edu)

## License

None

## Acknowledgements

This research was supported by [unknown]
