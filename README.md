## Overview
This repository provides a multi-omics analysis pipeline designed for performing differential analysis on bulk sequencing and proteomics data. Each omic type can be run independently.

## Supported Omics
| Omic type | Mapping / Processing | QC | Differential analysis |
|----------|---------------------|----|----------------------|
| Methylomics | Bismark | methylKit | methylKit |
| Chromatin accessibility (ATAC-seq) | Bowtie2 + MACS3 | ATACseqQC | DESeq2 |
| Transcriptomics (RNA-seq) | Bowtie2 | — | DESeq2 |
| Proteomics | MaxQuant | — | limma |

## Visualization
The pipeline generates standard and omic-specific visualizations, including:
- Principal Component Analysis (PCA)
- Volcano plot
- Heatmap
- Coverage and fragment size distribution plots (ATAC-seq)
- Coverage and methylation level distributions

## Status
This project is under active development as part of my PhD research.
Some analysis modules are still in progress.