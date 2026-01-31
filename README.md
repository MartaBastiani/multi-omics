## Overview
This repository provides a multi-omics analysis pipeline designed for bulk sequencing and proteomics data. Each omic type can be run independently.

## Supported Omics

| Omic type | Mapping / Processing | QC | Differential analysis |
|----------|---------------------|----|----------------------|
| Methylomics | Bismark | Bismark QC | methylKit |
| Chromatin accessibility (ATAC-seq) | Bowtie2 | ATACseqQC | DESeq2 |
| Transcriptomics (RNA-seq) | Bowtie2 | — | DESeq2 |
| miRNomics | COMPSRA | COMPSRA | limma |
| Proteomics | MaxQuant | — | limma |

## Visualization
The pipeline generates standard and omic-specific visualizations, including:
- Principal Component Analysis (PCA)
- Volcano plots
- Coverage and peak distribution plots (ATAC-seq)
- Methylation level distributions

## Status
This project is under active development as part of my PhD research.
Some analysis modules are still in progress.