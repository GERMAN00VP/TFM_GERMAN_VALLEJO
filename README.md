# TFM German Vallejo

This repository contains the code developed for my bioinformatics Master's final project:

**Transcriptomic Characterization of the Tumor Microenvironment and Chromatin Remodeler BPTF in Pancreatic Ductal Adenocarcinoma**

The project is subdivided into four directories, following the four objectives of my research.

## TCGA CLUSTERS ANALYSIS

This section presents the results of the RNA-Seq data analysis and survival analysis for a cohort of pancreatic ductal adenocarcinoma (PDAC) patients. The analysis is organized into 5 notebooks:

- **Clustering:** Unsupervised Gaussian classification of patients based on their gene expression.

- **Kaplan-Meier_TCGA_clusters:** Survival analysis of the clusters.

- **Differential_expression_and_Functional_analysis:** Differential expression and functional analysis between the clusters.

- **TumorDecon:** Deconvolution of immune cells and comparison of immune infiltration in both clusters.

- **WGCNA:** Analysis of co-expressed gene networks.

## sc_RNA_Seq

Re-analysis of a sc-RNA-Seq dataset obtained from Peng et al. (Cell Research, 2019, 29:725–738; [DOI](https://doi.org/10.1038/s41422-019-0195-y)). The data is not available due to its large size, but the code includes links to the repositories where it was obtained. It is divided into 5 notebooks:

- **Single_Cell_RNA_Seq_analysis:** Preprocessing, normalization, integration, and annotation of the sc-RNA-Seq experiment.

- **Differential_expression_analysis:** Analysis of differential and functional expression of ductal cells vs. malignant ductal cells and endothelial cells from tumor tissue vs. control.

- **TF_Pathways_activities:** Analysis of biological pathway activity and transcription factor activity in ductal cells and malignant ductal cells using decoupleR.

- **Complementary_analysis:** Generation of plots of interest.

- **Risk_signature_analysis:** Analysis of genes of interest from the risk signature.

## RNA-Seq

Analysis of the transcriptomic effects (RNA-Seq) of BPTF gene silencing in conjunction with TNFa treatment in PDAC model cell lines, with samples generated in our laboratory. It contains a notebook within the RNA-Seq pre-analysis folder with the bash commands used in the RNA-Seq analysis pipeline. It also contains another notebook (**RNA_Seq_Analysis**) with downstream analysis.

## Risk_signature

Obtaining a risk signature of BPTF-dependent genes. It contains two notebooks:

- **Risk_signature:** Obtaining and validating the risk signature.

- **Risk_signature_genes_analysis:** Analysis of the prognostic value of each gene in the signature individually.

## Utils:

This folder contains four Python scripts and a pickle file, usefull resources also developed during the project:

- **Volcanoplot:** Generates a customized volcano plot (including highlighted gene names) from the results of differential expression.

- **Functional_analysis:** Contains useful functions for functional analyses, result filtering, and plot generation.

- **gmt_tools:** Contains functions for reading and generating .gmt files, a format commonly used for storing molecular signatures.

- **gtf_dict_GRCh38.110:** Contains a dictionary with the correspondence of ENSEMBL_ID and gene symbol.

- **Risk signature:** Contains functions necessary to obtain the risk signature.
