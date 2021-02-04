# Gene signature of children with severe respiratory syncytial virus infection
 
[Clyde Dapat](https://orcid.org/0000-0002-7616-4680)<sup>1</sup>, Satoru Kumaki<sup>2</sup>, Hiroki Sakurai<sup>3</sup>, Hidekazu Nishimura<sup>4</sup>, Hannah Karen Mina Labayo<sup>1</sup>, Michiko Okamoto<sup>1</sup>, Mayuko Saito<sup>1</sup>, and Hitoshi Oshitani<sup>1</sup>

<sup>1</sup>Department of Virology, Tohoku University Graduate School of Medicine, 2-1 Seiryo-machi, Aoba-ku, Sendai, Japan 980-8575

<sup>2</sup>Department of Pediatrics, Sendai Medical Center, 11-12 Miyagino 2-chome Miyagino-ku Sendai, Miyagi Prefecture, Japan 983-8520

<sup>3</sup>Department of General Pediatrics, Miyagi Children’s Hospital, 3-17 Ochiai 4-chome Aoba-ku Sendai, Miyagi Prefecture, Japan 989-3126

<sup>4</sup>Virus Research Center, Sendai Medical Center, 11-12 Miyagino 2-chome Miyagino-ku Sendai, Miyagi Prefecture, Japan 983-8520

## Abstract

**Background**:  The limited treatment options for children with severe respiratory syncytial virus (RSV) infection highlights the need for a comprehensive understanding of the host cellular response during infection. We aimed to identify host genes that are associated with severe RSV disease and to identify drugs that can be repurposed for the treatment of severe RSV infection. 

**Methods**: We examined clinical data and blood samples from 37 hospitalized children (29 mild and 8 severe) with RSV infection.  We tested RNA from blood samples using next-generation sequencing to profile global mRNA expression and identify cellular processes.

**Results**: Retractions, decreased breath sounds, and tachypnea were associated with disease severity. We observed upregulation of genes related to neutrophil, inflammatory response, blood coagulation and downregulation of genes related to T cell response in children with severe RSV.  Using network-based approach, 43 drugs were identified that are predicted to interact with the gene products of these differentially expressed genes.

**Conclusions**: These results suggest that the changes in the expression pattern in the innate and adaptive immune responses may be associated with RSV clinical severity. Compounds that target these cellular processes can be repositioned as candidate drugs in the treatment of severe RSV. 


## Repository information
This repository accompanies the paper entitled "*Gene signature of children with severe respiratory syncytial virus infection*" **[(Pediatric Research 2021 Jan doi: 10.1038/s41390-020-01347-9)](https://pubmed.ncbi.nlm.nih.gov/33510411/)**. 

The codes and data provided here are for reproducing the differential gene expression analysis from the paper.

The RNA-seq data can be accessed at NCBI Gene Expression Omnibus database with accession number **[GSE155925](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155925)**.


## Software dependencies
The R script was created and run using R version 4.0.3. R packages needed for differential gene expression analysis are listed in the script. 

## Codes
### *differential_expression_analysis.R*
This file contains the R script for gene expression analysis.

## Data
### *raw_counts_matrix.txt*
This file contains the unnormalized RNA-seq count data.

### *patient_data.csv*
This file contains data on sex, age, severity.