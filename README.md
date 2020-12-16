# Gene signature of children with severe respiratory syncytial virus infection
 
[Clyde Dapat](https://orcid.org/0000-0002-7616-4680)<sup>1</sup>, Satoru Kumaki<sup>2</sup>, Hiroki Sakurai<sup>3</sup>, Hidekazu Nishimura<sup>4</sup>, Hannah Karen Mina Labayo<sup>1</sup>, Michiko Okamoto<sup>1</sup>, Mayuko Saito<sup>1</sup>, and Hitoshi Oshitani<sup>1</sup>

<sup>1</sup>Department of Virology, Tohoku University Graduate School of Medicine, 2-1 Seiryo-machi, Aoba-ku, Sendai, Japan 980-8575
<sup>2</sup>Department of Pediatrics, Sendai Medical Center, 11-12 Miyagino 2-chome Miyagino-ku Sendai, Miyagi Prefecture, Japan 983-8520
<sup>3</sup>Department of General Pediatrics, Miyagi Children’s Hospital, 3-17 Ochiai 4-chome Aoba-ku Sendai, Miyagi Prefecture, Japan 989-3126
<sup>4</sup>Virus Research Center, Sendai Medical Center, 11-12 Miyagino 2-chome Miyagino-ku Sendai, Miyagi Prefecture, Japan 983-8520

## Repository information
This repository accompanies the paper entitled "*Gene signature of children with severe respiratory syncytial virus infection*" submitted to the journal **Pediatric Research**. 

The code and data provided here are for reproducing the differential gene expression analysis from the paper.

The RNA-seq data can be accessed at NCBI Gene Expression Omnibus database with accession number **GSE155925**.


## Software dependencies
The R script was created and run using R version 4.0.3. R packages needed for differential gene expression analysis are listed in the script. 

## Code
### *differential_expression_analysis.R*
This file contains the R script for gene expression analysis.

## Data
### *raw_counts_matrix.txt*
This file contains the unnormalized RNA-seq count data.

### *patient_data.csv*
This file contains data on sex, age, severity.

