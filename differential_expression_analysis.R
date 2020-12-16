### R script for differential expression analysis 

## Load R packages
library(DESeq2)
library(vsn)
library(pcaExplorer)
library(ggplot2)
library(sva)
library(IHW)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(circlize)
library(EnrichmentBrowser)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(tmod)


### 1. PREPARING RNA-SEQ COUNT DATA AND PATIENT INFORMATION
## Load patient data
patients <- read.csv("patient_data.csv")
dim(patients)
colnames(patients)
head(patients)
# Assign sample code to rownames
rownames(patients) <- patients$Accession

## Convert variable as factor
# Sex (0=female, 1=male)
patients$Sex <- factor(patients$Sex, levels = c(0,1))
# Age (in months: 0="0to<12", 1="12to<24", 2="24to<36", 3="36to<60")
patients$Age <- factor(patients$Age, levels = c(0,1,2,3))
# RSV season (0=first season, 1=second season)
patients$Batch <- factor(patients$Batch, levels = c(0,1))
# Hospital (0=hospital A, 1=hospital B)
patients$Hospital <- factor(patients$Hospital, levels = c(0,1))
# WHO LRTI (1="RTI", 2="LRTI", 3="severe LRTI", 4="very severe LRTI")
patients$WHO_LRTI <- factor(patients$WHO_LRTI, levels = c(1,2,3,4))
# Severity (0=mild, 1=severe)
patients$Severity <- factor(patients$Severity, levels = c("mild", "severe"))

## Load unnormalized RNA-seq count data
counts <- read.table("raw_counts_matrix.txt", header = TRUE)
dim(counts)
colnames(counts)
head(counts)[1:4]

# Check if sample names match in columns of counts and rows in patient
all(rownames(patients) %in% colnames(counts))
# Check if order also matches
all(rownames(patients) == colnames(counts))


### 2. PREPARING DESEQ OBJECT 
## Construct the DESeqDataSet model (Unadjusted)
rsv <- DESeqDataSetFromMatrix(countData = counts,
                              colData = patients,
                              design = ~ Severity)

# Check levels of factor to ensure that the first level is the reference
factor(rsv$Severity)
# Check the counts in each group
table(rsv$Severity)


### 3. PRINCIPAL COMPONENT ANALYSIS (PCA)
## Transform raw counts to log2 scale 
## using variance stabilizing transformation (VST)
## for clustering analysis and visualization
vsd <- vst(rsv, blind = FALSE)
head(assay(vsd), 3)

# Visualize the effect of transformation on the variance
meanSdPlot(assay(vsd))

# Perform PCA using pcaExplorer package 
# to check for clustering of samples with respect to severity, age, sex, etc. 
pcaExplorer(dds=rsv, dst=vsd)

# Generate PCA plot
pcaData <- pcaplot(vsd, intgroup = "Severity", ntop = 13399, returnData = TRUE)
colnames(pcaData)
ggplot(pcaData, aes(x=PC1, y=PC2)) + 
  scale_color_manual(values = c("black","black")) + 
  geom_point(shape=21, size = 8, aes(fill = Severity)) +
  scale_fill_manual(name = "group", labels = c("mild", "severe"), values =c("#ffffbf","#2b83ba"))  + 
  # axis labels
  labs(x = "PC1 (20% variance)", y = "PC2 (17% variance)") +
  # modify legends
  theme_bw(base_size = 28) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1)) 


### 4. DIFFERENTIAL EXPRESSION ANALYSIS
## Perform differential gene expression using DESeq2 package
rsv <- DESeq(rsv)
res <- results(rsv)
summary(res)

## Perform batch correction using SVA package
dat <- counts(rsv, normalized = TRUE)
idx <- rowMeans(dat) >1
dat <- dat[idx, ]
dim(dat)
# Create a full model matrix
mod <- model.matrix(~ Severity, colData(rsv))
# Create a null model using only the intercept
mod0 <- model.matrix(~ 1, colData(rsv))
# Estimate using surrogate variables
svseq <- svaseq(dat, mod, mod0, n.sv = 2)
svseq$sv
# Add surrogate variables to the data
rsv$SV1 <- svseq$sv[ ,1]
rsv$SV2 <- svseq$sv[ ,2]

# Modify the design matrix to adjust for age, sex, and batch effects
design(rsv) <- ~ SV1 + SV2 + Sex + Age + Hospital + Batch + Severity 

## Rerun DESeq with batch correction
rsv_sv <- DESeq(rsv)
resva <- results(rsv_sv)
summary(resva)

## p-value correction using independent hypothesis weighting (IHW) package
res_rsv_IHW <- results(rsv_sv, filterFun = ihw)
summary(res_rsv_IHW)

## Create a volcano plot
res_rsv_IHW_df <- as.data.frame(res_rsv_IHW)
summary(res_rsv_IHW_df$padj)
summary(res_rsv_IHW_df$log2FoldChange)
res_df <- res_rsv_IHW_df[res_rsv_IHW_df$padj < 0.8, ]
# Create custom color scheme for up and down regulated genes
# Set the base color as 'grey'
keyvals <- rep('grey50', nrow(res_df))
# Set the label/base name as 'Mid'
names(keyvals) <- rep('Not significant', nrow(res_df))
# Set keyvals for pval and FC - upregulated
keyvals[which(res_df$log2FoldChange > 0.01 & res_df$padj < 0.05)] <- 'red'
names(keyvals)[which(res_df$log2FoldChange > 0.01 & res_df$padj < 0.05)] <- 'Upregulated'
# Set keyvals for pval and FC - downregulated
keyvals[which(res_df$log2FoldChange < -0.01 & res_df$padj < 0.05)] <- 'forestgreen'
names(keyvals)[which(res_df$log2FoldChange < -0.01 & res_df$padj < 0.05)] <- 'Downregulated'
unique(names(keyvals))
unique(keyvals)
keyvals[1:20]
table(keyvals)
EnhancedVolcano(res_df,
                lab = rownames(res_df),
                xlim = c(-3, 3),
                ylim = c(0, 6),
                x = 'log2FoldChange',
                y = 'padj',
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                selectLab = c(""),
                pointSize = 4,
                xlab = bquote(~Log[2]~ '(fold change)'),
                ylab = bquote(~-Log[10]~ '(adjusted p-value)'),
                axisLabSize = 28,
                pCutoff = 0.05,
                FCcutoff = FALSE,
                colCustom = keyvals,
                hlineType = 'dashed',
                colAlpha = 0.8)

## Generate a heatmap
# Subset genes with padj<0.05
res_df_05 <- res_rsv_IHW_df[res_rsv_IHW_df$padj < 0.05, ]
# Create list of genes
gene_list <- rownames(res_df_05)
# Subset dataset
vsd_05 <- vsd[gene_list, ]
# Extract the expression matrix
mat <- assay(vsd_05)
mat <- as.matrix(mat)
# Scale for expression for each gene
base_mean <- rowMeans(mat) 
mat_norm <- t(apply(mat, 1, scale))
colnames(mat_norm) <- colnames(mat)
summary(mat_norm)
# Annotation for heatmap
table(rsv$Severity)
df = data.frame(group=rsv$Severity)
group.cols <- c("#2b83ba", "#ffffbf") #  blue and yellow
group.cols.assigned <- setNames(group.cols, unique(as.character(df$group)))
group.cols.assigned <- setNames(group.cols, c("severe", "mild"))
ha <- HeatmapAnnotation(df=df, col = list(group=group.cols.assigned))
Heatmap(mat_norm, name = "Expression values",
        col = colorRamp2(c(-1.5, 0, 1.5), c("forestgreen", "white", "red")),
        show_row_names = FALSE, 
        show_column_names = FALSE,
        row_dend_reorder = TRUE, 
        column_dend_reorder = TRUE,
        clustering_distance_rows = "pearson",
        clustering_distance_columns = "maximum",
        clustering_method_rows = "complete",
        clustering_method_columns = "complete",
        show_row_dend = FALSE,
        show_column_dend = TRUE,
        top_annotation = ha
)


### 5. FUNCTIONAL ENRICHMENT ANALYSIS 
# Create a SummarizedExperiment object
rowdata <- res_rsv_IHW_df
colnames(rowdata) <- c("baseMean", "FC", "FCse", "DESeq2.STAT", "PVAL", "ADJ.PVAL", "weight")
se <- SummarizedExperiment(assays = assay(rsv),
                           rowData = rowdata,
                           colData = colData(rsv))
# Convert gene symbol to Entrez ID
se <- idMap(se, org = "hsa", from = "SYMBOL", to = "ENTREZID")

# Create a GROUP variable
se$GROUP <- ifelse(se$Severity == "severe", 1, 0)

# Load gene sets from GO and KEGG databases
go.gs <- getGenesets(org="hsa", db="go")
kegg.gs <- getGenesets(org="hsa", db="kegg")

# Perform gene set enrichment analysis
gsea_res <- sbea(method = "gsea", 
                 se = se,
                 gs = go.gs,
                 perm = 1000,
                 alpha = 0.05)
gsRanking(gsea_res)

# Signaling pathway analysis
hsa.grn <- compileGRN(org = "hsa", db = "kegg")
spia_res <- nbea(method = "spia", 
                 se = se,
                 gs = kegg.gs,
                 grn = hsa.grn)
gsRanking(spia_res)


# Modular analysis using tmod package
res_tmod_df <- res_rsv_IHW_df
ord <- order(res_tmod_df$padj)
res_tmod_df$Symbol <- rownames(res_tmod_df)
#colnames(res_tmod_df) <- "Symbol"
res_tmod <- list()
res_tmod$tmod <- tmodCERNOtest(res_tmod_df$Symbol[ord], mset = "LI")
tmodSummary(res_tmod)
tmodPanelPlot(res_tmod, text.cex = 0.8)


