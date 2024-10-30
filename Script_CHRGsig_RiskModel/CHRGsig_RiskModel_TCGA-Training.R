# Version 3 SUmmary: previous versions could no tbe validated due to missing variables in the valdiation set. so the solution which will be applied here is that we will increase the sample size of the dataset by combining GBM+LGG then partitioning the data into training (70%) and test sets (30%)

# clear environment
rm(list=ls())

#### load libraries needed for the script to run ####

# package to perform data manipulation
# and visualization
library(tidyverse)
library(tidyr)
library(dplyr)

# for loading phenotype data
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)

# transpose
library(data.table)

# DEG analysis

#BiocManager::install("DESeq2")
library(DESeq2)

#BiocManager::install("edgeR")
library(edgeR)

library(limma)
library(VennDiagram) # to visualise DEGs

# cox analysis
library(ezcox)

#survival analysis and lasso
library(survival)
library(survminer)
library(glmnet)
library(timeROC)

# nomogram
library(rms)

# for %>% function
library(magrittr)

# for visulaisation
library(ggplot2)

####---- I. create GBM + LGG dataset a.k.a. glioma sets----####

### A. Get phenotype (TCGAbiolinks) data for GBM 

# GBM and LGG molecular subtypes
gbm.path.subtypes <- TCGAquery_subtype(tumor = "gbm")
lgg.path.subtypes <- TCGAquery_subtype(tumor = "lgg")

# Prep phenotype data
i <- sapply(gbm.path.subtypes, is.factor)
gbm.path.subtypes[i] <- lapply(gbm.path.subtypes[i], as.character)

i <- sapply(lgg.path.subtypes, is.factor)
lgg.path.subtypes[i] <- lapply(lgg.path.subtypes[i], as.character)

# make subset of pheno data to contain only needed columns
gbm.path.subtypes.cp <- gbm.path.subtypes[c("patient", "Histology", "Grade", "Age..years.at.diagnosis.", "Gender", "Survival..months.", "Karnofsky.Performance.Score", "IDH.status", "X1p.19q.codeletion", "IDH.codel.subtype", "MGMT.promoter.status", "TERT.promoter.status", "ATRX.status", "ESTIMATE.stromal.score", "ESTIMATE.immune.score", "ESTIMATE.combined.score", "Transcriptome.Subtype")]

lgg.path.subtypes.cp <- lgg.path.subtypes[c("patient", "Histology", "Grade", "Age..years.at.diagnosis.", "Gender", "Survival..months.", "Karnofsky.Performance.Score", "IDH.status", "X1p.19q.codeletion", "IDH.codel.subtype", "MGMT.promoter.status", "TERT.promoter.status", "ATRX.status", "ESTIMATE.stromal.score", "ESTIMATE.immune.score", "ESTIMATE.combined.score", "Transcriptome.Subtype")]

# create big dataframe of GBM + LGG phenotype data = glioma.pheno
gbm.path.subtypes.cp
lgg.path.subtypes.cp

glioma.pheno <- rbind(gbm.path.subtypes.cp, lgg.path.subtypes.cp)
dim(glioma.pheno)

# rename column names
glioma.pheno <- glioma.pheno %>% 
  dplyr::rename(Patient_ID = patient,
                Age_Years = Age..years.at.diagnosis.,
                KPS = Karnofsky.Performance.Score,
                Codel = X1p.19q.codeletion,
                IDH.codel.status = IDH.codel.subtype,
                MGMT_methy = MGMT.promoter.status,
                TERT.status = TERT.promoter.status)

save(glioma.pheno, file="OP/Pipeline OP5/glioma.surv.RData")

#### B. Get survival data (TCGA GDC) data for GBM and LGG ####
gbm.surv.file <- "/Users/tonijue/Documents/Copper in glioma/Copper in glioma/TCGA LGG 20220523/Scripts/Version 3_2023/Files/TCGA-GBM.survival.tsv"
lgg.surv.file <- "/Users/tonijue/Documents/Copper in glioma/Copper in glioma/TCGA LGG 20220523/Clinical Datasets/TCGA-LGG.survival.tsv"

gbm.surv <- read.table(gbm.surv.file, header = TRUE, stringsAsFactors = FALSE) # 649 rows/patients
lgg.surv <- read.table(lgg.surv.file, header = TRUE, stringsAsFactors = FALSE) # 533 rows/patients

# create big dataframe of GBM + LGG phenotype data = glioma.pheno
gbm.surv 
dim(gbm.surv)
# [1] 649   4

lgg.surv
dim(lgg.surv)
# [1] 533   4

glioma.surv <- rbind(gbm.surv, lgg.surv) # 1182 rows/patients
dim(glioma.surv)
# [1] 1182    4

save(glioma.surv, file ="OP/Pipeline OP5/glioma.surv.RData")

#### C. Get htseq_counts [unit: log2(count+1)] RNAseq data (GDC TCGA) with  for GBM and LGG, backtransform to raw counts, get unique gene IDs by aggregate. This is needed for differential gene expression analysis ####
gbm.log2.file <- "/Users/tonijue/Documents/Copper in glioma/Copper in glioma/TCGA LGG 20220523/Scripts/Version 3_2023/Files/TCGA-GBM.htseq_counts.tsv.gz" # htseq_counts [unit: log2(count+1)]
lgg.log2.file <- "/Users/tonijue/Documents/Copper in glioma/Copper in glioma/TCGA LGG 20220523/Scripts/Version 3_2023/Files/TCGA-LGG.htseq_counts.log2(count+1).tsv.gz" # htseq_counts [unit: log2(count+1)]

gbm.log2 <- read.table(gbm.log2.file, header = TRUE, stringsAsFactors = FALSE, row.names = "Ensembl_ID") 
# 60488 genes/ 173 patients
lgg.log2 <- read.table(lgg.log2.file, header = TRUE, stringsAsFactors = FALSE, row.names = "Ensembl_ID") 
# 60488 genes / 529 patients

# merge log2(count+1) data
glioma.log2 <- merge(gbm.log2, lgg.log2, by = "row.names", sort = FALSE) 
# 60488 genes / 703 patients

#### Continue C. here --- backtransform to get raw counts ####
glioma.log2 <- glioma.log2 %>% dplyr::rename(Ensembl_ID = Row.names)

glioma.log2 <- tibble::column_to_rownames(glioma.log2, var = "Ensembl_ID") # Apply rownames_to_column

n = ncol(glioma.log2)
tmp = 2^glioma.log2[, 1:n]
glioma.raw = tmp - 1

#> nrow(glioma.raw)
#[1] 60488
#> ncol(glioma.raw)
#[1] 702

# Subset RNAseq data that only contains those with survival data
# tumour.pt.ID <- list(LGG.subtypes$patient) # makes a list with one element.
tumour.pt.ID <- as.factor(glioma.surv$X_PATIENT) # makes list with 1182 elements.
tumour.sample.ID <- as.factor(glioma.surv$sample)

# replace - with . to match the colnames in the counts dataframe.
tumour.sample.ID <- str_replace_all(tumour.sample.ID, "-", ".")

# get intersect of the patient IDs from the glioma.raw dataframe and the tumour patient sample IDs
tmp <- colnames(glioma.raw)
tmp2 <- intersect(tmp, tumour.sample.ID)

# subset patient IDs that were only in the survival dataframe
#glioma.raw.cp <- glioma.raw %>% select(c(all_of(tmp2)))
#glioma.log2.cp <- glioma.log2 %>% select(c(all_of(tmp2)))

glioma.raw.cp <- glioma.raw[,tmp2]
glioma.log2.cp <- glioma.log2[,tmp2]

#> ncol(glioma.log2)
#[1] 702
#> ncol(glioma.log2.cp)
#[1] 691

# Convert Ensembl IDs for both log2 and raw datasets

EnsemblID_to_Genes <- function(counts_df, genecode_map){
  # convert rownames in counts dataframe into a column
  counts_df$Ensembl_ID <- rownames(counts_df)
  counts_df <- counts_df %>% dplyr::select(Ensembl_ID, everything())
  
  # match ensembl ID of genecode dataframe to pt.raw.counts dataframe
  counts_df <- merge(counts_df, genecode_map, by.x = "Ensembl_ID", by.y = "id")
  counts_df <- counts_df %>% dplyr::select(gene, everything())
  counts_df <- counts_df[ ,-2]
  
  return(counts_df)
  }

genecode.file <- "/Users/tonijue/Documents/Copper in glioma/Copper in glioma/TCGA LGG 20220523/Scripts/Version 3_2023/Version 3_2023_rebuttal/IP/TCGA.gencode.v22.annotation.gene.probeMap"

genecode <- read.table(genecode.file, header = TRUE, sep = "\t")
genecode.tmp <- genecode[, c("id", "gene")]

#convert ensemblID to genes
glioma.raw.cp2 <- EnsemblID_to_Genes(glioma.raw.cp, genecode.tmp)
glioma.log2.cp2 <- EnsemblID_to_Genes(glioma.log2.cp, genecode.tmp)

save(glioma.raw.cp2, file="OP/Pipeline OP5/glioma.raw.cp2.RData")
save(glioma.log2.cp2, file="OP/Pipeline OP5/glioma.log2.cp2.RData")

# get unique gene names log2(counts+1) df
get_unique_genes <- function(counts_df){
  dt <- as.data.table(counts_df) 
  setkey(dt, gene)
  
  library(future.apply)
  result <- dt[, future_lapply(.SD, mean), by = gene]
  
  return(as.data.frame(dt))
  
 # counts_df <- aggregate(counts_df, by=list(counts_df$gene), FUN=mean)
  #counts_df$gene <- counts_df$Group.1
  #counts_df <- counts_df[ ,-1]
  
  #return(counts_df)
}

#glioma.raw.unique <- get_unique_genes(glioma.raw.cp2)
#> nrow(glioma.raw.unique)
#[1] 58387
dt <- as.data.table(glioma.raw.cp2) 
setkey(dt, gene)

library(future.apply)
result <- dt[, future_lapply(.SD, mean), by = "gene"]

glioma.raw.unique <- as.data.frame(result)

# get unique gene names log2(countts+1) df
# gene IDs must be in rownames which will be used later for the unicox analysis
dt <- as.data.table(glioma.log2.cp2) 
setkey(dt, gene)

result <- dt[, future_lapply(.SD, mean), by = "gene"]

glioma.log2.unique <- as.data.frame(result)

#glioma.log2.unique <- get_unique_genes(glioma.log2.cp2)
#glioma.log2.unique <- tibble::column_to_rownames(glioma.log2.unique, "gene")
#> ncol(glioma.log2.unique)
#[1] 691
#> nrow(glioma.log2.unique)
#[1] 58387

save(glioma.raw.unique, file ="OP/Pipeline OP5/glioma.raw.unique.RData")
save(glioma.log2.unique, file ="OP/Pipeline OP5/glioma.log2.unique.RData")

write.table(glioma.log2.unique, "/Users/tonijue/Documents/Copper in glioma/Copper in glioma/TCGA LGG 20220523/Scripts/Version 3_2023/glioma.log2.uniqueGenes.txt", sep="\t", row.names=TRUE, quote=FALSE)

write.table(glioma.raw.unique, "/Users/tonijue/Documents/Copper in glioma/Copper in glioma/TCGA LGG 20220523/Scripts/Version 3_2023/glioma.raw.uniqueGenes.txt", sep="\t", row.names=FALSE, quote=FALSE)

####---- II. Partition data into Training (70%) and Test (30%) Sets----####
# Reference: https://www.statology.org/train-test-split-r/. Use BaseR method

# transpose dataframe so that patient IDs are as rows and gene names are as columns
rownames(glioma.raw.unique) <- NULL
glioma.raw.unique <- column_to_rownames(glioma.raw.unique, var = "gene")

glioma.raw.unique.T <- transpose(glioma.raw.unique)

#redefine row and column names
rownames(glioma.raw.unique.T) <- colnames(glioma.raw.unique)
colnames(glioma.raw.unique.T) <- rownames(glioma.raw.unique)

#make this example reproducible
set.seed(100)

#Use 70% of dataset as training set and remaining 30% as testing set
sample <- sample(c(TRUE, FALSE), nrow(glioma.raw.unique.T), replace=TRUE, prob=c(0.7,0.3))
save(sample, file="OP/Pipeline OP5/Pipeline5_TCGA_sampleID_partition.RData")

glioma.train.raw  <- glioma.raw.unique.T[sample, ] # 478 samples
glioma.test.raw   <- glioma.raw.unique.T[!sample, ] # 213 samples

# get sample IDs of train set and test set
glioma.train.IDs <- rownames(glioma.train.raw) 
glioma.test.IDs <- rownames(glioma.test.raw)

save(glioma.train.IDs, file="OP/Pipeline OP5/Pipeline5_TCGA_glioma.train.IDs.RData")
save(glioma.test.IDs, file="OP/Pipeline OP5/Pipeline5_TCGA_glioma.test.IDs.RData")

# create train and test set for log2 counts dataframe based on IDs taken above
glioma.log2.unique <- column_to_rownames(glioma.log2.unique, var="gene")

tmp <- transpose(glioma.log2.unique)

colnames(tmp) <- rownames(glioma.log2.unique)
rownames(tmp) <- colnames(glioma.log2.unique)

glioma.log2.unique <- tmp

glioma.train.log2 <- glioma.log2.unique[glioma.train.IDs,]

glioma.test.log2 <- glioma.log2.unique[glioma.test.IDs,]

save(glioma.train.log2, file = "OP/Pipeline OP5/TCGA_glioma.train.log2counts.RData")
save(glioma.test.log2, file = "OP/Pipeline OP5/TCGA_glioma.test.log2counts.RData")

####---- III. Prepare normal GTEX samples ----####
GTEX.file <- "/Users/tonijue/Documents/Copper in glioma/Copper in glioma/TCGA LGG 20220523/Scripts/Version 3_2023/Files/gtex_gene_expected_count_log2(count+1).gz"
gtex.log2 <- read.table(GTEX.file, header = TRUE)
dim(gtex.log2)
# [1] 60498  7846

gtex.log2 <- column_to_rownames(gtex.log2, var="sample")

GTEX.pheno.file <- "/Users/tonijue/Documents/Copper in glioma/Copper in glioma/TCGA LGG 20220523/Scripts/Version 3_2023/Files/GTEX_phenotype.gz"
gtex.pheno <- read.table(GTEX.pheno.file, header = TRUE, sep = "\t")
dim(gtex.pheno)
# [1] 9783    6

# A. Create subset dataframe to contain only brain cortex samples

# first filter out other rows not in brain
gtex.brain <- gtex.pheno[gtex.pheno$X_primary_site == "Brain", ]
gtex.cortex <- gtex.brain[gtex.brain$body_site_detail..SMTSD. == "Brain - Cortex", ]

# make list of column IDs
gtex.cortex$Sample <- str_replace_all(gtex.cortex$Sample, "-", ".") # replace - with .
gtex.sample.ID <- as.character(gtex.cortex$Sample) # makes list with 133 elements.
#gtex.sample.ID <- append(gtex.sample.ID, "sample", 0)

# subset patient IDs that were only in the subtypes dataframe

# optional check: checking column names to keep are present in df ####
cols_to_keep <- gtex.sample.ID

columns_present <- cols_to_keep %in% names(gtex.log2)
# Create a named logical vector showing which columns are present
result <- setNames(columns_present, cols_to_keep)

# Print the result
print(result)

# Optionally, print a summary message
cat("\nSummary:\n")
cat("Columns present:", sum(columns_present), "\n") # 105
cat("Columns missing:", sum(!columns_present), "\n") # 28

# If you want to see which columns are missing
missing_columns <- cols_to_keep[!columns_present]
if (length(missing_columns) > 0) {
  cat("\nMissing columns:\n")
  print(missing_columns)
}

# remove missing column from gtex.sample.ID
cols_to_keep <- gtex.sample.ID[columns_present]

# subset gtex.log2 df to keep only cortex samples present
cortex.log2 <- gtex.log2[ ,cols_to_keep] 

# B. Back transform log2(counts+1). Ref: https://assets.researchsquare.com/files/pex-751/v1/1892ff4d-9dd8-42ef-8f64-b103bb6fdb63.pdf

# assign gene column as rowname
cortex.raw = (2^cortex.log2)-1

# C. Convert Ensmebl Ids to Gene names
gtex.genecode.file <- "/Users/tonijue/Documents/Copper in glioma/GTEX/probeMap_gencode.v23.annotation.gene.probemap"
gtex.genecode <- read.table(gtex.genecode.file, header = TRUE, sep = "\t")
gtex.genecode.tmp <- gtex.genecode[, 1:2]

cortex.raw.cp <- EnsemblID_to_Genes(cortex.raw, gtex.genecode.tmp)
#> ncol(cortex.raw.cp)
#[1] 106
#> nrow(cortex.raw.cp)
#[1] 60498

# D. get unique gene names
cortex.raw.unique <- get_unique_genes(cortex.raw.cp)
#> ncol(cortex.raw.unique)
#[1] 106
#> nrow(cortex.raw.unique)
#[1] 58581

dt <- as.data.table(cortex.raw.cp) 
setkey(dt, gene)

result <- dt[, future_lapply(.SD, mean), by = "gene"]

cortex.raw.unique <- as.data.frame(result)

save(cortex.raw.unique, file="OP/Pipeline OP5/cortex.raw.unique.RData")

write.table(cortex.raw.unique, "/Users/tonijue/Documents/Copper in glioma/Copper in glioma/TCGA LGG 20220523/Scripts/Version 3_2023/cortex.raw.uniqueGenes.txt", sep="\t", row.names=FALSE, quote=FALSE)

####---- IV. create big data frame that has all GTEX and TCGA raw counts data AND a metadata ----####

### Create metadata for as input differential gene expression analysis

# 1. Create sorted list off GTEX and TCGA train set from the counts data frames
cortex.sample.ID <- colnames(cortex.raw.unique)
cortex.sample.ID <- cortex.sample.ID[-1] #this removes the first element "gene"
cortex.sample.ID <- sort(cortex.sample.ID)

glioma.train.ID <- sort(glioma.train.IDs)

# 2. Make a list that says Normal or Tumor
Normal = replicate(length(cortex.sample.ID), "Normal")
Batch1 = replicate(length(cortex.sample.ID), "GTEX")

Tumor = replicate(length(glioma.train.ID), "Tumor")
Batch2 = replicate(length(glioma.train.ID), "TCGA")

# 3. Create metadata dataframe for combined GTEX and TCGA datasets
normal.cortex.meta = do.call(rbind, Map(data.frame, Sample_ID = cortex.sample.ID, Condition = Normal, Batch = Batch1))

tumor.glioma.meta = do.call(rbind, Map(data.frame, Sample_ID = glioma.train.ID, Condition = Tumor, Batch = Batch2))

# concatenate the two dataframes
cortex.glioma.meta <- rbind(normal.cortex.meta, tumor.glioma.meta)

# remove display of rownames
rownames(cortex.glioma.meta) <- NULL

save(cortex.glioma.meta, file = "OP/Pipeline OP5/1 cortex.glioma.meta-DEGanalysis.Rdata")

# 4. Create biq Raw counts dataframe that contains normal and tumor samples
# REf: https://towardsdatascience.com/merge-data-frames-in-python-r-725c0f874147

# transpose dataframe so that patient IDs are as columns and gene names are as rows this is to prepare dataframe for binding.
glioma.train.raw.T <- transpose(glioma.train.raw)

#redefine row and column names
rownames(glioma.train.raw.T) <- colnames(glioma.train.raw)
colnames(glioma.train.raw.T) <- rownames(glioma.train.raw)

cortex.raw.unique <- column_to_rownames(cortex.raw.unique, var="gene")

# merge raw counts dataframes. 
cortex.glioma.raw <- merge(cortex.raw.unique, glioma.train.raw.T, by=0)

# make sure that the first column (gene names) have been converted into rownames so that it dooesn't get affected by the following code below.
cortex.glioma.raw <- column_to_rownames(cortex.glioma.raw, var="Row.names")

i <- sapply(cortex.glioma.raw, is.character) # identifies which column is character type
cortex.glioma.raw[i] <- lapply(cortex.glioma.raw[i], as.numeric) # converts those identified as character type into numeric.

save(cortex.glioma.raw, file = "OP/Pipeline OP5/2 cortex.glioma.raw-DEGanalysis.Rdata")

####---- V. Differential analysis ----####

genes_for_DEG <- cortex.glioma.raw

###DESeq2####

# assign gene column as rowname
expr <- genes_for_DEG[rowMeans(genes_for_DEG)>1, ]
expr <- floor(expr) 

expr <- expr[ ,sort(names(expr))]

expr <- as.data.frame(expr)

colData <- cortex.glioma.meta

dds <- DESeqDataSetFromMatrix(countData = round(expr),
                              colData = colData,
                              design = ~ Condition)

dds2 <- DESeq(dds)

resultsNames(dds2)

res <- results(dds2, name = "Condition_Tumor_vs_Normal")

res <- res[order(res$padj), ]

head(res)

resdata <- merge(as.data.frame(res),
                 as.data.frame(counts(dds2,normalized = TRUE)),by = "row.names",sort = FALSE)

save(resdata, file = 'OP/Pipeline OP5/3 TCGA_train-DESeq2data.Rdata')

DEG = as.data.frame(res)

DEG_DESeq2 = na.omit(DEG) 

logFC_cutoff <- 1

DEG_DESeq2$change = as.factor(
  ifelse(DEG_DESeq2$padj < 0.05 & abs(DEG_DESeq2$log2FoldChange) > logFC_cutoff,
         ifelse(DEG_DESeq2$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT')
)

table(DEG_DESeq2$change)
#  DOWN   NOT    UP 
#  6489 16252 10145 

save(DEG_DESeq2, file="OP/Pipeline OP5/4a DEG_DESeq2.RData")
#write.csv(DEG_DESeq2, file = "/Users/tonijue/Documents/Copper in glioma/Copper in glioma/TCGA LGG 20220523/Scripts/Version 3_2023/DEG_DESeq2_20233110.csv", row.names = TRUE)

#> table(DEG_DESeq2$change) 

#DOWN   NOT    UP 
#6404 16208 10274 # original in submission

####edgeR####
group_list <- cortex.glioma.meta[ ,2]

exprSet <- genes_for_DEG[ ,sort(names(genes_for_DEG))]

dge <- DGEList(counts = exprSet, group=group_list)
design <- model.matrix(~0 + group_list)
rownames(design) <- colnames(dge)
colnames(design) <- levels(group_list)
design

keep_gene <- rowSums(cpm(dge) > 1 ) >= 2
table(keep_gene)
# FALSE  TRUE 
# 32394 25399 

dge <- dge[ keep_gene, , keep.lib.sizes = FALSE]
dge$samples$lib.size <- colSums(dge$counts)
dge <- calcNormFactors(dge)
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

fit <- glmFit(dge, design)
fit2 <- glmLRT(fit, contrast = c(-1,1))
DEG <- topTags(fit2, n = nrow(exp))
DEG <- as.data.frame(DEG)
DEG$change = as.factor(
  ifelse(DEG$FDR< 0.05 & abs(DEG$logFC) > logFC_cutoff,
         ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)

table(DEG$change)

#> table(DEG$change)

#DOWN   NOT    UP 
#4941 13090  7368 # rerun and original in submission are same - GOOD!

DEG_edgeR <- DEG

save(DEG_edgeR, file = "OP/Pipeline OP5/4b DEG_edgeR.RData")

#write.csv(DEG_edgeR, file = "/Users/tonijue/Documents/Copper in glioma/Copper in glioma/TCGA LGG 20220523/Scripts/Version 3_2023/DEG_edgeR_20233110csv", row.names = TRUE)

####limma ####
exprSet

d0 <- DGEList(exprSet, group=group_list)
d0 <- calcNormFactors(d0)

cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d)

#keep filter selection to default based on tutorial.
#keep_gene <- rowSums(cpm(d0) > 1 ) >= 5
#table(keep_gene)
#d <- d0[keep_gene, , keep.lib.sizes = FALSE]
#dim(d)

#{
#  design <- model.matrix(~0 + factor(group_list) )
#  colnames(design) = levels(factor(group_list) )
#  rownames(design) = colnames(exprSet)
#}

#design

group_list <- factor(group_list)
counts <- d #exprSet 
design <- model.matrix(~0 + group_list)
rownames(design) <- colnames(counts)
colnames(design) <- levels(group_list)

head(design)
tail(design)

v <- voom(d, design)#, plot = T) # checks mean variance distribution. if need to get a "good plot" I suggest using the code in lines 409-413 to fix the curve.
# 20233110: This part of the code is throwing off an error ", plot = T)". Error is "Error in plot.new() : figure margins too large". Silenced it for the mean time. see if it chnages the result. No time to tease out the issue.

fit <- lmFit(v, design = design)
head(coef(fit))

contrast.matrix <- makeContrasts("Tumor-Normal", levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

DEG_limma <- topTable(fit2, coef = 1, n = Inf)
DEG_limma <- na.omit(DEG_limma)

DEG_limma$change <- as.factor(
  ifelse(DEG_limma$adj.P.Val < 0.05 & abs(DEG_limma$logFC) > logFC_cutoff,
         ifelse(DEG_limma$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)

table(DEG_limma$change)

#> table(DEG_limma$change)

#DOWN   NOT    UP 
#6932 15359  5403 # rerun and original in submission are same - GOOD!

DEG_limmavoom <- DEG_limma

save(DEG_limmavoom, file = "OP/Pipeline OP5/4c DEG_limavoom.Rdata")
#write.csv(DEG_limmavoom, file = "/Users/tonijue/Documents/Copper in glioma/Copper in glioma/TCGA LGG 20220523/Scripts/Version 3_2023/DEG_limmavoom_20233110.csv", row.names = TRUE)

logCPM <- cpm(dge, prior.count = 3, log = TRUE)

save(logCPM, file = "OP/Pipeline OP5/4d TCGA_logCPM_limma-voom_20233110.Rdata")

save(DEG_DESeq2, DEG_edgeR, DEG_limmavoom, group_list, file = "OP/Pipeline OP5/5  DEGs_grouplist_20233110.Rdata")


#### DEG gene####
load("/Users/tonijue/Documents/Copper in glioma/Copper in glioma/TCGA LGG 20220523/Scripts/Version 3_2023/6  DEGs_grouplist_20233110.Rdata") # start point to continue analysis so that you don't have to run line 1 to 342 again.

D_edge<-DEG_edgeR[which(DEG_edgeR$change!='NOT'), ]
D_limma<-DEG_limmavoom[which(DEG_limmavoom$change!='NOT'), ]
D_DES<-DEG_DESeq2[which(DEG_DESeq2$change!='NOT'), ]

edgeR <- rownames(D_edge)
dim(D_edge)
# [1] 12309     6

limma <- rownames(D_limma)
dim(D_limma)
# [1] 12335     7

DESeq2 <- rownames(D_DES)
dim(D_DES)
# [1] 16634     7

DEGintersect <- list(edgeR, DESeq2, limma)

venn.diagram(DEGintersect,
             category.names = c('edgeR (12,309)' , 'limma (12,335)' , 'DESeq2 (16,678)'),
             filename = "OP/Pipeline OP5/Images/TCGA DEG Venn.tiff",
             col = "black",
             fill = c("dodgerblue", "goldenrod1", "darkorange1"),
             alpha = 0.5,
             cex = 2.4,
             cat.col = "black",
             cat.cex = 1.8,
             cat.fontface = "bold",
             margin = 0.05,
             main = "DEG analysis",
             main.cex = 2.4,
             
             # New parameters for TIFF output and resolution
             imagetype = "tiff",
             resolution = 300,
             compression = "lzw",
             height = 3000,  # in pixels
             width = 3000,   # in pixels
             units = "px"
)

D_DES$gene <- rownames(D_DES)
D_edge$gene <- rownames(D_edge)
D_limma$gene <- rownames(D_limma)
e_D <- merge(x = D_edge,y = D_DES,by = 'gene')
ed_l <- merge(x = e_D,y = D_limma,by = 'gene')

DEGlist <- ed_l[ ,c(1, 15:21)]
dim(DEGlist)
# [1] 9756    8

#fix(DEGlist) # Opens up R editor. not sure what this is for

save(DEGlist,file = "OP/Pipeline OP5/6 TCGA DEG_list_9756.Rdata")

load("OP/Pipeline OP5/6 TCGA DEG_list_9756.Rdata")

#write.csv(DEGlist, file = "/Users/tonijue/Documents/Copper in glioma/Copper in glioma/TCGA LGG 20220523/Scripts/Version 3_2023/9 TCGA DEG_list_9760.csv", row.names = FALSE)

####---- VI. Univariate Cox Analysis ----####
load("~/Documents/Copper in glioma/Copper in glioma/TCGA LGG 20220523/Scripts/Version 3_2023/Version 3_2023_rebuttal/OP/clean_TCGA.WHO-CNS5_v2.RData")


# Univariate Cox Analysis must use log2(count+1) data

### A. First select differentially expressed copper genes from DEGlist
copper.genelist <- c('ABCB6', 'ANKRD9', 'SLC31A1', 'SLC31A2', 'PRND', 'CCDC22', 'APP', 'ARF1', 'MT2A', 'ATOX1', 'ATP7A', 'ATP7B', 'PRNP', 'SCO1', 'COX19', 'SCO2', 'CYP1A1', 'DAXX', 'BACE1', 'AOC1', 'MT1DP', 'HSF1', 'AQP1', 'AQP2', 'MT1A', 'MT1B', 'MT1E', 'MT1F', 'MT1G', 'MT1H', 'MT1M', 'MT1X', 'MT3', 'NFE2L2', 'MT1HL1', 'SNCA', 'MAP1LC3A', 'MT4', 'BECN1', 'COMMD1', 'XIAP', 'CUTC', 'STEAP2', 'STEAP3', 'STEAP4', 'SLC11A2', 'COX17', 'CP', 'FKBP4', 'HEPHL1', 'MMGT1', 'HEPH', 'PARK7', 'AANAT', 'IL1A', 'LCAT', 'LOXL2', 'MT-CO1', 'PAM', 'ATP5F1D', 'SOD1', 'SOD3', 'SORD', 'TFRC', 'CDK1', 'MOXD2P', 'MTCO2P12', 'COX11', 'LACC1', 'DBH', 'DCT', 'ALB', 'F5', 'F8', 'OR5AR1', 'ADNP', 'ATP13A2', 'MOXD1', 'GPC1', 'ANG', 'SUMF1', 'AOC2', 'SNAI3', 'APOA4', 'COA6', 'LOX', 'LOXL1', 'MT-CO2', 'ACR', 'P2RX4', 'CUTA', 'HAMP', 'S100A5', 'S100A12', 'S100A13', 'SNCB', 'SNCG', 'TP53', 'TYR', 'LOXL4', 'LOXL3', 'AOC3', 'RNF7', 'CCS', 'AP1S1', 'AP1B1', 'TMPRSS6', 'SPATA5', 'COG2', 'ATP6V0A2', 'ATP6AP1', 'ADAM10', 'AKT1', 'MTF2', 'FOXO1', 'FOXO3', 'STEAP1', 'GSK3B', 'APC', 'JUN', 'MAPT', 'MDM2', 'MT1JP', 'MT1L', 'MTF1', 'PIK3CA', 'XAF1', 'PTEN', 'CCND1', 'SP1', 'ADAM17', 'CASP3', 'ADAM9')

# getting the subset DataFrame after checking values if belonging to vector
copper.DEG <- DEGlist[DEGlist$gene %in% copper.genelist, ]

copper.DEGlist <- c(copper.DEG$gene)

### B. create a subset of log2.counts that contains only the copper.DEGlist
glioma.train.log2.copper <- glioma.train.log2[,copper.DEGlist]
dim(glioma.train.log2.copper)
# [1] 478  50

glioma.test.log2.copper <- glioma.test.log2[,copper.DEGlist]
dim(glioma.test.log2.copper)
# [1] 213  50

glioma.log2.copper <- glioma.log2.unique[,copper.DEGlist]
dim(glioma.log2.copper)
# [1] 691  50


### C. Need to recreate the coxdata table I had from python which included the expression matrix plus the survival times from the survival file.
clean.WHO.CNS5_v2

# Function to extract the desired part of the string
extract_id <- function(x) {
  parts <- strsplit(x, "\\.")
  paste(sapply(parts, function(p) paste(p[1:3], collapse = ".")), collapse = "")
}

# Apply the function to rownames and create a new column
glioma.train.log2.copper$Patient_ID <- NA
glioma.train.log2.copper$Patient_ID <- sapply(rownames(glioma.train.log2.copper), extract_id)
length(unique(glioma.train.log2.copper$Patient_ID))
# [1] 466

glioma.test.log2.copper$Patient_ID <- NA
glioma.test.log2.copper$Patient_ID <- sapply(rownames(glioma.test.log2.copper), extract_id)
length(unique(glioma.test.log2.copper$Patient_ID))
# [1] 206

# replace - with . in surv file to match the colnames in the counts dataframe.
tmp <- str_replace_all(glioma.surv$sample, "-", ".")
glioma.surv$sample <- tmp

tmp <- str_replace_all(glioma.surv$X_PATIENT, "-", ".")
glioma.surv$X_PATIENT <- tmp

glioma.surv <- column_to_rownames(glioma.surv, var="sample")

# merge data into dataframe with rnaseq data of copper DEGs and OS/OS.time
glioma.train.log2.coxdata <- merge(glioma.train.log2.copper, glioma.surv, by = 0)
glioma.test.log2.coxdata <- merge(glioma.test.log2.copper, glioma.surv, by = 0)
glioma.log2.coxdata <- merge(glioma.log2.copper, glioma.surv, by = 0)

## get unique patient IDs. note that some patients have duplicate samples
glioma.train.log2.coxdata.unique <- aggregate(glioma.train.log2.coxdata, by=list(glioma.train.log2.coxdata$X_PATIENT), FUN=mean)
glioma.train.log2.coxdata.unique <- glioma.train.log2.coxdata.unique[ ,-c(2, 53, 55)]
glioma.train.log2.coxdata.unique <- glioma.train.log2.coxdata.unique %>% dplyr::rename(Patient_ID = Group.1)

glioma.test.log2.coxdata.unique <- aggregate(glioma.test.log2.coxdata, by=list(glioma.test.log2.coxdata$X_PATIENT), FUN=mean)
glioma.test.log2.coxdata.unique <- glioma.test.log2.coxdata.unique[ ,-c(2, 53, 55)]
glioma.test.log2.coxdata.unique <- glioma.test.log2.coxdata.unique %>% dplyr::rename(Patient_ID = Group.1)

glioma.log2.coxdata.unique <- aggregate(glioma.log2.coxdata, by=list(glioma.log2.coxdata$X_PATIENT), FUN=mean)
glioma.log2.coxdata.unique <- glioma.log2.coxdata.unique[ ,-c(2, 54)]
glioma.log2.coxdata.unique <- glioma.log2.coxdata.unique %>% dplyr::rename(Patient_ID = Group.1)

#convert first column to rowname
glioma.train.log2.coxdata.unique <- column_to_rownames(glioma.train.log2.coxdata.unique, var="Patient_ID")

glioma.test.log2.coxdata.unique <- column_to_rownames(glioma.test.log2.coxdata.unique, var="Patient_ID")

glioma.log2.coxdata.unique <- column_to_rownames(glioma.log2.coxdata.unique, var="Patient_ID")

save(glioma.train.log2.coxdata.unique, file = "OP/Pipeline OP5/7 TCGA_train_log2_coxdata_limma-voom.Rdata")
save(glioma.test.log2.coxdata.unique, file = "OP/Pipeline OP5/8 TCGA_test_log2_coxdata_limma-voom.Rdata")
save(glioma.log2.coxdata.unique, file = "OP/Pipeline OP5/9 TCGA_log2_coxdata_limma-voom.Rdata")

### D. First Univariate Cox Analysis
load("OP/Pipeline OP5/7 TCGA_train_log2_coxdata_limma-voom.Rdata")
load("OP/Pipeline OP5/8 TCGA_test_log2_coxdata_limma-voom.Rdata")
load("OP/Pipeline OP5/9 TCGA_log2_coxdata_limma-voom.Rdata")


# Univariate COX Analysis
coxdata <- glioma.train.log2.coxdata.unique # TCGA train set

#coxdata <- glioma.test.log2.coxdata.unique # TCGA test set

#coxdata <- glioma.IDHwt.log2.copper # IDHwt set

#coxdata <- glioma.IDHm.log2.copper # IDHm set

# to check distibution of gene expression (Optional) - they should have a normal distribution
ggplot(coxdata, aes(x = MT1X, color = OS)) + geom_histogram(color = "black", fill = "white", bins = 30) 
ggplot(coxdata, aes(x = SLC31A1, color = OS)) + geom_histogram(color = "black", fill = "white", bins = 30) 
ggplot(coxdata, aes(x = MT1M, color = OS)) + geom_histogram(color = "black", fill = "white", bins = 30) 

cox_res <- ezcox(coxdata, time = "OS.time", status = "OS", covariates = copper.DEGlist, global_method = c("likelihood", "wald", "logrank"))

unicox <- as.data.frame(cox_res)

save(unicox, file="OP/Pipeline OP5/Train OP/unicox.RData")
write.csv(unicox, file = "OP/Pipeline OP5/Train OP/CSV/TCGA_train_unicox.csv", row.names=FALSE)

load("OP/Pipeline OP5/Train OP/unicox.RData")
####---- VI. Constructing and validating risk-score system ----####

######contructing risk-score system, LASSO model#######

unicox_dif <- unicox[which(unicox$global.pval < 0.05), ] # do not change
uni_gene <- as.character(unicox_dif$Variable) # do not chnage

# change line coxdata to test set for validation.
irg_expr <- coxdata[ ,c(1:50)] # must include all gene columns from coxdata data frame. 
# fix(coxdata) # opens R editor. not sure what for.

x <- as.matrix(irg_expr[ ,uni_gene]) # uni_gene is taken from the train dataset. 
y <- data.matrix(Surv(coxdata$OS.time,coxdata$OS)) 

###################### keep this unchanged from train set, line 569 to 601 ###################### 
# Coefficient profiles in the LASSO regression model.
fit0 <- glmnet(x, y, family = "cox", alpha = 1, nlambda = 1000)
plot(fit0)
plot(fit0, xvar="lambda", label=TRUE)

#### Save LASSO plots ####

# open a TIFF device to save the plot####
tiff(filename = "OP/Pipeline OP5/Train OP/Images/TCGA_train_log2counts - lasso_path_plot2.tiff",
     width = 10, height = 10, units = "in", res = 300)

# Now, create your plot
plot(fit0) # LASSO (Least Absolute Shrinkage and Selection Operator) regularization path plot

# Close the device
dev.off()

# open a TIFF device to save the plot####
tiff(filename = "OP/Pipeline OP5/Train OP/Images/TCGA_train_log2counts - lasso_path_plot-lambda2.tiff", width = 10, height = 10, units = "in", res = 300)

# Now, create your plot
plot(fit0, xvar="lambda", label=TRUE) # LASSO (Least Absolute Shrinkage and Selection Operator) regularization path plot

# Close the device
dev.off()


# Cross-validation for tuning parameter screening in the LASSO regression model. ####
set.seed(100)

cv.fit <- cv.glmnet(x, y,
                    family="cox",
                    maxit = 1000000,
                    alpha=1)
print(cv.fit)
# Call:  cv.glmnet(x = x, y = y, family = "cox", maxit = 1e+06, alpha = 1) 

# Measure: Partial Likelihood Deviance 

#      Lambda Index Measure     SE Nonzero
# min 0.03357    26   10.64 0.2180      19
# 1se 0.10253    14   10.82 0.2325      11

plot(cv.fit)

# open a TIFF device to save the plot ####
tiff(filename = "OP/Pipeline OP5/Train OP/Images/TCGA_train_log2counts - partial-likelihood-deviance2.tiff", width = 10, height = 10, units = "in", res = 300)

# Now, create your plot
plot(cv.fit)

# Close the device
dev.off()

save(cv.fit, file ="OP/Pipeline OP5/Train OP/TCGA_train_log2counts - cv-fit.RData")

# LASSO_gene table
fit <- glmnet(x, y, alpha = 1, family='cox',lambda=cv.fit$lambda.min)
plot(fit)
coef(fit)

Coefficients <- coef(fit, s = cv.fit$lambda.min)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index]

Active.Index
# [1]  4  5  7 11 15 16 17 18 19 20 22 24 27 30 31 37 40 43 44

Active.Coefficients
# [1] -0.007409754  0.155383500  0.073934367  0.217785923  0.032618756  0.038602690 -0.167637522  0.070728710
# [9]  0.099598934  0.065588082  0.029122490 -0.062421631  0.006042011  0.087743960 -0.050878675 -0.001973116
# [17]  0.056465834  0.103646340 -0.395985411

lasso_gene <- row.names(Coefficients)[Active.Index]
lasso_min <- data.frame(Active.Index,Active.Coefficients,lasso_gene)

save(lasso_min,file = 'OP/Pipeline OP5/Train OP/TCGA_train_log2counts - lasso_min.RData')

save(cv.fit,fit,lasso_gene,file = 'OP/Pipeline OP5/Train OP/TCGA_train_log2counts - lasso_model_min.Rdata')


##### validation of LASSO model ###### 
# NOTE: Should be run on a new dataset. new x and y should be assigned for new dataset.
# to validate with test set with new gene expression which is the x variable 
lasso.prob <- predict(cv.fit, 
                      newx = x,
                      s = c(cv.fit$lambda.min,cv.fit$lambda.1se)
)
re <- cbind(y, lasso.prob)
head(re)

##### timeROC curve #####

coxdata$Risk.Score <- as.numeric(lasso.prob[ ,1])

with(coxdata,
     ROC <<- timeROC(T = OS.time,
                     delta = OS,
                     marker = Risk.Score,
                     cause = 1,
                     weighting = "marginal",
                     times = c(365, 1095, 1825),
                     ROC = TRUE,
                     iid = TRUE)
)

ROC$AUC
#     t=365    t=1095    t=1825 
# 0.8875994 0.9344898 0.8769280 

set.seed(100)
confint(ROC)
# $CI_AUC
# 2.5% 97.5%
# t=365  85.36 92.16
# t=1095 90.16 96.74
# t=1825 82.16 93.23

# $CB_AUC
# 2.5% 97.5%
# t=365  84.72 92.80
# t=1095 89.54 97.36
# t=1825 81.11 94.28

# $C.alpha
# 95% 
# `2.330898 


{
  auc_365 = ROC$AUC[[1]]
  auc_1095 = ROC$AUC[[2]]
  auc_1825 = ROC$AUC[[3]]
  dat <- data.frame(tpr365 = ROC$TP[,1],
                    fpr365 = ROC$FP[,1],
                    tpr1095 = ROC$TP[,2],
                    fpr1095 = ROC$FP[,2],
                    tpr1825 = ROC$TP[,3],
                    fpr1825 = ROC$FP[,3])
  library(ggplot2)
  ggplot() +
    geom_line(data = dat,aes(x = fpr365, y = tpr365),color = "blue", size=0.75) +
    geom_line(data = dat,aes(x = fpr1095, y = tpr1095),color = "red", size=0.75)+
    geom_line(data = dat,aes(x = fpr1825, y = tpr1825),color = "orange", size=0.75)+
    geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
    theme_bw()+
    theme(axis.text.x  = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x=element_text(size = 30),
          axis.title.y = element_text(size = 30))+
    annotate("text",x = .75, y = .25,
             label = paste("AUC of 1 year = ",round(auc_365,2)),color = "blue",size=10)+
    annotate("text",x = .75, y = .15,
             label = paste("AUC of 3 years = ",round(auc_1095,2)),color = "red",size=10)+
    annotate("text",x = .75, y = .05,
             label = paste("AUC of 5 year = ",round(auc_1825,2)),color = "orange",size=10)+
    scale_x_continuous(name = "False Positive Rate")+
    scale_y_continuous(name = "True Positive Rate")
}

ggsave("OP/Pipeline OP5/Train OP/Images/TCGA_train_log2counts - ROC.tiff", plot = last_plot(), width = 900/72, height = 900/72, dpi = 300)


##### risk level and KMplot ######
cutoff <- median(as.numeric(lasso.prob[ ,1]))
cutoff
# [1] 0.5370904 - old
# [1] -0.5888477 - new TRAIN

coxdata$Risk.Score.Level <- ifelse(coxdata$Risk.Score> cutoff, 'High', 'Low')

table(coxdata$Risk.Score.Level)
# High  Low 
#  233  233 

fit.surv <- Surv(time = coxdata$OS.time, event = coxdata$OS)
km <- survfit(fit.surv ~ 1, data = coxdata)

#names(coxdata)[names(coxdata) == "Risk Score Level"] <- "RiskScoreLevel"

km_2 <- survfit(fit.surv ~ Risk.Score.Level, data = coxdata)
dat.survdiff <- survdiff(Surv(OS.time,OS) ~ Risk.Score.Level, data = coxdata)
p.val <- 1 - pchisq(dat.survdiff$chisq, length(dat.survdiff$n) - 1)

plot <- 
  ggsurvplot(km_2,
             legend = "top",
             legend.title = "Risk Score",
             legend.labs = c("High", "Low"),
             font.tickslab = c(18),
             linetype =  c("solid", "dashed"),
             xlab="Days",
             ylab="Overall Survival",
             pval=TRUE,
             pval.method=TRUE, 
             pval.size = 8,
             pval.coord = c(0, 0.03),
             pval.method.coord = c(0,0.1),
             risk.table = "abs_pct",
             risk.table.col= 'strata',
             #risk.table.fontsize = 6,
             fontsize = 5.5,
             palette= c("High" = "#ED0000FF", "Low" = "#00468BFF"),
             ggtheme = theme_classic2(base_size=22),
  )

# Extract the plot and table
plot_surv <- plot$plot
plot_table <- plot$table + theme(legend.position = "none")

# Arrange the plot and table together
plot_combined <- plot_grid(
  plot_surv, plot_table,
  ncol = 1,
  align = 'v',
  rel_heights = c(2, 1)
)

print(plot_combined)

# Save the combined plot ####
ggsave("OP/Pipeline OP5/Train OP/Images/TCGA_train_log2counts - KM Risk Score.tiff", plot = plot_combined, width = 900/72, height = 900/72, dpi = 300)


####### Risk score analysis, survival analysis and prognostic performance of a riskscore model based on differential expression genes #######
fp_dat <- data.frame(s = 1:length(lasso.prob[ ,1]), v = as.numeric(sort(lasso.prob[ ,1])))
fp_dat$`Risk Group` <- ifelse(fp_dat$v > cutoff, 'High', 'Low')

sur_dat <- data.frame(s = 1:length(lasso.prob[ ,1]),
                      t = coxdata[names(sort(lasso.prob[ ,1])), 'OS.time'],
                      OS = coxdata[names(sort(lasso.prob[ ,1])), 'OS'])

sur_dat$`Vital Status` <- ifelse(sur_dat$OS == 0, 'Alive', 'Death')
sur_dat$`Vital Status`  <- factor(sur_dat$OS, levels = c("Death", "Alive"))
exp_dat <- coxdata[names(sort(lasso.prob[ ,1])), lasso_gene]

save(exp_dat, file="OP/Pipeline OP5/Train OP/TCGA_train_log2counts - exp_dat.RData")

{
  #### distribution of risk score ####
  plot.point <- ggplot(fp_dat, aes(x = s, y = v)) +
    geom_point(aes(color = `Risk Group`), size = 2) +
    scale_colour_manual(values = c("#ED0000FF","#00468BFF")) +
    theme_bw() + 
    theme(axis.text = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 24),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 18)) +
    labs(x = "Patient ID(increasing risk score)", y = "Risk Score") +
    geom_hline(yintercept = cutoff, #median(fp_dat$v),
               colour = "black", linetype = "dotted", size=0.8) +
    geom_vline(xintercept = sum(fp_dat$`Risk Group` == "Low"), colour = "black", linetype = "dotted",size = 0.8)
  
  print(plot.point)
  
  ##### distribution of survival time #####
  plot.sur <- ggplot(sur_dat, aes(x = s, y = t)) +
    geom_point(aes(col = `Vital Status`), size = 2) +
    theme_bw() +
    theme(axis.text = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 24),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 18)) +
    scale_colour_manual(values = c("#ED0000FF", "#00468BFF")) +
    labs(x = "Patient ID(increasing risk score)", y = "Survival time(year)") +
    geom_vline(xintercept = sum(fp_dat$`Risk Group` == "Low"), colour = "black",
               linetype = "dotted", size = 0.8)
  
  print(plot.sur)
  
  ##### heatmap of genes expression levels ####
  
  mycolors2 <- colorRampPalette(c("#0000ff", "white", "#ff0000"), bias = 1.2)(100)
  tmp = t(scale(exp_dat))
  #tmp<-dat_test[,lasso_gene]
  #tmp<-t(scale(tmp))
  tmp[tmp > 1] = 1
  tmp[tmp < -1] = -1
  plot.h = pheatmap(tmp,
                    show_rownames = T,
                    show_colnames = F,
                    legend = F,
                    col= mycolors2,
                    cluster_cols = F,
                    fontsize = 14,
                    #cluster_rows = F
  )
  
  library(cowplot)
  combined <- plot_grid(plot.point, plot.sur, plot.h$gtable,
                        #          labels = c("A","B","C"),
                        align = 'v',
                        axis = "bt",
                        ncol = 1)  
}

print(combined)

# Save the combined plot ####
ggsave("OP/Pipeline OP5/Train OP/Images/TCGA_train_log2counts - distribution and heatmap.tiff", plot = combined, width = 1100/72, height = 751/72, dpi = 300)


##### multiple COX regression ##### 
# prepare coxdata for merging - assign rownames as a column so you can use the merge() function
#coxdata <- tibble::rownames_to_column(coxdata, "Patient_ID")

# prepare pheno data so values in Patient ID column match the coxdata syntax
# replace - with . to match the colnames in the counts dataframe.
#tmp <- str_replace_all(glioma.pheno$Patient_ID, "-", ".")
#glioma.pheno$Patient_ID <- tmp

clean.WHO.CNS5_v2

# merge pheno data with coxdata
pheno.coxdata <- merge(coxdata, clean.WHO.CNS5_v2, by = 0, x.all=TRUE, y.all=FALSE)
pheno.coxdata <- pheno.coxdata[,-(56:57)]

new_dat_cox <- pheno.coxdata

new_dat_cox[new_dat_cox == "unknown"] <- NA

new_dat_cox <- droplevels.data.frame(new_dat_cox)

new_dat_cox <- new_dat_cox %>% 
  dplyr::rename(`Age (Years)` = Age_years,
                MGMTp.Status = MGMT.promoter.status.x,
                TERTp.Status = TERT.promoter.status,
                WHO.CNS5.Diagnosis = WHO_CNS5_diagnosis,
                Histology = WHO_CNS5_histology,
                Grade = WHO_CNS5_grade,
                Codel.Status = WHO_CNS_1p.19q._codeletion,
                IDHm.Status = WHO_CNS5_IDHstatus
  )

# fix up table to have ages in years instead of day
new_dat_cox$Age <- ifelse(new_dat_cox$`Age (Years)` > 39, 'Over 39', 'Under 39')
new_dat_cox$age_level2 <- ifelse(new_dat_cox$Age == 'Over 39', '1', '0')

table(new_dat_cox$Age, useNA="ifany")
#  Over 39 Under 39     <NA> 
#.     267      159       39 

# get mean and SD of Age
new_dat_cox$`Age (Years)` <- as.numeric(new_dat_cox$`Age (Years)`)
mean(new_dat_cox$`Age (Years)`, na.rm = TRUE)
# [1] 47.65728 - Train set
sd(new_dat_cox$`Age (Years)`, na.rm = TRUE)
# [1] 15.52152 - Train set

new_dat_cox$Gender <- as.character(new_dat_cox$Gender) #1 female,2 male

table(new_dat_cox$Gender, useNA="ifany")
# female   male   <NA> 
#.  174    252     39 


new_dat_cox$Grade <- as.character(new_dat_cox$WHO_CNS5_grade)

table(new_dat_cox$Grade, useNA="ifany")
# Grade 2 Grade 3 Grade 4    <NA> 
#.    132     109     166      58 

table(new_dat_cox$IDHm.Status, useNA="ifany")
#  IDH-mutant IDH-wildtype         <NA> 
#.        256          152           57 

table(new_dat_cox$MGMTp.Status, useNA="ifany")
#   Methylated Unmethylated         <NA> 
#          326          112           27 

table(new_dat_cox$Codel.Status, useNA="ifany")
# Codel Non-codel      <NA> 
# 117       344         4 

table(new_dat_cox$Histology, useNA="ifany")
#       astrocytoma      glioblastoma  oligoastrocytoma oligodendroglioma -OLD
#            115               115                76               121 - Train set
# Astrocytoma      Glioblastoma Oligodendroglioma              <NA> - NEW
#         150               152               106                57 

#new_dat_cox$KPS_level<-ifelse(new_dat_cox$KPS > 80,'High','Low')
#new_dat_cox$KPS_level2<-ifelse(new_dat_cox$KPS_level == 'High','1','0')
#sum(is.na(new_dat_cox$KPS_level2))
# [1] 203 - Train set
#table(new_dat_cox$KPS_level2)
# 0   1 
# 120 143  - Train set

new_dat_cox$Risk.Score.Level2 <- ifelse(new_dat_cox$Risk.Score.Level == 'High','1','0')
table(new_dat_cox$Risk.Score.Level2, useNA="ifany")
#   0   1 
# 233 233 

# clean new_dat_cox df. remove unwanted columns.
tmp <- new_dat_cox[, c(1:57, 71, 83:86, 89)]

new_dat_cox <- tmp
new_dat_cox <- column_to_rownames(new_dat_cox, var="Row.names")

save(new_dat_cox, file = "OP/Pipeline OP5/Train OP/TCGA_train_log2counts - new_dat_cox.RData")

write.csv(new_dat_cox,'OP/Pipeline OP5/Train OP/CSV/TCGA_train_log2counts - new_dat_cox.csv',row.names = FALSE)
###################### Get summary of clinical variables in new_dat_cox ###################### 
lapply(new_dat_cox[,54:62], table, useNA = "always")

new_dat_cox %>%
  group_by(Histology, Grade) %>%
  summarise(count = n())

###################### Hazards ratyion analysis for each gene in lasso_gene ###################### 
lasso_gene

lasso_res = ezcox(new_dat_cox,
             covariates = lasso_gene,
             time = 'OS.time.x',
             status = 'OS.x',
             global_method = c("likelihood", "wald", "logrank")
)

lasso_res

write.csv(lasso_res,'OP/Pipeline OP5/Train OP/CSV/TCGA_train_log2counts - lasso-gene_uniCOX.csv',row.names = FALSE)
###################### keep this unchanged from train set, line 806 to 825 ###################### 

k<-c('Gender', 'Grade', 'Age', 'Risk.Score.Level', 'Codel.Status', 'IDHm.Status',  'MGMTp.Status')

res3 = ezcox(new_dat_cox,
            covariates = k,
            time = 'OS.time.x',
            status = 'OS.x',
            global_method = c("likelihood", "wald", "logrank")
)

res3

write.csv(res3,'OP/Pipeline OP5/Train OP/CSV/TCGA_train_log2counts - clinical_uniCOX.csv',row.names = FALSE)

# multiple Cox regession analysis
# multi <- c('Gender', 'Grade', 'KPS_level', 'age_level', 'Risk_Score_level', 'Codel', 'IDH.status',  'MGMT_methy')
#formula_multicox <- as.formula(paste0('Surv(OS.time,OS)~', paste(multi,sep = '',collapse = '+')))
formula_multicox <- as.formula(paste0('Surv(OS.time.x,OS.x)~', paste(k, sep = '',collapse = '+')))
multi_cox <- coxph(formula_multicox,data = new_dat_cox)
summary(multi_cox)
# Call:
#   coxph(formula = formula_multicox, data = new_dat_cox)

# n= 377, number of events= 141 
# (88 observations deleted due to missingness)

#                             coef exp(coef) se(coef)      z Pr(>|z|)    
# Gendermale                0.1498    1.1616   0.1833  0.817 0.413832    
# GradeGrade 3              0.6582    1.9312   0.3413  1.928 0.053836 .  
# GradeGrade 4              1.9601    7.1003   0.5084  3.856 0.000115 ***
# AgeUnder 39              -0.7810    0.4579   0.2828 -2.762 0.005753 ** 
# Risk.Score.LevelLow      -0.8405    0.4315   0.3267 -2.573 0.010086 *  
# Codel.StatusNon-codel     0.2893    1.3355   0.3421  0.846 0.397649    
# IDHm.StatusIDH-wildtype   0.6090    1.8386   0.4644  1.311 0.189690    
# MGMTp.StatusUnmethylated  0.4007    1.4929   0.2093  1.915 0.055502 .  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#proportional hazards over time - Global Schoenfeld Test ####
test.ph <- cox.zph(multi_cox)
ggcoxzph(test.ph)

# open a TIFF device to save the plot 
tiff(filename = "OP/Pipeline OP5/Train OP/CSV/TCGA_train_log2counts - Global Schoenfeld test.tiff",
     width = 10, height = 10, units = "in", res = 300)

# Now, create your plot
ggcoxzph(test.ph)

# Close the device
dev.off()

save(new_dat_cox, file = "OP/Pipeline OP5/Train OP/CSV/TCGA_train_log2counts - new_dat_cox.RData")

####----VII. Develepmont and evaluation of the Nomogram----####

nom_dat_cox <- new_dat_cox

print(sapply(nom_dat_cox, function(x) sum(is.na(x)))) # Check for NAs

#nom_dat_cox <- nom_dat_cox[ ,-1] # don't need patient ID

#nom_dat_cox <- column_to_rownames(nom_dat_cox, var = "Patient_ID") 

names(nom_dat_cox)

# do this for test set
dd <- datadist(nom_dat_cox)
options(datadist="dd")
#nom_cox <- cph(Surv(OS.time,OS) ~ Grade + age_level2 + Risk_Score_level2, data = new_dat_cox, x = T, y = T, surv = T)
nom_cox <- cph(Surv(OS.time.x,OS.x) ~ Grade + Age + Risk.Score.Level, data = nom_dat_cox, x = T, y = T, surv = T)

###################### keep above unchanged from train set ###################### 

# plot nomogram ####
tiff(filename = "OP/Pipeline OP5/Train OP/Images/TCGA_train_log2counts - Nomogram3.tiff",
     width = 13, height = 10, units = "in", res = 300)

{
  survival <- Survival(nom_cox)
  survival1 <- function(x)survival(365,x)
  survival3 <- function(x)survival(1095,x)
  survival5 <- function(x)survival(1825,x)
  
  
  nom <- nomogram(nom_cox, lp = F ,
                  fun = list(survival1,survival3,survival5),
                  fun.at = c(0.1,seq(0.1,0.9,by = 0.2), 0.9),
                  funlabel = c("1-year survival",'3-year survival','5-year survival'))
#  plot(nom, xfrac=.2,
#       cex.axis = 1.05,
#       cex = 1.5,
#       force.label = F,
#       tcl = 0.8,
#       lmgp = 0.1,
#       vnames="labels",
#       col.grid=gray(c(0.85,0.95)))
  
  plot(nom, xfrac=0.5,  # Increased xfrac to give more space for labels
       cex.axis = 1.8,
       cex.var = 1.8,   # Increase size of variable names
       cex = 1.8,       # Increase overall text size
       lmgp = 0.5,      # Adjust label margin
       tcl = -0.5,      # Adjust tick mark length
       label.every = 1, # Label every tick mark
       force.label = TRUE,
       vnames="labels",
       col.grid=gray(c(0.85,0.95)))
  
}

# Close the device
dev.off()

# 1-year calibration curve ####
coxm_1 <- cph(Surv(OS.time.x,OS.x) ~ Grade + Age + Risk.Score.Level,
              data = nom_dat_cox,
              surv = T,x = T,y = T,
              time.inc = 365)

cal_1 <- calibrate(coxm_1,
                   u = 365,
                   cmethod = 'KM',
                   m = 60, # previously 60
                   B = 1000) # previously 1000

#par(mar = c(7,4,4,3), cex = 1.0)

# 3-year calibration curve ####
coxm_3 <- cph(Surv(OS.time.x,OS.x) ~ Grade + Age + Risk.Score.Level,
              data = nom_dat_cox,
              surv = T, x = T, y = T,
              time.inc = 365*3)

cal_3 <- calibrate(coxm_3, u = 365 * 3, cmethod = 'KM', m = 60, B = 1000)

# 5-year calibration curve ####
coxm_5 <- cph(Surv(OS.time.x,OS.x) ~ Grade + Age + Risk.Score.Level,
              data = nom_dat_cox,
              surv = T, x = T, y = T,
              time.inc = 365 * 5)

cal_5 <- calibrate(coxm_5, u = 365 * 5, cmethod = 'KM', m = 60, B = 1000)

# plot calibration curves ####

# open a TIFF device to save the plot
tiff(filename = "OP/Pipeline OP5/Train OP/Images/TCGA_train_log2counts - Nomogram Calibration.tiff",
     width = 10, height = 10, units = "in", res = 300)

# Set global graphical parameters
par(cex.axis = 1.4,    # Size of axis text
    cex.lab = 1.5,     # Size of axis labels
    cex.main = 1.6,    # Size of title (if you add one)
    cex.sub = 1.2,     # Size of subtitle (if you add one)
    mar = c(6, 5, 4, 2) + 0.1)  # Increase left margin (second value)
        
# Now, create your plot
plot(cal_1, lwd = 2, lty = 1,
     errbar.col = c(rgb(0, 118, 192, maxColorValue = 255)),
     xlab = 'Nomogram-Predicted Probability',
     ylab = 'Actual Proportion',
     col = '#42b540', # c(rgb(0,70,139, maxColorValue = 255)),
     xlim = c(0, 1), ylim = c(0, 1),
     add = FALSE)

plot(cal_3, lwd = 2, lty = 1,
     errbar.col = c(rgb(0, 118, 192, maxColorValue = 255)),
     xlab = 'Nomogram-Predicted Probability',
     ylab = 'Actual Proportion',
     col = '#ed0000', # c(rgb(237,0,0, maxColorValue = 255)),
     xlim = c(0,1),ylim = c(0,1),
     add = TRUE)

plot(cal_5, lwd = 2, lty = 1,
     errbar.col = c(rgb(0, 118, 192, maxColorValue = 255)),
     xlab = 'Nomogram-Predicted Probability',
     ylab = 'Actual Proportion',
     col = '#00468b', #c(rgb(66,181,64, maxColorValue = 255)),
     xlim = c(0, 1),ylim = c(0, 1),
     add = TRUE) 

# Add legend
legend("bottomright", 
       legend = c("1-year survival", "3-year survival", "5-year survival"), 
       col = c("#42b540", "#ed0000", "#00468b"), 
       lwd = 2,
       cex = 1.5)  # Adjust this value to change legend text size

# Close the device
dev.off()


# C-index
f <- coxph(Surv(OS.time.x,OS.x) ~ Grade + Age + Risk.Score.Level, data = nom_dat_cox)
sum.surv <- summary(f)
c_index <- sum.surv$concordance
c_index
#          C      se(C) 
# 0.82502128 0.01291766 

# timeROC for nomogram model
pred_nom <- predict(f, newdata = nom_dat_cox)
nom_dat_cox$pred_nom <- as.numeric(pred_nom)

with(nom_dat_cox,
     ROC <<- timeROC(T = OS.time.x,
                     delta = OS.x,
                     marker = pred_nom,
                     cause = 1,
                     weighting = "marginal",
                     times = c(365, 1095, 1825),
                     ROC = TRUE,
                     iid = TRUE)
)

ROC$AUC
#     t=365    t=1095    t=1825 
# 0.8515814 0.9498082 0.8646947  

set.seed(100)
confint(ROC)
#$CI_AUC
#2.5% 97.5%
#t=365  81.75 88.57
#t=1095 92.11 97.85
#t=1825 80.04 92.90

#$CB_AUC
#2.5% 97.5%
#t=365  81.00 89.32
#t=1095 91.48 98.48
#t=1825 78.63 94.31

#$C.alpha
#95% 
#2.390582  

{
  auc_365 = ROC$AUC[[1]]
  auc_1095 = ROC$AUC[[2]]
  auc_1825 = ROC$AUC[[3]]
  dat<-data.frame(tpr365=ROC$TP[,1],
                  fpr365=ROC$FP[,1],
                  tpr1095=ROC$TP[,2],
                  fpr1095=ROC$FP[,2],
                  tpr1825=ROC$TP[,3],
                  fpr1825=ROC$FP[,3])
  library(ggplot2)
  ggplot() +
    geom_line(data = dat,aes(x = fpr365, y = tpr365),color = "blue") +
    geom_line(data = dat,aes(x = fpr1095, y = tpr1095),color = "red")+
    geom_line(data = dat,aes(x = fpr1825, y = tpr1825),color = "orange")+
    geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
    theme_bw()+
    theme(axis.text.x  = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x=element_text(size = 30),
          axis.title.y = element_text(size = 30))+
    annotate("text",x = .75, y = .25,
             label = paste("AUC of 1 year = ",round(auc_365,2)),color = "blue",size=10)+
    annotate("text",x = .75, y = .15,
             label = paste("AUC of 3 years = ",round(auc_1095,2)),color = "red",size=10)+
    annotate("text",x = .75, y = .05,
             label = paste("AUC of 5 year = ",round(auc_1825,2)),color = "orange",size=10)+
    scale_x_continuous(name = "False Positive Rate")+
    scale_y_continuous(name = "True Positive Rate")
}

ggsave("OP/Pipeline OP5/Train OP/Images/TCGA_train_log2counts - Nomogram ROC.tiff", plot = last_plot(), width = 900/72, height = 751/72, dpi = 300)

new_tcga_nom <- nom_dat_cox

save(new_tcga_nom, file = 'OP/Pipeline OP5/Train OP/TCGA_train_log2counts - new_tcga_nom.Rdata')


####----VIII. Stratified Analysis----####

# expression levels of identified genes

View(new_tcga_nom)

new_tcga_nom$`Risk Group` <- new_tcga_nom$Risk.Score.Level
new_tcga_nom$Risk.Score.Level <- NULL

# Create stack bar graphs to recreate what you have in the paper.

# Histology ####

# First, filter out NA values and calculate percentages
tmp <- new_tcga_nom %>%
  filter(!is.na(Histology)) %>%
  group_by(`Risk Group`, Histology) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

# Then create the plot
library("ggsci")

g <- ggplot(tmp, aes(x = `Risk Group`, y = count, fill = Histology)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = sprintf("%.2f%%", percentage)), 
            position = position_stack(vjust = 0.5),
            size = 6) +  # Adjust size as needed
  scale_fill_lancet() +
  labs(x = "Risk Group", y = "No. of Patient Samples") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "right"
  )

# Display the plot
print(g)

ggsave("OP/Pipeline OP5/Train OP/Images/TCGA_train_log2counts - Nomogram ROC.tiff", plot = last_plot(), width = 900/72, height = 751/72, dpi = 300)




# Code below is from old script #####

### Histology

dat <- new_tcga_nom
tmp1 <- subset(dat, dat$Histology == 'Astrocytoma')
tmp1 <- data.frame(tmp1$Histology, tmp1$Risk.Score)
colnames(tmp1) <- c('Histology', 'Risk Score')
tmp2 <- subset(dat,dat$Histology == 'Oligodendroglioma')
tmp2 <- data.frame(tmp2$Histology, tmp2$Risk.Score)
colnames(tmp2) <- c('Histology', 'Risk Score')
tmp3 <- subset(dat,dat$Histology == 'Glioblastoma')
tmp3 <- data.frame(tmp3$Histology, tmp3$Risk.Score)
colnames(tmp3) <- c('Histology', 'Risk Score')

dat_t <- rbind(tmp1, tmp2, tmp3)

# open a TIFF device to save the plot
tiff(filename = "OP/Pipeline OP5/Train OP/Images/TCGA_train_log2counts - Strat_histology.tiff",
     width = 12, height = 10, units = "in", res = 300)

p <- ggboxplot(dat_t, x = "Histology", y = 'Risk Score',
               color = "Histology", palette = "lancet",
               add = "jitter") + 
  stat_compare_means(size = 6) +
  theme(
    axis.title = element_text(size = 18, face = "bold"),  # Increase axis label size
    axis.text = element_text(size = 18),                  # Increase axis text size
    legend.title = element_text(size = 18, face = "bold"),# Increase legend title size
    legend.text = element_text(size = 18),                # Increase legend text size
    plot.title = element_text(size = 18, face = "bold"),  # Increase plot title size if you have one
    text = element_text(size = 18),                        # Increase all other text sizes
    plot.margin = margin(t = 20, r = 20, b = 21, l = 21, unit = "pt")  # Increase margins
    ) +
  labs(x = "Histology", y = "Risk Score")                 # Explicitly set axis labels if needed

# Print the plot
print(p)

dev.off()

### IDH Mutation

tmp1 <- subset(dat,dat$IDH.status == 'Mutant')
tmp1 <- data.frame(tmp1$IDH.status, tmp1$Risk_score)
colnames(tmp1) <- c('IDH.status', 'risk_score')
tmp2 <- subset(dat,dat$IDH.status == 'WT')
tmp2 <- data.frame(tmp2$IDH.status, tmp2$Risk_score)
colnames(tmp2) <- c('IDH.status', 'risk_score')
dat_IDH <- rbind(tmp1, tmp2)

p2 = ggboxplot(dat_IDH, x = "IDH.status", y = "risk_score",
               color = "IDH.status", palette = "lancet",
               add = "jitter")+ stat_compare_means()
p2

### Age

tmp1 <- subset(dat,dat$age_level == 'Under 39')
tmp1 <- data.frame(tmp1$age_level, tmp1$Risk_score)
colnames(tmp1) <- c('age_level', 'risk_score')
tmp2 <- subset(dat,dat$age_level == 'Over 39')
tmp2 <- data.frame(tmp2$age_level, tmp2$Risk_score)
colnames(tmp2) <- c('age_level', 'risk_score')
dat_age <- rbind(tmp1, tmp2)

p3 = ggboxplot(dat_age, x = "age_level", y = "risk_score",
               color = "age_level", palette = "lancet",
               add = "jitter")+ stat_compare_means()
p3

### Tumour Grade
# Should Change glioblastoma grade to 4 instead of NA
n = nrow(dat)

dat$Grade <- as.character(dat$Grade) # change this column to characters so following code can work

#### assign grade 4 for glioblastoma. probably better to be moved up in the script.
dat <- within(dat, Grade[Grade == 'NA' & Histology == 'glioblastoma'] <- '4')

tmp1 <- subset(dat, dat$Grade == '2')
tmp1 <- data.frame(tmp1$Grade, tmp1$Risk_score)
colnames(tmp1) <- c('Grade', 'risk_score')
tmp2 <- subset(dat,dat$Grade == '3')
tmp2 <- data.frame(tmp2$Grade, tmp2$Risk_score)
colnames(tmp2) <- c('Grade', 'risk_score')
tmp3 <- subset(dat,dat$Grade == '4')
tmp3 <- data.frame(tmp3$Grade, tmp3$Risk_score)
colnames(tmp3) <- c('Grade', 'risk_score')
dat_grade <- rbind(tmp1, tmp2, tmp3)

p4 = ggboxplot(dat_grade, x = "Grade", y = "risk_score",
               color = "Grade", palette = "lancet",
               add = "jitter")+ stat_compare_means()
p4

### MGMT Methylation

tmp1 <- subset(dat, dat$MGMT_methy == 'Methylated')
tmp1 <- data.frame(tmp1$MGMT_methy, tmp1$Risk_score)
colnames(tmp1) <- c('MGMT_methy', 'risk_score')
tmp2 <- subset(dat,dat$MGMT_methy == 'Unmethylated')
tmp2 <- data.frame(tmp2$MGMT_methy, tmp2$Risk_score)
colnames(tmp2) <- c('MGMT_methy', 'risk_score')
dat_MGMT <- rbind(tmp1, tmp2)

p5 = ggboxplot(dat_MGMT, x = "MGMT_methy", y = "risk_score",
               color = "MGMT_methy", palette = "lancet",
               add = "jitter")+ stat_compare_means()
p5

### 1p/19q codeletion

tmp1 <- subset(dat, dat$Codel == 'codel')
tmp1 <- data.frame(tmp1$Codel, tmp1$Risk_score)
colnames(tmp1) <- c('Codel', 'risk_score')
tmp2 <- subset(dat,dat$Codel == 'non-codel')
tmp2 <- data.frame(tmp2$Codel, tmp2$Risk_score)
colnames(tmp2) <- c('Codel', 'risk_score')
dat_codel <- rbind(tmp1, tmp2)

p6 = ggboxplot(dat_codel, x = "Codel", y = "risk_score",
               color = "Codel", palette = "lancet",
               add = "jitter")+ stat_compare_means()
p6

### forest plot for unicox

fp_df <- unicox[ ,c(1,8,9,10,12)]
fp_df$Index <- 1:nrow(fp_df)

ggplot(data=fp_df, aes(y=Index, x=HR,
                             xmin=lower_95, 
                             xmax=upper_95)) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  scale_y_continuous(name="", breaks=fp_df$Index, labels=unicox$Variable) +
  labs(title='Hazards Ratio by Gene', x='Hazards Ratio (HR)', y = 'Copper Genes') +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme_minimal()

### to understand if MT1X expression on survival
md <- median(as.numeric(new_tcga_nom$MT1X))
new_tcga_nom$MT1X_level <- ifelse(new_tcga_nom$MT1X > md, 'MT1X High', 'MT1X Low')
table(new_tcga_nom$MT1X_level)

### to understand if there's a correlation betweem IDH.status and MT1X expression
new_tcga_nom$IDH.status_MT1X <- paste(new_tcga_nom$IDH.status, new_tcga_nom$MT1X_level)
table(new_tcga_nom$IDH.status_MT1X)

## to understand the correlation between IDH.status and Risk Score
new_tcga_nom$IDH.status_RS <- paste(new_tcga_nom$IDH.status, new_tcga_nom$Risk_Score_level)
table(new_tcga_nom$IDH.status_RS)

tmp <- new_tcga_nom[!is.na(new_tcga_nom$IDH.status),]
table(tmp$IDH.status_RS)


fit.surv <- Surv(time = tmp$OS.time, event = tmp$OS)
km <- survfit(fit.surv ~ 1, data = tmp)
km_2 <- survfit(fit.surv ~ IDH.status_RS, data = tmp)
dat.survdiff <- survdiff(Surv(OS.time,OS) ~ IDH.status_RS, data = tmp)
p.val <- 1 - pchisq(dat.survdiff$chisq, length(dat.survdiff$n) - 1)

ggsurvplot(km_2,
           legend = "top",
           legend.title = "MT1X expression",
           legend.labs = c('Mutant High', 'Mutant Low', 'WT High', 'WT Low'),
           font.tickslab = c(18),
           #linetype =  c("solid", "dashed"),
           xlab="Days",
           ylab="Overall Survival",
           pval=TRUE,
           pval.method=TRUE, 
           pval.size = 8,
           pval.coord = c(0, 0.03),
           pval.method.coord = c(0,0.1),
           risk.table = "abs_pct",
           risk.table.col= 'strata',
           #risk.table.fontsize = 6,
           fontsize = 5.5,
           palette= "lancet",
           ggtheme = theme_classic2(base_size=22),
)

### create a scatter plot to look at risk score distribution and immune score

x <- new_tcga_nom$ESTIMATE.immune.score
y <- new_tcga_nom$Risk_score

# Plot with main and axis titles
# Change point shape (pch = 19) and remove frame.
plot(x, y, main = "Risk Score vs ESTIMATE Immune Score",
     xlab = "ESTIMATE Immune Score", ylab = "Risk Score",
     pch = 19, frame = FALSE)

# Add regression line
plot(x, y, main = "Risk Score vs ESTIMATE Immune Score",
     xlab = "ESTIMATE Immune Score", ylab = "Risk Score",
     pch = 19, frame = FALSE)
abline(lm(y ~ x, data = new_tcga_nom), col = "blue")

# Scatter plot by groups ("cyl")
library(car)
scatterplot(Risk_score ~ ESTIMATE.immune.score | IDH.status, data = new_tcga_nom, 
            smoother = FALSE, grid = FALSE, frame = FALSE)

# Prepare the data set
x <- new_tcga_nom$ESTIMATE.immune.score
y <- new_tcga_nom$ESTIMATE.combined.score
z <- new_tcga_nom$Risk_score

grps <- as.factor(new_tcga_nom$IDH.status)
# Plot
library(scatterplot3d)
colors <- c("#999999", "#E69F00", "#56B4E9")
scatterplot3d(x, y, z, pch = 16, color = colors[grps],
              grid = TRUE, box = FALSE, xlab = "ESTIMATE.immune.score", 
              ylab = "ESTIMATE.combined.score", zlab = "Risk_score")

x <- new_tcga_nom$ESTIMATE.stromal.score
y <- new_tcga_nom$ESTIMATE.combined.score
z <- new_tcga_nom$Risk_score

grps <- as.factor(new_tcga_nom$IDH.status)

# Plot
library(scatterplot3d)
colors <- c("#999999", "#E69F00", "#56B4E9")
scatterplot3d(x, y, z, pch = 16, color = colors[grps],
              grid = TRUE, box = FALSE, xlab = "ESTIMATE.stromal.score", 
              ylab = "ESTIMATE.combined.score", zlab = "Risk_score")

#create a coxdata file with only the LASSO genes in it. Rownames = genes, colnames = Patient.ID
cox_tmp <- glioma.log2.coxdata.unique[, c(lasso_gene)]
cox_tmp <- data.frame(t(cox_tmp))

write.csv(cox_tmp,file = '/Users/tonijue/Documents/Copper in glioma/Copper in glioma/TCGA LGG 20220523/Scripts/Version 3_2023/20231113_lasso gene expression.csv',row.names = TRUE)

# create df of lass gene expression that has only high risk scrores
HighRisk_lassoGenes_trainSet <- new_dat_cox %>% filter(Risk_Score_level == "High")

HighRisk_PtID <- c(HighRisk_lassoGenes_trainSet$Patient_ID)

HighRisk_train <- glioma.raw.unique[,names(glioma.raw.unique) %in% HighRisk_PtID]  #worked!

#### check if any of the lasso_genes are transcription factors####

tf.file <- "/Users/tonijue/Documents/Copper in glioma/Copper in glioma/TCGA LGG 20220523/Scripts/Version 3_2023/Files/TF_gene_names.txt"

tf <- read.delim(tf.file, header = FALSE, row.names = FALSE, sep = "\t", dec = ".")

head(tf)

tf_check <- tf[, 'V1'] %in% lasso_gene

summary(tf_check)

