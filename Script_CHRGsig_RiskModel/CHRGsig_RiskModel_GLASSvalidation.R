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

load("~/Documents/Copper in glioma/Copper in glioma/TCGA LGG 20220523/Scripts/Version 3_2023/Version 3_2023_rebuttal/OP/Pipeline OP4/GLASS_log2counts_geneIDs.RData") # data was log2 transformed. we have yet to see whether this is enough since there is quite a difference between GLASS and TCGA data based on the PCA plot.

####---- VI. Univariate Cox Analysis ----####
dim(glass_aggregated_log2)
# [1] 25826   374

view(head(glass_aggregated_log2))

#glass_adjusted <- column_to_rownames(glass_aggregated_log2, var="gene")
glass_adjusted <-  glass_aggregated_log2
tmp <- transpose(glass_adjusted)

colnames(tmp) <- rownames(glass_adjusted)
rownames(tmp) <- colnames(glass_adjusted)

glass_adjusted.T <- tmp


# Univariate Cox Analysis must use log2(count+1) data

### A. First select differentially expressed copper genes from DEGlist
copper.genelist <- c('ABCB6', 'ANKRD9', 'SLC31A1', 'SLC31A2', 'PRND', 'CCDC22', 'APP', 'ARF1', 'MT2A', 'ATOX1', 'ATP7A', 'ATP7B', 'PRNP', 'SCO1', 'COX19', 'SCO2', 'CYP1A1', 'DAXX', 'BACE1', 'AOC1', 'MT1DP', 'HSF1', 'AQP1', 'AQP2', 'MT1A', 'MT1B', 'MT1E', 'MT1F', 'MT1G', 'MT1H', 'MT1M', 'MT1X', 'MT3', 'NFE2L2', 'MT1HL1', 'SNCA', 'MAP1LC3A', 'MT4', 'BECN1', 'COMMD1', 'XIAP', 'CUTC', 'STEAP2', 'STEAP3', 'STEAP4', 'SLC11A2', 'COX17', 'CP', 'FKBP4', 'HEPHL1', 'MMGT1', 'HEPH', 'PARK7', 'AANAT', 'IL1A', 'LCAT', 'LOXL2', 'MT-CO1', 'PAM', 'ATP5F1D', 'SOD1', 'SOD3', 'SORD', 'TFRC', 'CDK1', 'MOXD2P', 'MTCO2P12', 'COX11', 'LACC1', 'DBH', 'DCT', 'ALB', 'F5', 'F8', 'OR5AR1', 'ADNP', 'ATP13A2', 'MOXD1', 'GPC1', 'ANG', 'SUMF1', 'AOC2', 'SNAI3', 'APOA4', 'COA6', 'LOX', 'LOXL1', 'MT-CO2', 'ACR', 'P2RX4', 'CUTA', 'HAMP', 'S100A5', 'S100A12', 'S100A13', 'SNCB', 'SNCG', 'TP53', 'TYR', 'LOXL4', 'LOXL3', 'AOC3', 'RNF7', 'CCS', 'AP1S1', 'AP1B1', 'TMPRSS6', 'SPATA5', 'COG2', 'ATP6V0A2', 'ATP6AP1', 'ADAM10', 'AKT1', 'MTF2', 'FOXO1', 'FOXO3', 'STEAP1', 'GSK3B', 'APC', 'JUN', 'MAPT', 'MDM2', 'MT1JP', 'MT1L', 'MTF1', 'PIK3CA', 'XAF1', 'PTEN', 'CCND1', 'SP1', 'ADAM17', 'CASP3', 'ADAM9')

# getting the subset DataFrame after checking values if belonging to vector
copper.DEG <- DEGlist[DEGlist$gene %in% copper.genelist, ]

copper.DEGlist <- c(copper.DEG$gene)

### B. create a subset of log2.counts that contains only the copper.DEGlist
glass.log2.copper <- glass_adjusted.T[,copper.DEGlist]
dim(glass.log2.copper)
# [1] 373  50

### C. Need to recreate the coxdata table I had from python which included the expression matrix plus the survival times from the survival file.
dim(clean.GLASS.CNS5.TP)

# Function to extract the desired part of the string
extract_id <- function(x) {
  parts <- strsplit(x, "\\.")
  paste(sapply(parts, function(p) paste(p[1:4], collapse = ".")), collapse = "")
}

# Apply the function to rownames and create a new column
glass.log2.copper$Patient_ID <- NA
glass.log2.copper$Patient_ID <- sapply(rownames(glass.log2.copper), extract_id)
length(unique(glass.log2.copper$Patient_ID))
# [1] 367

# Filter rows where Patient_ID ends with ".TP"
glass.TP.log2.copper <- glass.log2.copper %>%
  filter(str_detect(Patient_ID, "\\.TP$"))
dim(glass.TP.log2.copper)
# [1] 156  51

glass.TP.log2.copper.unique <- aggregate(glass.TP.log2.copper, by=list(glass.TP.log2.copper$Patient_ID), FUN=mean)
glass.TP.log2.copper.unique <- glass.TP.log2.copper.unique[ ,-52]
glass.TP.log2.copper.unique <- glass.TP.log2.copper.unique %>% dplyr::rename(Patient_ID = Group.1)

rownames(glass.TP.log2.copper.unique) <- NULL
glass.TP.log2.copper.unique <- column_to_rownames(glass.TP.log2.copper.unique, var="Patient_ID")

# check seq distrubtion for each gene.
df <- glass.TP.log2.copper.unique

# View distribution of RNAseq data in glass.TP.log2.copper.unique
# open a TIFF device to save the plot
tiff(filename = "OP/Pipeline OP5/GLASS Test OP/Images/GLASS_test_log2transform - copper_DEG.tiff",
     width = 20, height = 10, units = "in", res = 300)

# Set up the plotting area to accommodate 50 plots (5 rows, 10 columns)
par(mfrow = c(5, 10))

# Loop through each column in the dataframe and plot it
for (col_name in colnames(df)) {
  MASS::truehist(df[[col_name]], main = col_name, nbins = 12)
}

dev.off()


# replace - with . in surv file to match the colnames in the counts dataframe.
tmp <- str_replace_all(clean.GLASS.CNS5.TP$sample_barcode, "-", ".")
clean.GLASS.CNS5.TP$sample_barcode <- tmp

rownames(clean.GLASS.CNS5.TP) <- NULL
clean.GLASS.CNS5.TP <- column_to_rownames(clean.GLASS.CNS5.TP, var="sample_barcode")

# take only needed columns
cols_to_keep <- c("case_sex", "case_age_diagnosis_years", "case_vital_status", "case_overall_survival_mo", "mgmt_methylation", "WHO_CNS5_IDHstatus", "WHO_CNS5_grade", "WHO_CNS5_histology", "Codel")

glass.pheno <- clean.GLASS.CNS5.TP[,cols_to_keep]

glass.pheno <- glass.pheno %>% 
  dplyr::rename(Gender = case_sex,
                `Age (Years)` = case_age_diagnosis_years,
                OS.time_mos = case_overall_survival_mo,
                OS = case_vital_status,
                MGMTp.Status = mgmt_methylation,
                IDHm.Status = WHO_CNS5_IDHstatus,
                Grade = WHO_CNS5_grade,
                Histology = WHO_CNS5_histology,
                Codel.Status = Codel
                )

glass.pheno <- glass.pheno %>% drop_na(OS.time_mos)
glass.pheno <- glass.pheno[!glass.pheno$OS.time_mos=="unknown",]

glass.pheno$`Age (Years)`<-as.numeric(glass.pheno$`Age (Years)`) # remove NAs - 16 lines
glass.pheno$OS.time_mos<-as.numeric(glass.pheno$OS.time_mos) # definitely remove NAs - 18 lines

glass.surv <- glass.pheno[,c("OS.time_mos", "OS")]

glass.surv$OS <- ifelse(glass.surv$OS =="dead", 1, 0)

str(glass.surv)

# convert OS.time_mos to days. Mean month length:30.436875
glass.surv$OS.time <- glass.surv$OS.time_mos * 30.436875

glass.surv <- glass.surv[,-1] # remove OS.time_mos column



# merge data into dataframe with rnaseq data of copper DEGs and OS/OS.time
glass.log2.coxdata <- merge(glass.TP.log2.copper.unique, glass.surv, by = 0)
dim(glass.log2.coxdata)
# [1] 153  53

save(glass.log2.coxdata, file = "OP/Pipeline OP5/GLASS Test OP/GLASS_test_log2transform - coxdata.RData")



### D. First Univariate Cox Analysis

# Univariate COX Analysis
#coxdata <- glioma.train.log2.coxdata.unique # TCGA train set

#coxdata <- glioma.test.log2.coxdata.unique # TCGA test set

coxdata <-  glass.log2.coxdata # GLASS Z-score transformed log2counts test set
coxdata <- column_to_rownames(coxdata, var="Row.names")

#coxdata <- glioma.IDHwt.log2.copper # IDHwt set

#coxdata <- glioma.IDHm.log2.copper # IDHm set

# to check distibution of gene expression (Optional) - they should have a normal distribution
#ggplot(coxdata, aes(x = MT1X, color = OS)) + geom_histogram(color = "black", fill = "white", bins = 30) 
#ggplot(coxdata, aes(x = SLC31A1, color = OS)) + geom_histogram(color = "black", fill = "white", bins = 30) 
#ggplot(coxdata, aes(x = MT1M, color = OS)) + geom_histogram(color = "black", fill = "white", bins = 30) 

####---- VI. Constructing and validating risk-score system ----####

######contructing risk-score system, LASSO model#######
load("OP/Pipeline OP5/Train OP/unicox.RData")

unicox_dif <- unicox[which(unicox$global.pval < 0.05), ] # do not change. Use this from TCGA Train set.
uni_gene <- as.character(unicox_dif$Variable) # do not chnage. Use this from TCGA Train set.

# change line coxdata to test set for validation.
irg_expr <- coxdata[ ,c(1:50)] # must include all gene columns from coxdata data frame. Change for TCGA Validation Set.
# fix(coxdata) # opens R editor. not sure what for.

x <- as.matrix(irg_expr[ ,uni_gene]) # uni_gene is taken from the train dataset. 
y <- data.matrix(Surv(coxdata$OS.time,coxdata$OS)) 

###################### keep this unchanged from train set ###################### 
#cv.fit <- cv.glmnet(x, y,
#                    family="cox",
#                    maxit = 1000000,
#                    alpha=1)
#print(cv.fit)
# Call:  cv.glmnet(x = x, y = y, family = "cox", maxit = 1e+06, alpha = 1) 

# Measure: Partial Likelihood Deviance 

#      Lambda Index Measure     SE Nonzero
# min 0.03357    26   10.64 0.2180      19
# 1se 0.10253    14   10.82 0.2325      11

# Cv.fit based on training set

#Active.Index
# [1]  4  5  7 11 15 16 17 18 19 20 22 24 27 30 31 37 40 43 44

#Active.Coefficients
# [1] -0.007409754  0.155383500  0.073934367  0.217785923  0.032618756  0.038602690 -0.167637522  0.070728710
# [9]  0.099598934  0.065588082  0.029122490 -0.062421631  0.006042011  0.087743960 -0.050878675 -0.001973116
# [17]  0.056465834  0.103646340 -0.395985411

##### validation of LASSO model ###### 
load('OP/Pipeline OP5/Train OP/TCGA_train_log2counts - lasso_model_min.Rdata')

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
# 0.8875994 0.9344898 0.8769280 - TCGA TRAINING ROC 
# 0.8817041 0.8775450 0.9107595 - TCGA VALIDATION ROC
# 0.6313236 0.8177750 0.8461348  - GLASS VALIDATION ROC
set.seed(100)
confint(ROC)
# $CI_AUC - TCGA TRAINING ROC 
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

# $CI_AUC - TCGA VALIDATION ROC
#         2.5% 97.5%
# t=365  81.42 94.92
# t=1095 80.33 95.18
# t=1825 84.70 97.45

# $CB_AUC
#         2.5% 97.5%
# t=365  80.22 96.13
# t=1095 79.00 96.51
# t=1825 83.57 98.59

# $C.alpha
#      95% 
# 2.310094 


# $CI_AUC  - GLASS VALIDATION ROC
#         2.5% 97.5%
# t=365  49.64 76.62
# t=1095 73.80 89.76
# t=1825 73.42 95.81

# $CB_AUC
#         2.5% 97.5%
# t=365  46.80 79.46
# t=1095 72.12 91.44
# t=1825 71.07 98.16

# $C.alpha
#      95% 
# 2.371966 


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

ggsave("OP/Pipeline OP5/GLASS Test OP/Images/GLASS_test_log2transform - ROC.tiff", plot = last_plot(), width = 900/72, height = 900/72, dpi = 300)


##### risk level and KMplot ######
#cutoff <- median(as.numeric(lasso.prob[ ,1])) # KEEP UNCHANGED. REMAINS CONSTANT BASED ON TCGA TRAINING SET
cutoff
# [1] 0.5370904 - old
# [1] -0.5888477 - new TCGA TRAINING SET

coxdata$Risk.Score.Level <- ifelse(coxdata$Risk.Score> cutoff, 'High', 'Low')

table(coxdata$Risk.Score.Level)
# High  Low 
#  233  233 - TCGA TRAINING SET
#  103  103 - TCGA TEST SET
#  124   29 - GLASS TEST SET

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
ggsave("OP/Pipeline OP5/GLASS Test OP/Images/GLASS_test_log2transform - KM Risk Score.tiff", plot = plot_combined, width = 900/72, height = 900/72, dpi = 300)


####### Risk score analysis, survival analysis and prognostic performance of a riskscore model based on differential expression genes #######
fp_dat <- data.frame(s = 1:length(lasso.prob[ ,1]), v = as.numeric(sort(lasso.prob[ ,1])))
fp_dat$`Risk Group` <- ifelse(fp_dat$v > cutoff, 'High', 'Low')

sur_dat <- data.frame(s = 1:length(lasso.prob[ ,1]),
                      t = coxdata[names(sort(lasso.prob[ ,1])), 'OS.time'],
                      OS = coxdata[names(sort(lasso.prob[ ,1])), 'OS'])

sur_dat$`Vital Status` <- ifelse(sur_dat$OS == 0, 'Alive', 'Death')
#sur_dat$`Vital Status`  <- factor(sur_dat$OS, levels = c("Death", "Alive"))
exp_dat <- coxdata[names(sort(lasso.prob[ ,1])), lasso_gene]

save(exp_dat, file="OP/Pipeline OP5/GLASS Test OP/GLASS_test_log2transform - exp_dat.RData")

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
ggsave("OP/Pipeline OP5/GLASS Test OP/Images/GLASS_test_log2transform - distribution and heatmap.tiff", plot = combined, width = 1100/72, height = 751/72, dpi = 300)


##### multiple COX regression ##### 
# prepare coxdata for merging - assign rownames as a column so you can use the merge() function
#coxdata <- tibble::rownames_to_column(coxdata, "Patient_ID")

# prepare pheno data so values in Patient ID column match the coxdata syntax
# replace - with . to match the colnames in the counts dataframe.
#tmp <- str_replace_all(glioma.pheno$Patient_ID, "-", ".")
#glioma.pheno$Patient_ID <- tmp

dim(glass.pheno)
# [1] 312   9

dim(coxdata)
# [1] 206  54 - TCGA TEST SET
# [1] 153  54 - GLASS TEST SET


# merge pheno data with coxdata
pheno.coxdata <- merge(coxdata, glass.pheno, by = 0, x.all=TRUE, y.all=FALSE)
dim(pheno.coxdata)
# [1] 153  64

names(pheno.coxdata)

pheno.coxdata <- pheno.coxdata[,-c(58:59)] # removes duplicate OS.y and OS.time_mos columns

new_dat_cox <- pheno.coxdata

new_dat_cox[new_dat_cox == "unknown"] <- NA

new_dat_cox <- droplevels.data.frame(new_dat_cox)

new_dat_cox$`Age (Years)` <- as.numeric(new_dat_cox$`Age (Years)`)

#new_dat_cox <- new_dat_cox %>% 
#  dplyr::rename(#`Age (Years)` = Age_years,
#                MGMTp.Status = MGMT.promoter.status.x,
#                TERTp.Status = TERT.promoter.status,
#                WHO.CNS5.Diagnosis = WHO_CNS5_diagnosis,
#                Histology = WHO_CNS5_histology,
#                Grade = WHO_CNS5_grade,
#                Codel.Status = WHO_CNS_1p.19q._codeletion,
#                IDHm.Status = WHO_CNS5_IDHstatus
#  )

# fix up table to have ages in years instead of day
new_dat_cox$Age <- ifelse(new_dat_cox$`Age (Years)` > 39, 'Over 39', 'Under 39')
new_dat_cox$age_level2 <- ifelse(new_dat_cox$Age == 'Over 39', '1', '0')

table(new_dat_cox$Age, useNA="ifany")
#  Over 39 Under 39     <NA> 
#.     267      159       39 - TCGA TRAIN
#      124       63       18 - TCGA TEST
#      126       27 - GLASS TEST

# get mean and SD of Age
#new_dat_cox$`Age (Years)` <- as.numeric(new_dat_cox$`Age (Years)`)
mean(new_dat_cox$`Age (Years)`, na.rm = TRUE)
# [1] 47.65728 - Train set
# [1] 47.54011 - TCGA TEST SET
# [1] 51.56863 - GLASS TEST SET
sd(new_dat_cox$`Age (Years)`, na.rm = TRUE)
# [1] 15.52152 - Train set
# [1] 14.61411 - TCGA TEST SET
# [1] 13.17881 - GLASS TEST SET

#new_dat_cox$Gender <- as.character(new_dat_cox$Gender) #1 female,2 male

table(new_dat_cox$Gender, useNA="ifany")
# female   male   <NA> 
#.  174    252     39 - TRAIN SET
#     83    104     18 - TCGA TEST SET
#     53    100 - GLASS TEST SET

#new_dat_cox$Grade <- as.character(new_dat_cox$WHO_CNS5_grade)

new_dat_cox$Grade[new_dat_cox$Grade == 'Unclassified'] <- NA
table(new_dat_cox$Grade, useNA="ifany")
# Grade 2 Grade 3 Grade 4    <NA> 
#.    132     109     166      58 - TRAIN SET
#      55      55      70      25 - TCGA TEST SET
#      15       9     128       1 - GLASS TEST SET

table(new_dat_cox$IDHm.Status, useNA="ifany")
#  IDH-mutant IDH-wildtype         <NA> 
#.        256          152           57 - TRAIN SET
#          120           61           24 - TCGA TEST SET
#           31          112           10 - GLASS TEST SET

new_dat_cox$MGMTp.Status <- ifelse(new_dat_cox$MGMTp.Status == 'M', 'Methylated',
                                   ifelse(new_dat_cox$MGMTp.Status == 'U', 'Unmethylated', 
                                          NA)
)

table(new_dat_cox$MGMTp.Status, useNA="ifany")
#   Methylated Unmethylated         <NA> 
#          326          112           27 - TRAIN SET
#          148           48            9 - TCGA TEST SET
#           41           40           72 - GLASS TEST SET

new_dat_cox$Codel.Status <- ifelse(new_dat_cox$Codel.Status == 'codel', 'Codel',
                                   ifelse(new_dat_cox$Codel.Status == 'non-codel', 'Non-codel', 
                                          NA)
)

table(new_dat_cox$Codel.Status, useNA="ifany")
# Codel Non-codel      <NA> 
# 117       344         4 - TRAIN SET
#  48       155         2 - TCGA TEST SET
#   9       128        16 - GLASS TEST SET

new_dat_cox$Histology[new_dat_cox$Histology == 'Unclassified'] <- NA

table(new_dat_cox$Histology, useNA="ifany")
#       astrocytoma      glioblastoma  oligoastrocytoma oligodendroglioma -OLD
#            115               115                76               121 - Train set
# Astrocytoma      Glioblastoma Oligodendroglioma              <NA> - NEW
#         150               152               106                57 - TRAIN SET
#          77                61                43                24 - TCGA TEST SET
#          52                84                 9                 8 - GLASS TEST SET


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
# 233 233 - TRAIN SET
# 102 103 - TCGA TEST SET
#  29 124 - GLASS TEST SET

# clean new_dat_cox df. remove unwanted columns.
#tmp <- new_dat_cox[, c(1:57, 71, 83:86, 89)]
#dim(tmp)
# [1] 465   63 - TRAIN SET
# [1] 205  63 - TCGA TEST SET

#new_dat_cox <- tmp
new_dat_cox <- column_to_rownames(new_dat_cox, var="Row.names")

save(new_dat_cox, file = "OP/Pipeline OP5/GLASS Test OP/GLASS_test_log2transform - new_dat_cox.RData")
load("OP/Pipeline OP5/GLASS Test OP/GLASS_test_log2transform - new_dat_cox.RData")
write.csv(new_dat_cox, 'OP/Pipeline OP5/GLASS Test OP/CSV/GLASS_test_log2transform - new_dat_cox.csv',row.names = FALSE)
###################### Get summary of clinical variables in new_dat_cox ###################### 
lapply(new_dat_cox[,54:62], table, useNA = "always")

new_dat_cox %>%
  group_by(Histology, Grade) %>%
  summarise(count = n())

####----VII. Validate Nomogram Model----####

nom_dat_cox <- new_dat_cox

print(sapply(nom_dat_cox, function(x) sum(is.na(x)))) # Check for NAs

# DO TTHIS FOR VALIDATION ####
# 1-year calibration curve ####
coxm_1 <- cph(Surv(OS.time,OS.x) ~ Grade + Age + Risk.Score.Level,
              data = nom_dat_cox,
              surv = T,x = T,y = T,
              time.inc = 365)

cal_1 <- calibrate(coxm_1,
                   u = 365,
                   cmethod = 'KM',
                   m = 10, # previously 60
                   B = 500) # previously 1000

#par(mar = c(7,4,4,3), cex = 1.0)

# 3-year calibration curve ####
coxm_3 <- cph(Surv(OS.time,OS.x) ~ Grade + Age + Risk.Score.Level,
              data = nom_dat_cox,
              surv = T, x = T, y = T,
              time.inc = 365*3)

cal_3 <- calibrate(coxm_3, u = 365 * 3, cmethod = 'KM', m = 10, B = 200)

# 5-year calibration curve ####
coxm_5 <- cph(Surv(OS.time,OS.x) ~ Grade + Age + Risk.Score.Level,
              data = nom_dat_cox,
              surv = T, x = T, y = T,
              time.inc = 365 * 5)

cal_5 <- calibrate(coxm_5, u = 365 * 5, cmethod = 'KM', m = 10, B = 200)

# plot calibration curves ####

# open a TIFF device to save the plot
tiff(filename = "OP/Pipeline OP5/GLASS Test OP/Images/GLASS_test_log2transform - Nomogram Calibration2.tiff",
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
f <- coxph(Surv(OS.time,OS.x) ~ Grade + Age + Risk.Score.Level, data = nom_dat_cox)
sum.surv <- summary(f)
c_index <- sum.surv$concordance
c_index
#          C      se(C) 
# 0.82502128 0.01291766 - TRAIN SET
# 0.83902401 0.01643817 - TCGA TEST SET
# 0.6768994 0.0225137 - GLASS TEST SET

# timeROC for nomogram model
pred_nom <- predict(f, newdata = nom_dat_cox)
nom_dat_cox$pred_nom <- as.numeric(pred_nom)

with(nom_dat_cox,
     ROC <<- timeROC(T = OS.time,
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
# 0.8515814 0.9498082 0.8646947  - TRAIN SET
# 0.8721805 0.9080091 0.9634069 - TCGA TEST SET
# 0.6211064 0.8061585 0.8940603 - GLASS TEST SET

confint(ROC)
# $CI_AUC - TRAIN SET
#         2.5% 97.5%
# t=365  81.75 88.57
# t=1095 92.11 97.85
# t=1825 80.04 92.90

# $CB_AUC
#.        2.5% 97.5%
# t=365  80.92 89.39
# t=1095 91.42 98.54
# t=1825 78.48 94.46

# $C.alpha
#      95% 
# 2.434855  

# $CI_AUC - TCGA TEST SET
#         2.5%  97.5%
# t=365  83.50  90.94
# t=1095 84.88  96.72
# t=1825 92.49 100.20

# $CB_AUC
#         2.5%  97.5%
# t=365  82.71  91.72
# t=1095 83.63  97.97
# t=1825 91.67 101.01

# $C.alpha
#    95% 
# 2.3747 

# $CI_AUC - GLASS TEST SET
#         2.5% 97.5%
# t=365  54.40 69.83
# t=1095 73.17 88.06
# t=1825 80.17 98.64

# $CB_AUC
#         2.5%  97.5%
# t=365  52.78  71.45
# t=1095 71.61  89.62
# t=1825 78.23 100.58

# $C.alpha
#      95% 
# 2.371474 

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

ggsave("OP/Pipeline OP5/GLASS Test OP/Images/GLASS_test_log2transform - Nomogram ROC.tiff", plot = last_plot(), width = 900/72, height = 751/72, dpi = 300)

new_glass_nom <- nom_dat_cox

save(new_glass_nom, file = 'OP/Pipeline OP5/GLASS Test OP/GLASS_test_log2transform - new_glass_nom.Rdata')

#### ------------------- MODIFY HEATMAP TO INCLUDE CLINICAL ANNOTATIONS -------------------------------####

library(pheatmap)

tmp = t(scale(exp_dat))
tmp[tmp > 1] = 1
tmp[tmp < -1] = -1

df <- as.data.frame(new_glass_nom) #log2(counts+1)

names(df)

# Fix column names for annotation - troubleshooting problem that annotation does not appear above heatmap
# rename column names in prep for figure
df <- df %>% 
  dplyr::rename(`Risk Group` = Risk.Score.Level,
                `Risk Score` = Risk.Score,
                `MGMTp Status` = MGMTp.Status,
                `1p/19q Codel Status` = Codel.Status,
                #`WHO CNS5 Classification` = WHO_CNS_diagnosis.., # fix for Tcga set
                #`WHO CNS5 Classification` = WHO_CNS5_diagnosis, # fix for cgga sett
                `IDHm Status` = IDHm.Status
  )

str(df)

# convert CGGA_ID column to rownames####
#df <- column_to_rownames(df, var = "CGGA_ID")

# remove pred_nom column
#df <- df[,-60]

# This step is only needed if column have factor levels (for training set. Not needed for Validation set) ####
df$Histology <- droplevels(df$Histology)
df$Histology <- as.character(df$Histology)

df$`IDHm Status` <- droplevels(df$`IDHm Status`)
df$`IDHm Status` <- as.character(df$`IDHm Status`)

df$`MGMTp Status` <- droplevels(df$`MGMTp Status`)
df$`MGMTp Status` <- as.character(df$`MGMTp Status`)

df$`1p/19q Codel Status` <- droplevels(df$`1p/19q Codel Status`)
df$`1p/19q Codel Status` <- as.character(df$`1p/19q Codel Status`)

str(df)

# NOTE: Each data structure should be characters. for code below to work.

#### Generate column annotations ####
annotation <- df[ , c('OS.x', 'Grade', 'Histology', 'Gender', 'Age', "MGMTp Status", '1p/19q Codel Status', 'IDHm Status', 'Risk Group')] 

# Change "OS" [0,1] to "Vital Status" [Alive, Dead]####
annotation$`Vital Status` <- ifelse(annotation$OS.x == 0, 'Alive', 'Dead')
annotation <- annotation %>% dplyr::select(`Vital Status`, everything())
annotation <- annotation[,-2] # Remove OS so it doesn't appear in image.


ann_colors <- list(
  `Risk Group` = c(High = "#ff0000",
                   Low = "#00468B"),
  `IDHm Status` = c("IDH-wildtype" = "#fdaf91", 
                    "IDH-mutant" = "#AD002A"),
  `1p/19q Codeletion Status` = c(Codel = "#0099B4FF", 
                                 "Non-codel" = "#1B191999"),
  `MGMTp Status` = c(Methylated = "#AD002A99", 
                     Unmethylated = "#0099B4FF"),
  Age = c("Over 39" = "#7E6148B2", 
          "Under 39" = "#91D1C2B2"),
  Gender = c(female = "#E64B35B2", 
             male = "#3C5488B2"),
  Histology = c(Astrocytoma = "#0099B499", 
                Glioblastoma = "#AD002A99", 
                Oligodendroglioma = "#925E9F99"),
  Grade = c("Grade 2" = "#00468B99", 
            "Grade 3" = "#ED000099", 
            "Grade 4" = "#42B54099"),
  `Vital Status` = c(Alive = "#42b540",
                     Dead = "#925e9f")
)

# troubleshooting annotations again! ####
#annotation <- df[ , c("OS", 'Grade', 'Histology', "Gender", "Age", "MGMTp Status", "1p/19q Codel Status", "IDHm Status", "Risk Score Group")]

# Change "OS" [0,1] to "Vital Status" [Alive, Dead]
#annotation$`Vital Status` <- ifelse(annotation$OS == 0, 'Alive', 'Dead')
#annotation <- annotation %>% select(`Vital Status`, everything())
#annotation <- annotation[,-2] # Remove OS so it doesn't appear in image.

ann_colors <- list(
  `Risk Group` = c(High = "#ff0000",
                   Low = "#00468B"),
  `IDHm Status` = c("IDH-wildtype" = "#fdaf91", 
                    "IDH-mutant" = "#AD002A"),
  `1p/19q Codel Status` = c(Codel = "#0099B4FF", 
                            "Non-codel" = "#1B191999"),
  `MGMTp Status` = c(Methylated = "#AD002A99", 
                     Unmethylated = "#0099B4FF"),
  Age = c("Over 39" = "#7E6148B2", 
          "Under 39" = "#91D1C2B2"),
  Gender = c(female = "#E64B35B2", 
             male = "#3C5488B2"),
  Histology = c(Astrocytoma = "#0099B499", 
                Glioblastoma = "#AD002A99", 
                Oligodendroglioma = "#925E9F99"),
  Grade = c("Grade 2" = "#00468B99", 
            "Grade 3" = "#ED000099", 
            "Grade 4" = "#42B54099"),
  `Vital Status` = c(Alive = "#42b540",
                     Dead = "#925e9f")
)

# plot heatmap ####
heatmap_plot <-
  pheatmap(tmp, 
           #main = "Training Set (TCGA)",
           main = "TCGA Test Set",
           col = colorRampPalette(c("blue", "white", "red"))(50),
           #annotation_row = gene_functions_df,
           annotation_col = annotation, 
           annotation_colors = ann_colors, 
           cluster_cols = F,
           show_colnames = F,
           show_rownames = T,
           fontsize = 14,
           treeheight_row = 0, # remove dendoogram from rows of heatmap
           #fontsize_col = 20,
           cluster_rows = F,
           #cutree_rows = 2,
           
           # Adjust cell size to increase overall plot size
           #cellwidth = 1,  # Increase this value to make the plot wider
           #cellheight = 10, # Increase this value to make the plot taller
           
           # Add margin to the plot
           margins = c(10, 10, 10, 10),  # top, right, bottom, left margins
           
           # Optional: Add gaps between groups of rows or columns
           # gaps_row = c(10, 20),  # Add gaps after row 10 and 20
           # gaps_col = c(50, 100), # Add gaps after column 50 and 100
           #gaps_row = train_row_order  # Use the extracted row order
           
           # Adjust legend position
           #legend_labels = c("Low", "Medium", "High"),
           #legend = TRUE,
           #legend_breaks = c(-2, 0, 2),
           #legend_position = "bottom"
  )

print(heatmap_plot)

# open a TIFF device to save the plot
tiff(filename = "OP/Pipeline OP5/GLASS Test OP/Images/GLASS_test_log2transform - clinical heatmap output no cluster.tiff",
     width = 12, height = 10, units = "in", res = 300)
print(heatmap_plot)
dev.off()


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

# Invesitigate why there's low AUC in GLASS 1-year predictioon ---------------------------
# Assuming your dataframe is named 'df' and the column is named 'col'
# If different, replace 'df' with your dataframe name and 'col' with your column name

count_365_or_less <- sum(new_dat_cox$OS.time <= 365, na.rm = TRUE)
count_366_to_1095 <- sum(new_dat_cox$OS.time > 365 & new_dat_cox$OS.time <= 1095, na.rm = TRUE)
count_1096_or_more <- sum(new_dat_cox$OS.time >= 1096, na.rm = TRUE)

# Print the results
cat("Rows with col <= 365:", count_365_or_less, "\n")
cat("Rows with 366 <= col <= 1095:", count_366_to_1095, "\n")
cat("Rows with col >= 1096:", count_1096_or_more, "\n")

total_non_na <- sum(!is.na(df$col))
cat("Percentage <= 365:", round(count_365_or_less / total_non_na * 100, 2), "%\n")
cat("Percentage 366-1095:", round(count_366_to_1095 / total_non_na * 100, 2), "%\n")
cat("Percentage >= 1096:", round(count_1096_or_more / total_non_na * 100, 2), "%\n")


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

