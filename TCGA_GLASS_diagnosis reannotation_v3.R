# Required packages
library(TCGAbiolinks)  # For TCGA data queries
library(dplyr)         # For data manipulation
library(data.table)    # For data manipulation
library(vtable)        # For summary statistics
library(tidyverse)     # For data wrangling
library(reshape2)      # For restructuring and transforming data 
library(readr)         # For fast and efficient reading of text data files
library(stringr)       # For string manipulation


####--- TASK: Re-classify TCGA dataset to WHO CNS5 nomenclature ----####
# Reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9820617/

# Fetch molecular subtypes for GBM and LGG
gbm.path.subtypes <- TCGAquery_subtype(tumor = "gbm")
lgg.path.subtypes <- TCGAquery_subtype(tumor = "lgg")

# Convert factor columns to character
i <- sapply(gbm.path.subtypes, is.factor)
gbm.path.subtypes[i] <- lapply(gbm.path.subtypes[i], as.character)
i <- sapply(lgg.path.subtypes, is.factor)
lgg.path.subtypes[i] <- lapply(lgg.path.subtypes[i], as.character)

# Combine GBM and LGG phenotype data
glioma.pheno <- rbind(gbm.path.subtypes, lgg.path.subtypes)

# Generate summary statistics
st(glioma.pheno, out='csv', file='OP/full_TCGAorig_pheno_summary_v3.csv')

# Read WHO CNS5 classification data
fn4 = "IP/TCGAreclass_ijms-2057006_input.csv"  # Update path as needed
WHO.CNS5 <- read.csv(fn4, header = TRUE)
st(WHO.CNS5, out='csv', file='OP/full_TCGAijms_pheno_summary_v3.csv')

# Merge original TCGA annotations with WHO CNS5 classifications
TCGA.WHO.CNS5 <- merge(glioma.pheno, WHO.CNS5, by.x = "patient", by.y = "Patient_ID")

# Read survival data from UCSC Xena
gbm.surv.file <- "IP/TCGA-GBM.survival.tsv"  # Update path as needed
lgg.surv.file <- "IP/TCGA-LGG.survival.tsv"  # Update path as needed

gbm.surv <- read.table(gbm.surv.file, header = TRUE, stringsAsFactors = FALSE)
lgg.surv <- read.table(lgg.surv.file, header = TRUE, stringsAsFactors = FALSE)

# Combine GBM and LGG survival data
glioma.surv <- rbind(gbm.surv, lgg.surv)
glioma.surv2 <- glioma.surv
glioma.surv2$patient <- glioma.surv2$X_PATIENT
glioma.surv2 <- glioma.surv2[,c(5,4,2)]

# Function to identify duplicates
clean_get_dups <- function(df = NULL, var = '', freq = 1) {
  n_occur <- data.frame(table(df[[var]]))
  dups <- n_occur$Freq > freq
  df_dups <- df[df[[var]] %in% n_occur$Var1[dups], ]
  return(df_dups)
}

# Check for duplicates
dup_check <- clean_get_dups(glioma.surv2, "patient", 1)

# Remove duplicates, keeping first occurrence
glioma.surv2 <- glioma.surv2[!duplicated(glioma.surv2$patient), ]

# Merge survival data with annotation data
TCGA.WHO.CNS52 <- TCGA.WHO.CNS5 %>% 
  left_join(glioma.surv2, by = join_by(patient))

# Define columns to keep
cols_to_keep <- c(
  # Patient Information
  "patient", "OS", "OS.time", "Gender", "Age..years.at.diagnosis..x", 
  "Karnofsky.Performance.Score",
  
  # Molecular Characteristics
  "ABSOLUTE.purity", "ESTIMATE.stromal.score", "ESTIMATE.immune.score",
  "ESTIMATE.combined.score", "Transcriptome.Subtype", "Mutation.Count",
  "Percent.aneuploidy",
  
  # Clinical Classifications
  "TCGA.Histology", "TCGA.Grade",
  
  # Molecular Markers
  "IDH.status.x", "X1p.19q.codeletion.x", "IDH.codel.subtype",
  "MGMT.promoter.status.x", "Chr.7.gain.Chr.10.loss.x",
  "TERT.promoter.status.x", "TERT.expression.status",
  "ATRX.status.x", "BRAF.V600E.status", "BRAF.KIAA1549.fusion",
  
  # Copy Number Variations
  "EGFR_CNV", "CDKN2A_CNV", "CDKN2B_CNV",
  
  # WHO CNS5 Classifications
  "H3.K27.G34.status", "WHO_CNS5_diagnosis", "WHO_CNS5_histology",
  "WHO_CNS5_grade", "WHO_CNS5_IDHstatus", "WHO_CNS_1p.19q._codeletion",
  "WHO_CNS5_diagnosis.original.ijms.2057006"
)

# Extract relevant columns and rename them
clean.WHO.CNS5_v2 <- TCGA.WHO.CNS52[, cols_to_keep] %>%
  dplyr::rename(
    Patient_ID = patient,
    Age_years = Age..years.at.diagnosis..x,
    KPS = Karnofsky.Performance.Score,
    IDH.status = IDH.status.x,
    Codel = X1p.19q.codeletion.x,
    TERT.promoter.status = TERT.promoter.status.x,
    ATRX.status = ATRX.status.x,
    Chr.7.gain.Chr.10.loss = Chr.7.gain.Chr.10.loss.x
  )

# Handle missing values and clean up data
clean.WHO.CNS5_v2[is.na(clean.WHO.CNS5_v2)] <- 'unknown'
clean.WHO.CNS5_v2 <- clean.WHO.CNS5_v2[!(clean.WHO.CNS5_v2$OS == "unknown" & 
                                           clean.WHO.CNS5_v2$OS.time == "unknown"), ]

# Create WHO CNS diagnosis classifications
clean.WHO.CNS5_v2$WHO_CNS_diagnosis.. <- case_when(
  clean.WHO.CNS5_v2$WHO_CNS5_diagnosis %in% c(
    "Astrocytoma, IDH-mutant. Grade 2.",
    "Astrocytoma, IDH-mutant. Grade 3.",
    "Astrocytoma, IDH-mutant. Grade 4.",
    "Astrocytoma, IDH-mutant. Grade NA."
  ) ~ "Astrocytoma, IDH-mutant",
  
  clean.WHO.CNS5_v2$WHO_CNS5_diagnosis == "Glioblastoma, IDH-wildtype. Grade 4." ~ 
    "Glioblastoma, IDH-wildtype",
  
  clean.WHO.CNS5_v2$WHO_CNS5_diagnosis %in% c(
    "Oligodendroglioma, IDH-mutant, 1p/19q-codeleted. Grade 2.",
    "Oligodendroglioma, IDH-mutant, 1p/19q-codeleted. Grade 3."
  ) ~ "Oligodendroglioma, IDH-mutant, 1p/19q-codeleted",
  
  TRUE ~ "Unclassified"
)

# Remove patients under 18 years old
clean.WHO.CNS5_v2 <- clean.WHO.CNS5_v2[clean.WHO.CNS5_v2$Age_years >= 18, ]

# Save outputs
write.csv(clean.WHO.CNS5_v2, "OP/clean_TCGA.WHO-CNS5_v3.csv", row.names = FALSE)
save(clean.WHO.CNS5_v2, file = "OP/clean_TCGA.WHO-CNS5_v3.RData")

####--- TASK: Re-classify GLASS dataset to WHO CNS5 nomenclature ----####

# Import GLASS Clinical Surgeries data
fn5 <- "IP/GLASS_clinical_surgeries.csv"
GLASS_pheno <- read.csv(fn5, header = TRUE)
GLASS_pheno <- GLASS_pheno[, -c(1:2)]

# Import GLASS Clinical Cases data
fn6 <- "IP/GLASS_clinical_cases.csv"
GLASS_cases <- read.csv(fn6, header = TRUE)
GLASS_cases <- GLASS_cases[, -c(1:2)]

# Merge Clinical Surgeries and Cases data
glass.pheno.cases_merged <- merge(GLASS_pheno, GLASS_cases, all = TRUE)
glass.pheno.cases_merged[is.na(glass.pheno.cases_merged)] <- 'unknown'

# Remove TCGA samples
glass.pheno.cases_merged <- glass.pheno.cases_merged[!grepl("TCGA", glass.pheno.cases_merged$case_project), ]

# Rename columns to match TCGA nomenclature
glass.pheno.cases_merged <- glass.pheno.cases_merged %>%
  rename(
    TCGA.Histology = histology,
    TCGA.Grade = grade,
    IDH.status = idh_status,
    Codel = codel_status,
    MGMT.promoter.status = mgmt_methylation,
    Gender = case_sex,
    Age_years = case_age_diagnosis_years,
    Vital_status = case_vital_status,
    Survival_Months = case_overall_survival_mo
  )

# Initialize GLASS CNS5 dataset
GLASS_CNS5 <- glass.pheno.cases_merged

# Process IDH status annotations
names(GLASS_CNS5)[7] <- "IDH.status"
GLASS_CNS5$WHO_CNS5_IDHstatus <- case_when(
  GLASS_CNS5$IDH.status == 'IDHmut' ~ 'IDH-mutant',
  GLASS_CNS5$IDH.status == 'IDHwt' ~ 'IDH-wildtype',
  TRUE ~ 'unknown'
)

# Process Codeletion status
names(GLASS_CNS5)[8] <- "Codel"
GLASS_CNS5 <- GLASS_CNS5 %>%
  mutate(Codel = case_when(
    Codel == "noncodel" ~ "non-codel",
    TRUE ~ Codel
  ))

# Create IDH-Codel subtype classification
GLASS_CNS5 <- GLASS_CNS5 %>%
  mutate(IDH.codel.subtype = case_when(
    IDH.status == "IDHmut" & Codel == "codel" ~ "IDHmut-codel",
    IDH.status == "IDHmut" & Codel == "non-codel" ~ "IDHmut-non-codel",
    IDH.status == "IDHwt" & Codel == "non-codel" ~ "IDHwt",
    IDH.status == "IDHmut" & Codel == "unknown" ~ "Unclassified",
    IDH.status == "IDHwt" & Codel == "unknown" ~ "IDHwt",
    IDH.status == "unknown" & Codel == "codel" ~ "Unclassified",
    IDH.status == "unknown" & Codel == "non-codel" ~ "Unclassified",
    IDH.status == "unknown" & Codel == "unknown" ~ "unknown",
    TRUE ~ "Check"
  ))

# Save processed phenotype data
save(GLASS_CNS5, file = "OP/GLASS.pheno_v3.RData")

# Filter for primary tumor samples
GLASS_CNS5.TP <- GLASS_CNS5[grepl("-TP", GLASS_CNS5$sample_barcode), ]
save(GLASS_CNS5.TP, file = "OP/GLASS.TP.pheno_v3.RData")

# Process variant annotations
glass.anno <- read_csv(
  "IP/GLASS_variants_anno_20220531.csv.gz",
  show_col_types = FALSE
)

glass.passgeno <- read_csv(
  "IP/GLASS_variants_passgeno_20220531.csv.gz",
  show_col_types = FALSE
)

# Merge variant annotations with genotype data
glass.passgeno.anno <- merge(
  glass.passgeno, 
  glass.anno, 
  by = c("chrom", "start", "end", "alt"), 
  all.x = TRUE
)

# Add project information and filter TCGA samples
glass.passgeno.anno$case_project <- substr(glass.passgeno.anno$aliquot_barcode, 1, 4)
glass.passgeno.anno <- glass.passgeno.anno[!grepl("TCGA", glass.passgeno.anno$case_project), ]

# Save processed variant data
save(glass.passgeno.anno, file = "OP/GLASS.var.passgeno.anno_v3.RData")

# Filter primary tumor samples from variant data
glass.passgeno.anno.TP <- glass.passgeno.anno[grepl("-TP-", glass.passgeno.anno$aliquot_barcode), ]

# Create sample identification columns
glass.passgeno.anno.TP$case_barcode <- substr(glass.passgeno.anno.TP$aliquot_barcode, 1, 12)
glass.passgeno.anno.TP$sample_barcode <- substr(glass.passgeno.anno.TP$aliquot_barcode, 1, 15)

# Reorganize columns for better readability
glass.passgeno.anno.TP <- glass.passgeno.anno.TP %>% 
  select(c(5, 27, 26, 25, 11), everything())

save(glass.passgeno.anno.TP, file = "OP/GLASS.TP.var.passgeno.anno_v3.RData")

# Function to process mutation data
process_mutation_data <- function(data, gene, mutation_pattern, status_column) {
  # Filter for specific gene
  filtered_data <- data[data$gene_symbol == gene, ]
  filtered_data <- filtered_data[!is.na(filtered_data$sample_barcode), ]
  
  # Create mutation status column using with() to properly evaluate in data context
  filtered_data[[status_column]] <- with(filtered_data, 
                                         ifelse(
                                           case_when(
                                             # TERT promoter mutations
                                             gene == "TERT" ~ genome_change %in% c("g.chr5:1295228C>T", "g.chr5:1295250C>T"),
                                             # BRAF V600E mutation
                                             gene == "BRAF" ~ grepl("p.V600E", protein_change),
                                             # ATRX mutations
                                             gene == "ATRX" ~ grepl("g.chrX", genome_change),
                                             # IDH1 mutations
                                             gene == "IDH1" ~ protein_change %in% c("p.R132H", "p.R132C", "p.R132G", "p.R132S"),
                                             # IDH2 mutations
                                             gene == "IDH2" ~ protein_change == "p.R172K",
                                             # H3 mutations
                                             gene %in% c("H3F3A", "HIST1H3B", "HIST1H3C") ~ protein_change %in% c("p.K27M", "p.G34R"),
                                             # Default case
                                             TRUE ~ FALSE
                                           ),
                                           "Mutant",
                                           "Wildtype"
                                         )
  )
  
  # Create unique entries per sample
  unique_data <- filtered_data %>%
    group_by(sample_barcode) %>%
    summarize(
      !!status_column := case_when(
        any(!!sym(status_column) == "Mutant") ~ "Mutant",
        n_distinct(!!sym(status_column)) == 1 ~ first(!!sym(status_column)),
        TRUE ~ first(!!sym(status_column))
      )
    ) %>%
    ungroup()
  
  return(list(full = filtered_data, unique = unique_data))
}

# Process TERT mutations
tert_data <- process_mutation_data(
  glass.passgeno.anno.TP,
  'TERT',
  NULL,  # Not used anymore
  'TERTp.C228T.C250T.mutation'
)

# Process BRAF mutations
braf_data <- process_mutation_data(
  glass.passgeno.anno.TP,
  'BRAF',
  NULL,
  'BRAF.V600E.status'
)

# Process ATRX mutations
atrx_data <- process_mutation_data(
  glass.passgeno.anno.TP,
  'ATRX',
  NULL,
  'ATRX.status'
)

# Process IDH1 mutations
idh1_data <- process_mutation_data(
  glass.passgeno.anno.TP,
  'IDH1',
  NULL,
  'IDH1.R132seq.status'
)

# Process IDH2 mutations
idh2_data <- process_mutation_data(
  glass.passgeno.anno.TP,
  'IDH2',
  NULL,
  'IDH2.R172seq.status'
)

# Process H3 mutations
h3_data <- process_mutation_data(
  glass.passgeno.anno.TP,
  'H3F3A',  # Process each H3 gene separately if needed
  NULL,
  'H3.K27M.G34R.status'
)

# Combine mutation data
mutation_summary <- Reduce(function(x, y) {
  full_join(x, y, by = "sample_barcode")
}, list(
  tert_data$unique,
  braf_data$unique,
  atrx_data$unique,
  idh1_data$unique,
  idh2_data$unique,
  h3_data$unique
))

save(mutation_summary, file = "OP/GLASS.mutation.summary_v3.RData")

# Import gene copy number data
glass.geneCN <- read_csv(
  "IP/GLASS_variants_gene_copy_number.csv.gz",
  show_col_types = FALSE
)

# Create project identifier and remove TCGA samples
glass.geneCN$case_project <- substr(glass.geneCN$aliquot_barcode, 1, 4)
glass.geneCN <- glass.geneCN[!grepl("TCGA", glass.geneCN$case_project), ]
save(glass.geneCN, file = "OP/GLASS.var.geneCN_v3.RData")

# Filter for primary tumor samples and add identifiers
glass.geneCN.TP <- glass.geneCN[grepl("-TP-", glass.geneCN$aliquot_barcode), ] %>%
  mutate(
    case_barcode = substr(aliquot_barcode, 1, 12),
    sample_barcode = substr(aliquot_barcode, 1, 15)
  ) %>%
  select(aliquot_barcode, case_barcode, sample_barcode, case_project, everything())

save(glass.geneCN.TP, file = "OP/GLASS.TP.var.geneCN_v3.RData")

# Reshape data to wide format for specific genes
glass.geneCN.TP.wide <- reshape2::dcast(
  glass.geneCN.TP, 
  aliquot_barcode ~ gene_symbol, 
  value.var = "hlvl_call"
)

# Extract copy number data for specific genes
target_genes <- c("aliquot_barcode", "EGFR", "CDKN2A", "CDKN2B")
glass.geneCN.processed <- glass.geneCN.TP.wide[, target_genes] %>%
  mutate(
    case_barcode = substr(aliquot_barcode, 1, 12),
    sample_barcode = substr(aliquot_barcode, 1, 15),
    case_project = substr(aliquot_barcode, 1, 4)
  ) %>%
  select(aliquot_barcode, sample_barcode, case_barcode, case_project, everything())

# Remove duplicate samples, keeping first occurrence
glass.geneCN.final <- glass.geneCN.processed[!duplicated(glass.geneCN.processed$case_barcode), ] %>%
  select(-c(aliquot_barcode, case_barcode, case_project))

# Save processed copy number data
save(glass.geneCN.final, file = "OP/GLASS.TP.var.geneCN.small_df_v3.RData")

# Data dictionary for copy number values:
# -2 = likely homozygous (high-level) deletion
# -1 = likely deletion
#  0 = likely copy neutral
#  1 = likely amplification
#  2 = likely high-level amplification

#### investigate whether a sample in GLASS has chr7 gain and/or chr10 loss to match clean.WHO.CNS5 annotations ####
#### GLASS titan_seg ####
# TITAN segmentation output after using GATK copy number pipeline as input
# Each row represents a segment

# Import TITAN segmentation data
glass.titanseg <- read_csv("IP/GLASS_variants_titan_seg.csv", show_col_types = FALSE) %>%
  mutate(
    case_project = substr(pair_barcode, 1, 4),
    case_barcode = substr(pair_barcode, 1, 12),
    sample_barcode = substr(pair_barcode, 1, 15)
  )

# Filter GLASS project samples
glass.titanseg <- glass.titanseg[!grepl("TCGA", glass.titanseg$case_project), ]
save(glass.titanseg, file = "OP/GLASS.titan_segmentation_v3.RData")

# Filter primary tumor samples
glass.titanseg.TP <- glass.titanseg[grepl("-TP", glass.titanseg$sample_barcode), ]
save(glass.titanseg.TP, file = "OP/GLASS.TP.titan_segmentation_v3.RData")

# Function to classify CNV status
classify_cnv <- function(call) {
  case_when(
    call %in% c("HOMD", "DLOH") ~ "Loss",
    call %in% c("ALOH", "GAIN", "ASCNA", "BCNA", "UBCNA", "HLAMP") ~ "Gain",
    TRUE ~ "Neutral"
  )
}

# Process segmentation data with simplified CNV annotations
glass.titanseg.TP.CNsimple <- glass.titanseg.TP %>%
  mutate(call_simple = sapply(corrected_call, classify_cnv)) %>%
  group_by(sample_barcode, chrom) %>%
  summarise(
    median_copy_number = round(median(corrected_copy_number)),
    call_simple = names(which.max(table(call_simple)))
  ) %>%
  ungroup() %>%
  mutate(
    median_call = case_when(
      median_copy_number <= 1 ~ "loss",
      median_copy_number >= 3 ~ "gain",
      median_copy_number == 2 ~ "neutral"
    )
  )

# Process chromosome 7 gain and 10 loss status
glass.Chr7.Chr10.CNcall <- glass.titanseg.TP.CNsimple %>%
  filter(chrom %in% c("7", "10")) %>%
  pivot_wider(
    id_cols = sample_barcode,
    names_from = chrom,
    values_from = call_simple,
    names_prefix = "chr_"
  ) %>%
  mutate(
    Chr.7.gain.Chr.10.loss = case_when(
      chr_7 == "Gain" & chr_10 == "Loss" ~ "Gain chr 7 & loss chr 10",
      chr_7 %in% c("Gain", "Loss", "Neutral") & 
        chr_10 %in% c("Gain", "Loss", "Neutral") ~ "No combined CNA",
      TRUE ~ "unknown"
    )
  ) %>%
  select(sample_barcode, Chr.7.gain.Chr.10.loss)

save(glass.Chr7.Chr10.CNcall, file = "OP/GLASS.TP.Chr7-gain.Chr10loss_v3.RData")

#### 3.1f - MERGE all dataframes together ####
GLASS_CNS5.TP # main pheno df - merged GLASS clinical cases and clinical surgeries
glass.Chr7.Chr10.CNcall # chromosomal 7 gain and 10 loss - GLASS titan_seg df
glass.geneCN.final # copy number hlvl calls for EGFR, CDKN2A, CDKN2B - GLASS gene_copy_number df
atrx_data$unique # ATRX mutation - filtered from GLASS passgeno and anno dfs
braf_data$unique # BRAF mutation - filtered from GLASS passgeno and anno dfs
tert_data$unique # TERTp (C228T and C250T) mutation -  filtered from GLASS passgeno and anno dfs
h3_data$unique # H3.K27M.G34R mutations -  filtered from GLASS passgeno and anno dfs
idh1_data$unique # IDH1 (R132) mutations -  filtered from GLASS passgeno and anno dfs
idh2_data$unique # IDH2 (R172) mutations -  filtered from GLASS passgeno and anno dfs

# Define column groupings for final dataset organization
id = c("case_barcode", "sample_barcode", "case_project", "case_source")
clinical = c("Gender", "Age_years", "Vital_status", 
               "Survival_Months")
surgical = c("surgery_number", "surgical_interval_mo", "surgery_type", 
               "surgery_indication", "surgery_extent_of_resection", 
               "surgery_laterality", "surgery_location")
treatment = c("treatment_tmz", "treatment_tmz_cycles", "treatment_tmz_cycles_6",
                "treatment_tmz_cycles_notes", "treatment_concurrent_tmz",
                "treatment_radiotherapy", "treatment_radiation_dose_gy",
                "treatment_radiation_fractions", "treatment_radiation_other",
                "treatment_chemotherapy_other", "treatment_chemotherapy_other_cycles",
                "comments", "treatment_alkylating_agent")
who2007 = c("TCGA.Histology", "TCGA.Grade", "IDH.status", "Codel", "who_classification",
              "MGMT.promoter.status", "mgmt_methylation_method", "idh_codel_subtype")
molecular = c("Chr.7.gain.Chr.10.loss", "EGFR", "CDKN2A", "CDKN2B",
                "IDH1.R132seq.status", "IDH2.R172seq.status", "BRAF.V600E.status",
                "ATRX.status", "TERTp.C228T.C250T.mutation", "H3.K27M.G34R.status")
who2021 = c("WHO_CNS5_IDHstatus", "IDH.codel.subtype")

# Combine all lists in the desired order
all_cols <- c(id, clinical, surgical, treatment, who2007, molecular, who2021)

# Merge all datasets
integrated_glass_data <- GLASS_CNS5.TP %>%
  left_join(glass.Chr7.Chr10.CNcall, by = "sample_barcode") %>%
  left_join(glass.geneCN.final, by = "sample_barcode") %>%
  left_join(atrx_data$unique, by = "sample_barcode") %>%
  left_join(braf_data$unique, by = "sample_barcode") %>%
  left_join(tert_data$unique, by = "sample_barcode") %>%
  left_join(h3_data$unique, by = "sample_barcode") %>%
  left_join(idh1_data$unique, by = "sample_barcode") %>%
  left_join(idh2_data$unique, by = "sample_barcode")

#Rearrange columns and keep any additional columns at the end
integrated_glass_data_arranged <- integrated_glass_data[, all_cols]

# Save integrated dataset
save(integrated_glass_data_arranged, file = "OP/GLASS.integrated.dataset.RData")

# Generate summary statistics
data_summary <- integrated_glass_data_arranged %>%
  summarise(
    total_samples = n(),
    samples_with_molecular = sum(!is.na(Chr.7.gain.Chr.10.loss) | 
                                   !is.na(EGFR) | 
                                   !is.na(IDH1.R132seq.status)),
    samples_with_clinical = sum(!is.na(Age_years) & 
                                  !is.na(Vital_status)),
    samples_with_treatment = sum(!is.na(treatment_tmz) | 
                                   !is.na(treatment_radiotherapy))
  )

write_csv(data_summary, "OP/GLASS.data.summary.csv")

#### CONCLUSION ####

integrated_glass_data_arranged # Has all variables needed for re-annotation.

####---- TASK: reannotate GLASS dataset based on WHO.CNS5 nomenclature. ----####
clean.GLASS.CNS5.TP <- integrated_glass_data_arranged# main df for re-annotation

#### Subtasks:
#### 4.1. Remove samples less than 18 years old. Keep only patient >=18 years at diagnosis.

clean.GLASS.CNS5.TP <- clean.GLASS.CNS5.TP[clean.GLASS.CNS5.TP$Age_years >= 18, ]
# NOTE: only 1 patient was under 18 years of age which was excluded from the dataset.

# First, create the new column and initialize it with NA
clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis <- NA

# Now, apply the reannotation criteria

# 1) Astrocytoma, IDH-mutant if Mutant in IDH1/2seq and non-codel
table(clean.GLASS.CNS5.TP$IDH1.R132seq.status) # 153 Mutant samples based on IDH1.R132seq
table(clean.GLASS.CNS5.TP$IDH2.R172seq.status) # 0 Mutant samples based on IDH1.R172seq

clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis[
  is.na(clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis) &
    (clean.GLASS.CNS5.TP$IDH1.R132seq.status == "Mutant" | clean.GLASS.CNS5.TP$IDH2.R172seq.status == "Mutant") &
    clean.GLASS.CNS5.TP$Codel == "non-codel"
] <- "Astrocytoma, IDH-mutant"

# 2) Oligodendroglioma, IDH-mutant, 1p/19q-codeleted if Mutant in IDH1/2seq and non-codel
clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis[
  is.na(clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis) &
    (clean.GLASS.CNS5.TP$IDH1.R132seq.status == "Mutant" | clean.GLASS.CNS5.TP$IDH2.R172seq.status == "Mutant") &
    clean.GLASS.CNS5.TP$Codel == "codel"
] <- "Oligodendroglioma, IDH-mutant, 1p/19q-codeleted"

table(clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis, useNA="ifany")

# 3) Glioblastoma, IDH-wildtype
clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis[
  is.na(clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis) &
    (clean.GLASS.CNS5.TP$IDH1.R132seq.status == "Wildtype" |
       clean.GLASS.CNS5.TP$IDH2.R172seq.status == "Wildtype") &
    (
      (clean.GLASS.CNS5.TP$EGFR >= 1 & clean.GLASS.CNS5.TP$EGFR != "unknown") |
        (clean.GLASS.CNS5.TP$TERTp.C228T.C250T.mutation == "Mutant" & clean.GLASS.CNS5.TP$TERTp.C228T.C250T.mutation != "unknown") |
        (clean.GLASS.CNS5.TP$Chr.7.gain.Chr.10.loss == "Gain chr 7 & loss chr 10" & clean.GLASS.CNS5.TP$Chr.7.gain.Chr.10.loss != "unknown")
    ) # condition 4: and has at least one of these categories, explicitly excluding "unknown"
] <- "Glioblastoma, IDH-wildtype"

# 4) Assign unclassified to WHO_CNS5_diagnosis that didn't meet the above criteria

clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis[is.na(clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis)] <- 'Unclassified'

table(clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis, useNA="ifany")

#### Second pass reannotation after manually checking unclassified rows. ################################################

# 1) First capture all remaining IDH1.R132seq Mutant samples (n=21). Classify them as Astrocytoma, even though codel status is not known, IDH-mutant since original histology (i.e., Glioblastoma, Astrocytoma) does not show any oligodendroglial features.
clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis[
  clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis == "Unclassified" &
    (clean.GLASS.CNS5.TP$TCGA.Histology == "Glioblastoma" | clean.GLASS.CNS5.TP$TCGA.Histology == "Astrocytoma") & # Histology
    (clean.GLASS.CNS5.TP$IDH.status == "IDHmut" | clean.GLASS.CNS5.TP$IDH1.R132seq.status == "Mutant" | clean.GLASS.CNS5.TP$IDH2.R172seq.status == "Mutant") # molecular
] <- "Astrocytoma, IDH-mutant"

# 2) Oligodendroglioma, IDH-mutant, 1p/19q-codeleted
clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis[
  clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis == "Unclassified" &
    (clean.GLASS.CNS5.TP$IDH.status == "IDHmut" | clean.GLASS.CNS5.TP$IDH1.R132seq.status == "Mutant" | clean.GLASS.CNS5.TP$IDH2.R172seq.status == "Mutant") &
    clean.GLASS.CNS5.TP$Codel == "codel"
] <- "Oligodendroglioma, IDH-mutant, 1p/19q-codeleted"

# 3) Glioblastoma, IDH-wildtype
clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis[
  clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis == "Unclassified" & # condition 1: take all unclassified
    clean.GLASS.CNS5.TP$IDH.status == "IDHwt" & # condition 2: that are IDHwt
    (
      (clean.GLASS.CNS5.TP$EGFR >= 1 & clean.GLASS.CNS5.TP$EGFR != "unknown") |
        (clean.GLASS.CNS5.TP$TERTp.C228T.C250T.mutation == "Mutant" & clean.GLASS.CNS5.TP$TERTp.C228T.C250T.mutation != "unknown") |
        (clean.GLASS.CNS5.TP$Chr.7.gain.Chr.10.loss == "Gain chr 7 & loss chr 10" & clean.GLASS.CNS5.TP$Chr.7.gain.Chr.10.loss != "unknown")
    ) # condition 4: and has at least one of these categories, explicitly excluding "unknown"
] <- "Glioblastoma, IDH-wildtype"

#### 3rd pass of reannotation for unclassified rows. ######################################################################

# 1) Astrocytoma, IDH-mutant
clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis[
  clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis == "Unclassified" & # condition 1: take all unclassified
    clean.GLASS.CNS5.TP$IDH.status == "IDHmut" & # condition 2: that are IDHmut
    clean.GLASS.CNS5.TP$Codel == "non-codel" # condition 4: and is non-codeleted
] <- "Astrocytoma, IDH-mutant"

# 2) Oligodendroglioma, IDH-mutant, 1p/19q-codeleted
# View a few rows to check if they make sense
glass.CNS5.unclass <- clean.GLASS.CNS5.TP[clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis == "Unclassified", c("sample_barcode", "who_classification", "TCGA.Histology", "TCGA.Grade", "WHO_CNS5_diagnosis", "IDH.status", "IDH1.R132seq.status", "IDH2.R172seq.status", "Codel",  "EGFR", "TERTp.C228T.C250T.mutation", "Chr.7.gain.Chr.10.loss")]


# check if there are anymore codeleted IDHmutant samples in unclassified rows
table(glass.CNS5.unclass$IDH.status)
# IDHwt unknown 
# 42      43 
# NOTE: No need for third pass reannotation for oligodendroglioma as there are no more IDHmut unclassified samples

# 3) Glioblastoma, IDH-wildtype
clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis[
  clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis == "Unclassified" & # condition 1
    clean.GLASS.CNS5.TP$IDH.status == "IDHwt" & # condition 2
    (clean.GLASS.CNS5.TP$TCGA.Histology == "Glioblastoma" |
       clean.GLASS.CNS5.TP$TCGA.Histology == "Astrocytoma")
] <- "Glioblastoma, IDH-wildtype"

table(clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis)

# Astrocytoma, IDH-mutant                 133
# Glioblastoma, IDH-wildtype              124
# Oligodendroglioma, IDH-mutant, 1p/19q-codeleted         32                           
# Unclassified                            41

# Apply the classification
clean.GLASS.CNS5.reannotated <- classify_who_cns5(clean.GLASS.CNS5.TP)

# View final classification distribution
print(table(clean.GLASS.CNS5.reannotated$WHO_CNS5_diagnosis))

#### 4.3. Assign WHO_CNS5_grade (new column) ####
## Grade 4 - if:
# criterion 1: WHO_CNS5_diagnosis == Glioblastoma, IDH-wildtype.
# criterion 2: WHO_CNS5_diagnosis == Astrocytoma, IDH-mutant and CDKN2A or CDKN2B == -2.
## Grade 3 - if:
# criterion 1: WHO_CNS5_diagnosis == Astrocytoma, IDH-mutant and who_classification == grepl("Anaplastic")
# criterion 2: WHO_CNS5_diagnosis == Oligodendroglioma, IDH-mutant, 1p/19q-codeleted. and one of these other criteria has been satisfied: CDKN2A == -2, or who_classification == grepl("Anaplastic")
## Grade 2 - if:
# criterion 1: WHO_CNS5_diagnosis != Glioblastoma, IDH-wildtype.

# First, create the new column and initialize it with NA
clean.GLASS.CNS5.TP$WHO_CNS5_grade <- NA

# 1) Grade 4 - Glioblastoma
clean.GLASS.CNS5.TP$WHO_CNS5_grade[
  is.na(clean.GLASS.CNS5.TP$WHO_CNS5_grade) &
    clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis == "Glioblastoma, IDH-wildtype"
] <- "Grade 4"

# 2) Grade 4 - Astrocytoma, IDH-mutant
clean.GLASS.CNS5.TP$WHO_CNS5_grade[
  is.na(clean.GLASS.CNS5.TP$WHO_CNS5_grade) &
    clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis == "Astrocytoma, IDH-mutant" &
    (clean.GLASS.CNS5.TP$CDKN2A == -2 |
       clean.GLASS.CNS5.TP$CDKN2B == -2)
] <- "Grade 4"

# 3) Grade 3 - Astrocytoma, IDH-mutant
clean.GLASS.CNS5.TP$WHO_CNS5_grade[
  is.na(clean.GLASS.CNS5.TP$WHO_CNS5_grade) &
    clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis == "Astrocytoma, IDH-mutant" &
    (grepl("Anaplastic", clean.GLASS.CNS5.TP$who_classification)
    )
] <- "Grade 3"

# 4) Grade 3 -  Oligodendroglioma, IDH-mutant, 1p/19q-codeleted
clean.GLASS.CNS5.TP$WHO_CNS5_grade[
  is.na(clean.GLASS.CNS5.TP$WHO_CNS5_grade) &
    clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis == "Oligodendroglioma, IDH-mutant, 1p/19q-codeleted" &
    (grepl("Anaplastic", clean.GLASS.CNS5.TP$who_classification) |
       clean.GLASS.CNS5.TP$CDKN2A == -2)
] <- "Grade 3"

table(clean.GLASS.CNS5.TP$WHO_CNS5_grade, useNA="ifany")


# 5) If CDKN2A/B loss cannot be used for grading and WHO_CNS5_grade is NA. copy grading from original annotation.

# Grade IV to 4 
clean.GLASS.CNS5.TP$WHO_CNS5_grade[
  is.na(clean.GLASS.CNS5.TP$WHO_CNS5_grade) &
    clean.GLASS.CNS5.TP$TCGA.Grade == "IV" 
] <- "Grade 4"

# Grade III to 3 
clean.GLASS.CNS5.TP$WHO_CNS5_grade[
  is.na(clean.GLASS.CNS5.TP$WHO_CNS5_grade) &
    clean.GLASS.CNS5.TP$TCGA.Grade == "III" 
] <- "Grade 3"

# Grade II to 2 
clean.GLASS.CNS5.TP$WHO_CNS5_grade[
  is.na(clean.GLASS.CNS5.TP$WHO_CNS5_grade) &
    clean.GLASS.CNS5.TP$TCGA.Grade == "II" 
] <- "Grade 2"

clean.GLASS.CNS5.TP$WHO_CNS5_grade[is.na(clean.GLASS.CNS5.TP$WHO_CNS5_grade)] <- 'Unclassified'

table(clean.GLASS.CNS5.TP$WHO_CNS5_grade, useNA = "ifany")
# Grade 2      Grade 3      Grade 4 Unclassified 
#      86           16          207           21 

#### 4.4. Assign WHO_CNS5_IDH.status (new column)####
# First, let's check if the column exists, and if not, create it
if(!"WHO_CNS5_IDH.status" %in% names(clean.GLASS.CNS5.TP)) {
  clean.GLASS.CNS5.TP$WHO_CNS5_IDH.status <- NA
}

#### IDH mutant based on IDH1.R132seq or IDH2.R172seq data
clean.GLASS.CNS5.TP$WHO_CNS5_IDH.status[
  is.na(clean.GLASS.CNS5.TP$WHO_CNS5_IDH.status) &
    (clean.GLASS.CNS5.TP$IDH1.R132seq.status == "Mutant" | clean.GLASS.CNS5.TP$IDH2.R172seq.status == "Mutant")
] <- "IDH-mutant"

#### IDH wildtype based on IDH1.R132seq or IDH2.R172seq data
clean.GLASS.CNS5.TP$WHO_CNS5_IDH.status[
  is.na(clean.GLASS.CNS5.TP$WHO_CNS5_IDH.status) &
    (clean.GLASS.CNS5.TP$IDH1.R132seq.status == "Wildtype" | clean.GLASS.CNS5.TP$IDH2.R172seq.status == "Wildtype")
] <- "IDH-wildtype"

### Assign unclassified to WHO_CNS5_IDH.status that didn't meet the above criteria
clean.GLASS.CNS5.TP$WHO_CNS5_IDH.status[is.na(clean.GLASS.CNS5.TP$WHO_CNS5_IDH.status)] <- 'Unclassified'

table(clean.GLASS.CNS5.TP$WHO_CNS5_IDH.status, useNA="ifany")
# IDH-mutant IDH-wildtype Unclassified 
# 153           16          161

#### IDH mutant based on original status if seq data says otherwise
clean.GLASS.CNS5.TP$WHO_CNS5_IDH.status[
  clean.GLASS.CNS5.TP$WHO_CNS5_IDH.status == "Unclassified" &
    clean.GLASS.CNS5.TP$IDH.status == "IDHmut"
] <- "IDH-mutant"

#### IDH wildtype based on on original status if seq data says otherwise
clean.GLASS.CNS5.TP$WHO_CNS5_IDH.status[
  clean.GLASS.CNS5.TP$WHO_CNS5_IDH.status == "Unclassified" &
    clean.GLASS.CNS5.TP$IDH.status == "IDHwt"
] <- "IDH-wildtype"

table(clean.GLASS.CNS5.TP$WHO_CNS5_IDH.status, useNA="ifany")
# IDH-mutant IDH-wildtype         <NA> 
# 171        125                  34 

#### Create WHO_CNS5_classification column ####
clean.GLASS.CNS5.TP$WHO_CNS5_classification <- paste0(clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis, ". ", clean.GLASS.CNS5.TP$WHO_CNS5_grade, ".")

table(clean.GLASS.CNS5.TP$WHO_CNS5_classification, useNA = "ifany")

#### Create WHO_CNS5_histology column ####
clean.GLASS.CNS5.TP$WHO_CNS5_histology <- 
  ifelse(clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis == 'Glioblastoma, IDH-wildtype', 'Glioblastoma', 
         ifelse(clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis == 'Astrocytoma, IDH-mutant', 'Astrocytoma',
                ifelse(clean.GLASS.CNS5.TP$WHO_CNS5_diagnosis == 'Oligodendroglioma, IDH-mutant, 1p/19q-codeleted', 'Oligodendroglioma', 
                       'Unclassified')))

table(clean.GLASS.CNS5.TP$WHO_CNS5_histology)

save(clean.GLASS.CNS5.TP, file = "OP/clean_GLASS.TP.WHO-CNS5_v3.RData")
write.csv(clean.GLASS.CNS5.TP,"OP/clean_GLASS.TP.WHO-CNS5_v3.csv", row.names = FALSE)

#### create summary for notes. what is the frequency of the samples that were changed. ####

# TCGA summary
ijms <- WHO.CNS5[ , c("Patient_ID", "WHO_CNS5_diagnosis.original.ijms.2057006")]
tcga.summary <- clean.WHO.CNS5_v2[ , c("Patient_ID", "WHO_CNS5_histology", "WHO_CNS5_grade", "WHO_CNS5_IDHstatus", "WHO_CNS_1p.19q._codeletion")]
st(tcga.summary, out='csv', file='OP/full_TCGAijms_pheno_summary_v3.csv')

# GLASS summary
clean.GLASS.CNS5.TP

glass.summary <- clean.GLASS.CNS5.TP[ , c("case_barcode", "Gender", "Age_years", "Vital_status", "Survival_Months", "TCGA.Histology", "TCGA.Grade", "IDH.status", "Codel", "who_classification", "MGMT.promoter.status", "idh_codel_subtype", "WHO_CNS5_IDHstatus", "IDH.codel.subtype", "WHO_CNS5_diagnosis", "WHO_CNS5_grade", "WHO_CNS5_IDH.status", "WHO_CNS5_classification", "WHO_CNS5_histology")]

st(glass.summary, out='csv', file='OP/GLASS_pheno_summary_v3.csv')


