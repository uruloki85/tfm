library(cBioPortalData)
# library(AnVIL)

cbio <- cBioPortal()

study <- cBioDataPack(
  "paad_tcga_pan_can_atlas_2018",
  use_cache = TRUE,
  names.field = c("Hugo_Symbol", "Entrez_Gene_Id", "Gene"),
  ask = TRUE
)

study

# Metadata

colData(study)
names(colData(study))

table(study$OS_STATUS)
# table(study$OS_MONTHS)
table(study$DSS_STATUS)
# table(study$DSS_MONTHS)
table(study$DFS_STATUS)
# table(study$DFS_MONTHS)
table(study$PFS_STATUS)
# table(study$PFS_MONTHS)
table(study$AGE)
table(study$SEX)

# How many samples are present in each of the experiments
library(UpSetR)
upsetSamples(study)

# Explore missing values
missing_values <- data.frame(sum(is.na(study$OS_STATUS)))
colnames(missing_values) <- "OS_STATUS"
rownames(missing_values) <- "missing"
missing_values$OS_MONTHS <- sum(is.na(study$OS_MONTHS))
missing_values$DSS_STATUS <- sum(is.na(study$DSS_STATUS))
missing_values$DSS_MONTHS <- sum(is.na(study$DSS_MONTHS))
missing_values$DFS_STATUS <- sum(is.na(study$DFS_STATUS))
missing_values$DFS_MONTHS <- sum(is.na(study$DFS_MONTHS))
missing_values$PFS_STATUS <- sum(is.na(study$PFS_STATUS))
missing_values$PFS_MONTHS <- sum(is.na(study$PFS_MONTHS))
missing_values

####################################################
#              Correlation Analysis                #
####################################################

subacc <- study[, , "RNA_Seq_v2_mRNA_median_all_sample_Zscores"]
nrow(assay(subacc)) # 20531
ncol(assay(subacc)) # 177
# subacc <- intersectColumns(subacc)
# subacc <- intersectRows(subacc)

subacc_wide <-
  wideFormat(subacc,
             colDataCols = c("OS_MONTHS", "DSS_MONTHS", "PFS_MONTHS"),
  )
nrow(subacc_wide) # 177
ncol(subacc_wide) # 20535

subacc_wide <- fix_col_names(subacc_wide)
# colnames(subacc_wide)

subacc_wide_only_dead <-
  wideFormat(subacc,
             colDataCols = c("OS_MONTHS", "DSS_MONTHS", "PFS_MONTHS", "OS_STATUS"),
  )
subacc_wide_only_dead <-
  subacc_wide_only_dead[subacc_wide_only_dead$OS_STATUS == "1:DECEASED", ]
# Drop column
subacc_wide_only_dead$OS_STATUS <- NULL
nrow(subacc_wide_only_dead) # 92
ncol(subacc_wide_only_dead) # 20535

subacc_wide_only_dead <- fix_col_names(subacc_wide_only_dead)
# colnames(subacc_wide_only_dead)

# subacc_wide[0:3,0:5]
# subacc_wide[0:3,20531:20534]
# subacc_wide_orig[0:3,20531:20534]
# subacc_wide[0:3,-1]

# gene_signature_df <- read.csv(file = "signature_final.csv", header = TRUE)
# gene_signature <- gene_signature_df[['x']]

# Get signature for "Gene expression vs outcome"
signature_gene_exp_vs_outcome_df <-
  read.csv(file = "signature_gene_expression_vs_outcome.csv", header = TRUE)
signature_gene_exp_vs_outcome <-
  signature_gene_exp_vs_outcome_df[signature_gene_exp_vs_outcome_df$count >= 3,][['HGNC']]
length(signature_gene_exp_vs_outcome) #

# Get signature for "Treatment vs outcome"
signature_treatment_vs_outcome_df <-
  read.csv(file = "signature_treatment_vs_outcome.csv", header = TRUE)
signature_treatment_vs_outcome <-
  signature_treatment_vs_outcome_df[signature_treatment_vs_outcome_df$count >= 3,][['HGNC']]
length(signature_treatment_vs_outcome) #

# Get signature common to both groups
signature_final_df <- read.csv(file = "signature_final.csv", header = TRUE)
signature_final <- signature_final_df[['x']]
length(signature_final) #

###################################
# Calculate and draw correlations #
###################################

# Gene expression vs outcome
ggheatmap <- calculate_corr_and_heatmap(subacc_wide, 
                                        "All patients",
                                        signature_gene_exp_vs_outcome)
print(ggheatmap)
ggheatmap_dead <- calculate_corr_and_heatmap(subacc_wide_only_dead, 
                                             "Only dead patients",
                                             signature_gene_exp_vs_outcome)
print(ggheatmap_dead)

# Treatment vs outcome
ggheatmap <- calculate_corr_and_heatmap(subacc_wide, 
                                        "All patients",
                                        signature_treatment_vs_outcome)
print(ggheatmap)
ggheatmap_dead <- calculate_corr_and_heatmap(subacc_wide_only_dead, 
                                             "Only dead patients",
                                             signature_treatment_vs_outcome)
print(ggheatmap_dead)

# Final signature
ggheatmap <- calculate_corr_and_heatmap(subacc_wide, 
                                        "All patients",
                                        signature_final)
print(ggheatmap)
ggheatmap_dead <- calculate_corr_and_heatmap(subacc_wide_only_dead, 
                                             "Only dead patients",
                                             signature_final)
print(ggheatmap_dead)



####################################################
#           Kaplan-Meier Survival Analysis         #
####################################################
library(survival)
library(survminer)

# Surv(study$OS_MONTHS, as.numeric(as.factor(study$OS_STATUS)))

# Remove patients missing survival information
study_surv <- study[, complete.cases(study$OS_MONTHS, study$OS_STATUS), ]
study_surv_df <- as.data.frame(colData(study_surv))

# Convert OS_STATUS to numeric
study_surv_df$OS_STATUS <- as.numeric(as.factor(study$OS_STATUS))

# Survival analysis
fit_surv <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ PATH_N_STAGE, data = study_surv_df)
ggsurvplot(fit_surv, data = colData(study_surv), risk.table = TRUE)
