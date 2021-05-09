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

# Missing values
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

# Experiment data


experiments(study)

experiments(study[["RNA_Seq_v2_mRNA_median_all_sample_Zscores"]])
assay <- assay(study[["RNA_Seq_v2_mRNA_median_all_sample_Zscores"]])
assay
metadata(assay)

# Get expression data
eset <- study@ExperimentList@listData[["RNA_Seq_v2_mRNA_median_all_sample_Zscores"]]@assays@data@listData[[1]]


metadata(study)


###############
# Correlation #
###############

# maemerge <- mergeReplicates(intersectColumns(study))
# study_wide <- wideFormat(maemerge, colDataCols=c("OS_STATUS","PFS_STATUS"))

subacc <- study[, , "RNA_Seq_v2_mRNA_median_all_sample_Zscores"]
# subacc <- intersectColumns(subacc)
# subacc <- intersectRows(subacc)

subacc_wide <- wideFormat(subacc, colDataCols=c("OS_MONTHS","PFS_MONTHS"))
nrow(subacc_wide) # 177
ncol(subacc_wide) # 20534
# rownames(subacc_wide)
# colnames(subacc_wide)
# 
# subacc_wide[0:3,0:5]
# subacc_wide[0:3,-1]

# [,-1] to remove the first column which contains the sample name
mycors <- cor(as.matrix(subacc_wide[,-1]))
mycors <- abs(mycors)
# diag(mycors) <- NA

has.high.cor <- apply(mycors, 2, max, na.rm=TRUE) > 0.5
mycors <- mycors[has.high.cor, has.high.cor]
pheatmap::pheatmap(mycors)





# # Get the assays
# subacc.list <- assays(subacc)
# # Transpose rows to columns, so genes are in columns
# subacc.list <- lapply(subacc.list, t)
# 
# colnames(subacc.list[["RNA_Seq_v2_mRNA_median_all_sample_Zscores"]])
# rownames(subacc.list[["RNA_Seq_v2_mRNA_median_all_sample_Zscores"]])
# 
# corres <- cor(subacc.list[[1]])
# has.high.cor <- apply(corres, 2, max, na.rm=TRUE) > 0.5
# mycors <- corres[has.high.cor, has.high.cor]
# 
# library(pheatmap)
# pheatmap::pheatmap(mycors)



##################################
# Kaplan-Meier Survival Analysis #
##################################
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
