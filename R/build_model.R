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
# metadata(assay)

# Get expression data
eset <- study@ExperimentList@listData[["RNA_Seq_v2_mRNA_median_all_sample_Zscores"]]@assays@data@listData[[1]]


metadata(study)


##########################################
#              Correlation               #
##########################################

# maemerge <- mergeReplicates(intersectColumns(study))
# study_wide <- wideFormat(maemerge, colDataCols=c("OS_STATUS","PFS_STATUS"))

subacc <- study[, , "RNA_Seq_v2_mRNA_median_all_sample_Zscores"]
nrow(assay(subacc)) # 20531
ncol(assay(subacc)) # 177
# subacc <- intersectColumns(subacc)
# subacc <- intersectRows(subacc)

subacc_wide <- wideFormat(subacc, colDataCols=c("OS_MONTHS","PFS_MONTHS"), )
#subacc_wide_orig <- wideFormat(subacc, colDataCols=c("OS_MONTHS","PFS_MONTHS"), )

nrow(subacc_wide) # 177
ncol(subacc_wide) # 20534

# rownames(subacc_wide)
col_names <- colnames(subacc_wide)
#length((col_names))

# Remove the prefix with the experiment name in all the columns
new_col_names <- sub("RNA_Seq_v2_mRNA_median_all_sample_Zscores_", "", col_names)
#length((new_col_names))

# Set the new column names
colnames(subacc_wide) <- new_col_names

# subacc_wide[0:3,0:5]
# subacc_wide[0:3,20531:20534]
# subacc_wide_orig[0:3,20531:20534]
# subacc_wide[0:3,-1]

# gene_signature_df <- read.csv(file = "signature_final.csv", header = TRUE)
# gene_signature <- gene_signature_df[['x']]

gene_signature_df <- read.csv(file = "signature_gene_expression_vs_outcome.csv", header = TRUE)
gene_signature <- gene_signature_df[gene_signature_df$count >= 2, ][['HGNC']]

length(gene_signature)
genes_in_common <- intersect(gene_signature,colnames(subacc_wide))
length(genes_in_common)

cols_subset <- c(genes_in_common, c("OS_MONTHS","PFS_MONTHS"))
length(cols_subset)

# [,-1] to remove the first column which contains the sample name
# mycors <- cor(as.matrix(subacc_wide[,-1]))

# subacc_wide[,cols_subset]
mycors <- cor(as.matrix(subacc_wide[,cols_subset]))
# mycors <- abs(mycors)
diag(mycors) <- NA

# show only components that have at least one correlation greater than 0.5.
# has.high.cor <- apply(mycors, 2, max, na.rm=TRUE) > 0.5
# mycors <- mycors[has.high.cor, has.high.cor]
# pheatmap::pheatmap(mycors)

mycors[lower.tri(mycors)] <- NA
pheatmap::pheatmap(mycors, cluster_row =F, cluster_cols=F, na_col="white")

library(reshape2)
library(ggplot2)

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

mycors <- cor(as.matrix(subacc_wide[,cols_subset]))
# Reorder the correlation matrix
cormat <- reorder_cormat(mycors)
upper_tri <- get_upper_tri(mycors)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
# ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
#   geom_text(aes(Var2, Var1, label = value), color = "black", size = 1) +
#   geom_tile(color = "white") +
#   scale_fill_gradient2(low = "blue", high = "red", mid = "white",
#                        midpoint = 0, limit = c(-1,1),
#                        space = "Lab", name = "Pearson\ncorrelation") +
#   # theme_minimal()+ 
#   theme(
#     axis.text.x = element_text(angle = 45, vjust = 1, 
#                                size = 8, hjust = 1),
#     axis.text.y = element_text(size = 8),
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.border = element_blank(),
#     panel.background = element_blank(),
#     axis.ticks = element_blank(),
#     legend.justification = c(1, 0),
#     legend.position = c(0.6, 0.7),
#     legend.direction = "horizontal") +
#   guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
#                                title.position = "top", title.hjust = 0.5))
# # Print the heatmap
# print(ggheatmap)


ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = sprintf("%0.2f", round(value, digits = 2))), 
            color = "black", 
            size = 2) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, 
                               size = 8, hjust = 1),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

print(ggheatmap)
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
