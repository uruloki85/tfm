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

subacc_wide <- wideFormat(subacc, colDataCols=c("OS_MONTHS","PFS_MONTHS"), )
nrow(subacc_wide) # 177
ncol(subacc_wide) # 20534

fix_col_names <- function(subacc) {
  col_names <- colnames(subacc)
  # Remove the prefix with the experiment name in all the columns
  new_col_names <- sub("RNA_Seq_v2_mRNA_median_all_sample_Zscores_", "", col_names)
  # Set the new column names
  colnames(subacc) <- new_col_names
  return(subacc)
}

subacc_wide <- fix_col_names(subacc_wide)
# colnames(subacc_wide)

subacc_wide_only_dead <- wideFormat(subacc, colDataCols=c("OS_MONTHS","PFS_MONTHS", "OS_STATUS"), )
subacc_wide_only_dead <- subacc_wide_only_dead[subacc_wide_only_dead$OS_STATUS == "1:DECEASED",]
# Drop column
subacc_wide_only_dead$OS_STATUS <- NULL
nrow(subacc_wide_only_dead) # 92
ncol(subacc_wide_only_dead) # 20534

subacc_wide_only_dead <- fix_col_names(subacc_wide_only_dead)
# colnames(subacc_wide_only_dead)

# subacc_wide[0:3,0:5]
# subacc_wide[0:3,20531:20534]
# subacc_wide_orig[0:3,20531:20534]
# subacc_wide[0:3,-1]

# gene_signature_df <- read.csv(file = "signature_final.csv", header = TRUE)
# gene_signature <- gene_signature_df[['x']]

# Get signature for "Gene expression vs outcome"
signature_gene_exp_vs_outcome_df <- read.csv(file = "signature_gene_expression_vs_outcome.csv", header = TRUE)
signature_gene_exp_vs_outcome <- signature_gene_exp_vs_outcome_df[signature_gene_exp_vs_outcome_df$count >= 3, ][['HGNC']]
length(signature_gene_exp_vs_outcome) #

# Get signature for "Treatment vs outcome"
signature_treatment_vs_outcome_df <- read.csv(file = "signature_treatment_vs_outcome.csv", header = TRUE)
signature_treatment_vs_outcome <- signature_treatment_vs_outcome_df[signature_treatment_vs_outcome_df$count >= 3, ][['HGNC']]
length(signature_treatment_vs_outcome) #

# Get signature common to both groups
signature_final_df <- read.csv(file = "signature_final.csv", header = TRUE)
signature_final <- signature_final_df[signature_final_df$count >= 2, ][['HGNC']]
length(signature_final) #


library(reshape2)
library(ggplot2)

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
  return (cormat)
}

create_ggheatmap <- function(metled_corr, title) {
  ggheatmap <- ggplot(metled_corr, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name=paste(title,"Pearson","Correlation",sep="\n")) +
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
  return (ggheatmap)
}

calculate_corr_and_heatmap <- function(subacc, title, gene_signature) {
  # Find genes in common
  genes_in_common <- intersect(gene_signature,colnames(subacc))
  # length(genes_in_common)
  cols_subset <- c(genes_in_common, c("OS_MONTHS","PFS_MONTHS"))
  # length(cols_subset)
  
  mycors <- cor(as.matrix(subacc[,cols_subset]))
  # Reorder the correlation matrix
  cormat <- reorder_cormat(mycors)
  # cormat <- mycors
  cormat[lower.tri(cormat)] <- NA
  # Melt the correlation matrix
  melted_cormat <- melt(cormat, na.rm = TRUE)
  ggheatmap <- create_ggheatmap(melted_cormat,title)
  return (ggheatmap)
}


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
