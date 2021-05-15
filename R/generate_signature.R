###############################################################################
##                           Ranked file                                     ##
###############################################################################
library(data.table)

# series_treatment_vs_outcome <- "GSE45757", "GSE112282", "GSE14426", 
# discarded: GSE37645
# series_gene_exp_vs_outcome <- c("GSE28735", "GSE62452", "GSE21501", 
# "GSE62165", "GSE57495")

########################
# Treatment vs outcome #
########################

series1 <- read.csv(file = "GSE112282_common_genes.csv", header = TRUE)
series2 <- read.csv(file = "GSE45757_common_genes.csv", header = TRUE)
# series3 = read.csv(file = "GSE14426_common_genes.csv", header = TRUE)
series3 <- read.csv(file = "GSE14426_common_genes.csv", header = TRUE)
series4 <- read.csv(file = "GSE55643_common_genes.csv", header = TRUE)
series5 <- read.csv(file = "GSE77858_common_genes.csv", header = TRUE)
series6 <- read.csv(file = "GSE11838_common_genes.csv", header = TRUE)

common_treatment_vs_outcome <- rbindlist(mget(paste0("series", 1:6)))[, .N, HGNC]
common_treatment_vs_outcome[order(-N,HGNC)]
colnames(common_treatment_vs_outcome)[2] <- "count"
# Save as CSV
write.csv(common_treatment_vs_outcome[order(-count,HGNC)],
          file="signature_treatment_vs_outcome.csv", 
          quote = FALSE,
          row.names = FALSE)

##############################
# Gene expression vs outcome #
##############################

series20 <- read.csv(file = "GSE21501_common_genes.csv", header = TRUE)
series21 <- read.csv(file = "GSE28735_common_genes.csv", header = TRUE)
series22 <- read.csv(file = "GSE62165_common_genes.csv", header = TRUE)
series23 <- read.csv(file = "GSE71729_common_genes.csv", header = TRUE)
series24 <- read.csv(file = "GSE56560_common_genes.csv", header = TRUE)
# series25 <- read.csv(file = "GSE15471_common_genes.csv", header = TRUE)


common_g_expr_vs_outcome <- rbindlist(mget(paste0("series", 20:24)))[, .N, HGNC]
common_g_expr_vs_outcome[order(-N,HGNC)]
colnames(common_g_expr_vs_outcome)[2] <- "count"
# Save as CSV
write.csv(common_g_expr_vs_outcome[order(-count,HGNC)],
          file="signature_gene_expression_vs_outcome.csv", 
          quote = FALSE,
          row.names = FALSE)

################################################
# Number of genes grouping by number of series #
################################################

table(common_treatment_vs_outcome$count)
table(common_g_expr_vs_outcome$count)

##################
# Common to both #
##################
subset1 <- as.data.frame(common_treatment_vs_outcome[common_treatment_vs_outcome$count >= 2]$HGNC)
subset2 <- as.data.frame(common_g_expr_vs_outcome[common_g_expr_vs_outcome$count==2]$HGNC)
colnames(subset1)[1] <- "HGNC"
colnames(subset2)[1] <- "HGNC"

common <- rbindlist(mget(paste0("subset", 1:2)), use.names = FALSE)[, .N, HGNC]
common[order(-N,HGNC)]
colnames(common)[2] <- "count"
# Save as CSV
common_to_both_groups <- common[common$count == 2]
write.csv(common_to_both_groups[order(-count,HGNC)]$HGNC,
          file="signature_final.csv", 
          quote = FALSE,
          row.names = FALSE)
