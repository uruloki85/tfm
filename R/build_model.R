library(cBioPortalData)
# library(AnVIL)

cbio <- cBioPortal()

study <- cBioDataPack(
  "paad_tcga_pan_can_atlas_2018",
  use_cache = TRUE,
  names.field = c("Hugo_Symbol", "Entrez_Gene_Id", "Gene"),
  ask = TRUE
)

##################
# Get signatures #
##################

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

####################
# Prepare the data #
####################

subacc <- study[, , "RNA_Seq_v2_mRNA_median_all_sample_Zscores"]
nrow(assay(subacc)) # 20531
ncol(assay(subacc)) # 177
# subacc <- intersectColumns(subacc)
# subacc <- intersectRows(subacc)
my_data <- wideFormat(subacc, colDataCols=c("PFS_MONTHS"), )
# Drop column
my_data$primary <- NULL

my_data <- fix_col_names(my_data)

# Examine variable
min(my_data$PFS_MONTHS)
max(my_data$PFS_MONTHS)

# Histogram
hist(my_data$PFS_MONTHS, main="Histogram", xlab="PFS_MONTHS", col="light green")

# Median
median(my_data$PFS_MONTHS)

ncol(my_data) # 20532
nrow(my_data) # 177

find_cols_subset <- function(gene_signature, my_df) {
  # Find genes in common
  genes_in_common <- intersect(gene_signature,colnames(my_df))
  # length(genes_in_common)
  cols_subset <- c(genes_in_common, c("PFS_MONTHS"))
  # length(cols_subset)
  
  return(cols_subset)
}

# Subset the data to only those columns of interest
# my_data <- my_data[,cols_subset]

# Subset data to columns of interest
my_data_treatment_vs_outcome <- my_data[,find_cols_subset(signature_treatment_vs_outcome, my_data)]
ncol(my_data_treatment_vs_outcome) # 10
nrow(my_data_treatment_vs_outcome) # 177

my_data_gene_exp_vs_outcome <- my_data[,find_cols_subset(signature_gene_exp_vs_outcome, my_data)]
ncol(my_data_gene_exp_vs_outcome) # 18
nrow(my_data_gene_exp_vs_outcome) # 177

my_data_signature_final <- my_data[,find_cols_subset(signature_final, my_data)]
ncol(my_data_signature_final) # 18
nrow(my_data_signature_final) # 177


##################
# Missing values #
##################

# Explore missing values
sapply(my_data_treatment_vs_outcome, function(x) sum(is.na(x)))
sapply(my_data_gene_exp_vs_outcome, function(x) sum(is.na(x)))
sapply(my_data_signature_final, function(x) sum(is.na(x)))

######################
# Correlation matrix #
######################

#1st step
# Correlation matrix with all genes
my_cors2 <- cor(as.matrix(my_data))
my_cors2_filtered <- my_cors2["PFS_MONTHS",]
#abs(my_cors2_filtered) <= .90 & 
index <- which(abs(my_cors2_filtered) >= .15, arr.ind = T)
length(my_cors2_filtered[index])
# my_cors2_filtered[index]
filtered_gene_names <- names(my_cors2_filtered[index])
# filtered_gene_names <- filtered_gene_names[-length(filtered_gene_names)]

my_data <- my_data[,filtered_gene_names]

# 2nd step
# Correlation matrix with the genes selected in the previous step
my_cors3 <- cor(as.matrix(my_data[,filtered_gene_names]))
# Reorder the correlation matrix
cormat3 <- reorder_cormat(my_cors3)
cormat3[lower.tri(cormat3)] <- NA
# Melt the correlation matrix
melted_cormat3 <- melt(cormat3, na.rm = TRUE)
ggheatmap <- create_ggheatmap(melted_cormat3,"")
print(ggheatmap)

# Discard the diagonal of the matrix
diag(my_cors3) <- NA 
# Check if there is any value higher than 0.9
index <- which(abs(my_cors3) > .9, arr.ind = T) 
if (index > 0) {
  cbind.data.frame(gene1 = rownames(my_cors3)[index[,1]], # get the row name 
                   gene2 = colnames(my_cors3)[index[,2]]) # get the column name
  names(my_data[index])
}

###################################
#               Model             #
###################################
library(randomForest)

# my_data_saved <- my_data
# my_data <- my_data_saved

# Prepare the class variable using the median of PFS_MONTHS
median_pfs <- median(my_data$PFS_MONTHS)
my_data$class[my_data$PFS_MONTHS <= median_pfs] <- "0"
my_data$class[my_data$PFS_MONTHS > median_pfs] <- "1"
my_data$class <- as.factor(my_data$class)
# Check new column
my_data[,c("PFS_MONTHS","class")]
# Remove PFS_MONTHS
my_data$PFS_MONTHS <- NULL

my_data$class

my_matrix <- as.matrix(my_data)
my_df <- as.data.frame(my_matrix)
my_df$class <- as.factor(my_df$class)

# Set random seed to make results reproducible:
set.seed(42)

# Calculate the size of each of the data sets:
data_set_size <- floor(nrow(my_df) * .7)
# Generate a random sample of "data_set_size" indexes
indexes <- sample(1:nrow(my_df), size = data_set_size)

# Assign the data to the correct sets
training <- my_df[indexes,]
validation1 <- my_df[-indexes,]
nrow(training)
nrow(validation1)

rf_classifier <- randomForest(class ~ ., data=training, ntree=100, mtry=2)
rf_classifier
# varImpPlot(rf_classifier)

#predicting the class for the test data set
tree.pred1 <- predict(rf_classifier,validation1,type="class")

#creating confusion matrix
table(tree.pred1,validation1$class) 



