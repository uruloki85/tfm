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

find_genes_threshold <- function(my_df, threshold) {
  # Correlation matrix with all genes
  my_cors2 <- cor(as.matrix(my_df))
  my_cors2_filtered <- my_cors2["PFS_MONTHS",]
  #abs(my_cors2_filtered) <= .90 & 
  index <- which(abs(my_cors2_filtered) >= threshold, arr.ind = T)
  length(my_cors2_filtered[index])
  # my_cors2_filtered[index]
  filtered_gene_names <- names(my_cors2_filtered[index])
  # filtered_gene_names <- filtered_gene_names[-length(filtered_gene_names)]
  
  return(filtered_gene_names)
}

# Treatment vs outcome
#1st step
filtered_genes_treatment <- find_genes_threshold(my_data_treatment_vs_outcome, 
                                                            0.1)
filtered_genes_treatment

my_data_treatment_vs_outcome <-
  my_data_treatment_vs_outcome[, filtered_genes_treatment]
ncol(my_data_treatment_vs_outcome) # 5
nrow(my_data_treatment_vs_outcome) # 177

# 2nd step
ggheatmap <- calculate_corr_and_heatmap(my_data_treatment_vs_outcome, 
                                        "Treatment vs outcome",
                                        NULL,
                                        TRUE)
print(ggheatmap)

# Gene expression vs outcome
#1st step
filtered_genes_gene_exp <- find_genes_threshold(my_data_gene_exp_vs_outcome, 
                                                           0.1)
filtered_genes_gene_exp

my_data_gene_exp_vs_outcome <-
  my_data_gene_exp_vs_outcome[, filtered_genes_gene_exp]
ncol(my_data_gene_exp_vs_outcome) # 13
nrow(my_data_gene_exp_vs_outcome) # 177

# 2nd step
ggheatmap <- calculate_corr_and_heatmap(my_data_gene_exp_vs_outcome, 
                                        "Gene expression vs outcome",
                                        NULL,
                                        TRUE)
print(ggheatmap)

# Common signature
#1st step
filtered_genes_final_signature <- find_genes_threshold(my_data_signature_final, 
                                                       0.1)
filtered_genes_final_signature

my_data_signature_final <-
  my_data_signature_final[, filtered_genes_final_signature]
ncol(my_data_signature_final) # 12
nrow(my_data_signature_final) # 177

# 2nd step
ggheatmap <- calculate_corr_and_heatmap(my_data_signature_final, 
                                        "Common signature",
                                        NULL,
                                        TRUE)
print(ggheatmap)


############################################################
## Only useful if the is some correlation higher than 90% ##
############################################################

# Discard the diagonal of the matrix
diag(my_cors3) <- NA 
# Check if there is any value higher than 0.9
index <- which(abs(my_cors3) > .9, arr.ind = T) 
if (index > 0) {
  cbind.data.frame(gene1 = rownames(my_cors3)[index[,1]], # get the row name 
                   gene2 = colnames(my_cors3)[index[,2]]) # get the column name
  names(my_data[index])
}

##########################################
#               Prepare data             #
##########################################

# my_data_saved <- my_data
# my_data <- my_data_saved

init_class_col <- function(my_df) {
  # my_df <- my_data_treatment_vs_outcome
  
  # Prepare the class variable using the median of PFS_MONTHS
  median_pfs <- median(my_df$PFS_MONTHS)
  my_df$class[my_df$PFS_MONTHS <= median_pfs] <- "X0"
  my_df$class[my_df$PFS_MONTHS > median_pfs] <- "X1"
  # my_df$class <- as.factor(my_df$class)
  # Check new column
  # my_df[,c("PFS_MONTHS","class")]
  # Remove PFS_MONTHS
  my_df$PFS_MONTHS <- NULL
  
  # my_data$class
  
  my_df2 <- as.data.frame(my_df@listData)
  
  # my_matrix <- data.matrix(my_df) # as.matrix(my_df)
  # my_df2 <- as.data.frame(my_matrix)
  my_df2$class <- as.factor(my_df2$class)
  # my_df2
  
  return(my_df2)
}

# Initialize variable to predict: 0 (Bad prognosis), 1 (Good prognosis)
my_data_treatment_vs_outcome <- init_class_col(my_data_treatment_vs_outcome)
my_data_treatment_vs_outcome$class

my_data_gene_exp_vs_outcome <- init_class_col(my_data_gene_exp_vs_outcome)
my_data_gene_exp_vs_outcome$class

my_data_signature_final <- init_class_col(my_data_signature_final)
my_data_signature_final$class

# Set random seed to make results reproducible:
set.seed(42)

library(caTools)

# Split data into training and test
sample_treatment <- sample.split(my_data_treatment_vs_outcome$class, 
                                 SplitRatio = .6)
training_treatment <- subset(my_data_treatment_vs_outcome, 
                             sample_treatment == TRUE)
test_treatment <- subset(my_data_treatment_vs_outcome, 
                         sample_treatment == FALSE)
nrow(training_treatment)
nrow(test_treatment)
# Check distribution of values -> should be similar in both sets
table(training_treatment$class)
table(test_treatment$class)

# Split data into training and test
sample_gene_exp <- sample.split(my_data_gene_exp_vs_outcome$class, 
                                SplitRatio = .6)
training_gene_exp <- subset(my_data_gene_exp_vs_outcome, 
                            sample_gene_exp == TRUE)
test_gene_exp <- subset(my_data_gene_exp_vs_outcome, 
                        sample_gene_exp == FALSE)
nrow(training_gene_exp)
nrow(test_gene_exp)
table(training_gene_exp$class)
table(test_gene_exp$class)

# Split data into training and test
sample_final <- sample.split(my_data_signature_final$class, 
                             SplitRatio = .7)
training_common <- subset(my_data_signature_final, 
                          sample_final == TRUE)
test_common <- subset(my_data_signature_final, 
                      sample_final == FALSE)
nrow(training_common)
nrow(test_common)
table(training_common$class)
table(test_common$class)

#################################
#              Model            #
#################################

########################
## With caret package ##
########################

library(caret)
library(e1071)

set.seed(42)

length((filtered_genes_treatment)) # 5
# filtered_genes_gene_exp
train <- training_treatment
test <- test_treatment

# length((filtered_genes_gene_exp)) # 13
# # filtered_genes_gene_exp
# train <- training_gene_exp
# test <- test_gene_exp

# length((filtered_genes_final_signature)) # 12
# #filtered_genes_final_signature
# train <- training_common
# test <- test_common

# 0) Define the control
trControl <- trainControl(method = "repeatedcv",
                          # method = "cv",
                          number = 10,
                          repeats = 10,
                          search = "grid",
                          # summaryFunction = twoClassSummary,
                          # classProbs = TRUE,
                          # savePredictions = TRUE
                          )

# 1) Search best mtry
tuneGrid <- expand.grid(.mtry = c(1,2,3,4))
# tuneGrid <- expand.grid(.mtry = c(2,3,4,5,6,7,8,9,10,11,12))
# tuneGrid <- expand.grid(.mtry = c(2,3,4,5,6,7,8,9,10,11))
rf_default <- train(class~.,
                    data = train,
                    method = "rf",
                    metric = "Accuracy",
                    trControl = trControl,
                    tuneGrid = tuneGrid,
                    importance = TRUE)
print(rf_default)
best_mtry <- rf_default$bestTune$mtry 
best_mtry 
# 3 (56%)

# 2) Search best ntree
tuneGrid <- expand.grid(.mtry = best_mtry)
store_maxtrees <- list()
for (ntree in c(250, 300, 350, 400, 450, 500, 550, 
                600, 800, 1000, 2000, 3000)) {
  # Run the model
  rf_maxtrees <- train(class~.,
                      data = train,
                      method = "rf",
                      metric = "Accuracy",
                      trControl = trControl,
                      tuneGrid = tuneGrid,
                      ntree = ntree,
                      importance = TRUE)
  key <- toString(ntree)
  store_maxtrees[[key]] <- rf_maxtrees
}
results_tree <- resamples(store_maxtrees)
summary(results_tree) # 550 (54%)

# colnames(training_common)
# colnames(test_common)

# 3) Train with best settings
fit_rf <- train(class~.,
                train,
                method = "rf",
                metric = "Accuracy",
                # metric = "ROC",
                tuneGrid = tuneGrid,
                trControl = trControl,
                importance = TRUE,
                ntree = 250)
# fit_rf <- train(class~.,
#                 training_common,
#                 method = "rf",
#                 # metric = "Accuracy",
#                 metric = "ROC",
#                 tuneGrid = tuneGrid,
#                 trControl = trControl,
#                 importance = TRUE,
#                 ntree = 400)

# 4) Evaluate
prediction <- predict(fit_rf, test)
confusion_matrix <- confusionMatrix(prediction, test$class)
# prediction <- predict(fit_rf, my_data_gene_exp_vs_outcome)
# confusion_matrix <- confusionMatrix(prediction, my_data_gene_exp_vs_outcome$class)
confusion_matrix$overall
confusion_matrix$byClass
confusion_matrix$table

# Plot variable importance
varImp(fit_rf)
plot(varImp(fit_rf))

# rf_model_info <- getModelInfo(model = "rf")
# rf_model_info$rfRules


# library(pROC)
# 
# plot(roc(predictor = fit_rf$pred$X1, response = fit_rf$pred$obs))
