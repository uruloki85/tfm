library(randomForest)

# Find best mtry for each signature

# Treatment vs outcome
mtry_treatment <- tuneRF(my_data_treatment_vs_outcome[,1:ncol(my_data_treatment_vs_outcome)-1], 
                         my_data_treatment_vs_outcome$class, 
                         ntreeTry=500,
                         stepFactor=1.5, 
                         improve=0.001, 
                         trace=TRUE, 
                         plot=TRUE)
best_mtry_treatment <- mtry_treatment[mtry_treatment[, 2] == min(mtry_treatment[, 2]), 1]
mtry_treatment
best_mtry_treatment

# Gene expression vs outcome
mtry_gene_exp <- tuneRF(my_data_gene_exp_vs_outcome[,1:ncol(my_data_gene_exp_vs_outcome)-1], 
                        my_data_gene_exp_vs_outcome$class, 
                        ntreeTry=500,
                        stepFactor=1.5, 
                        improve=0.001, 
                        trace=TRUE, 
                        plot=TRUE)
best_mtry_gene_exp <- mtry_gene_exp[mtry_gene_exp[, 2] == min(mtry_gene_exp[, 2]), 1]
mtry_gene_exp
best_mtry_gene_exp

# Common signature
mtry_common <- tuneRF(my_data_signature_final[,1:ncol(my_data_signature_final)-1], 
                      my_data_signature_final$class, 
                      ntreeTry=1500,
                      stepFactor=1.5, 
                      improve=0.001, 
                      trace=TRUE, 
                      plot=TRUE)
best_mtry_common <- mtry_common[mtry_common[, 2] == min(mtry_common[, 2]), 1]
mtry_common
best_mtry_common # 4, 6

# Build model using this mtry

# Treatment vs outcome
rf_common <- randomForest(class ~ ., data=my_data_signature_final, ntree=1500, mtry=6, importance=TRUE)
rf_common
# varImpPlot(rf_classifier)

#predicting the class for the test data set
pred_common <- predict(rf_common,validation_common,type="class")

#creating confusion matrix
table(pred_common,validation_common$class) 
rf_common$importance
varImpPlot(rf_common)


# Performance
pred1 <- predict(rf_common, type = "prob")

library(ROCR)

perf = prediction(pred1[,2], training_common$class)
# 1. Area under curve
auc = performance(perf, "auc")
auc
# 2. True Positive and Negative Rate
pred3 = performance(perf, "tpr","fpr")
# 3. Plot the ROC curve
plot(pred3,main="ROC Curve for Random Forest",col=2,lwd=2)
abline(a=0,b=1,lwd=2,lty=2,col="gray")
