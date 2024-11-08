####Training codes for CRC-Healthy XGBoost classification and CRA detection model
library(xgboost)
library(caret)
library(ParBayesianOptimization)
library(e1071)
library(pROC)
library(DMwR)
library(ggplot2)
library(lattice)
library(pROC)
library(rBayesianOptimization)

trainData<-read.table("TrainingGroup_matrix_CRC_classification.txt")
trainLabel <- c(rep(0, 93), rep(1, 75))  

validationData<-read.table("ValidationGroup_matrix_CRC_classification.txt")
ValidationLabel <- c(rep(0, 32), rep(1, 32))  

# Adenoma samples
adenoma_samples <- read.table("Adenoma_matrix_CRC_classification.txt")

bayes_eval <- function(max_depth, eta, gamma, colsample_bytree, min_child_weight, subsample, alpha, lambda, scale_pos_weight, nrounds) {
  params <- list(
    objective = "binary:logistic",  
    eval_metric = "auc",            
    max_depth = as.integer(max_depth),
    eta = eta,
    gamma = gamma,
    colsample_bytree = colsample_bytree,
    min_child_weight = as.integer(min_child_weight),
    subsample = subsample,
    alpha = alpha,
    lambda = lambda,
    scale_pos_weight = scale_pos_weight
  )
  
  cv <- xgb.cv(
    params = params,
    data = xgb.DMatrix(data = trainData, label = trainLabel),
    nfold = 10,  
    nrounds = as.integer(nrounds),
    early_stopping_rounds = 50,  
    verbose = 0
  )
  
  list(Score = max(cv$evaluation_log$test_auc_mean), Pred = NULL)
}

bounds <- list(
  max_depth = c(3L, 10L),
  eta = c(0.001, 0.3),  
  gamma = c(0, 5),
  colsample_bytree = c(0.5, 1),
  min_child_weight = c(1L, 10L),
  subsample = c(0.5, 1),
  alpha = c(0, 10),
  lambda = c(0, 10),
  scale_pos_weight = c(1, 10),  
  nrounds = c(100L, 1000L)  
)

opt_results <- BayesianOptimization(
  FUN = bayes_eval,
  bounds = bounds,
  init_points = 10,
  n_iter = 30  
)

best_params <- list(
  objective = "binary:logistic",
  eval_metric = "auc",  
  max_depth = as.integer(opt_results$Best_Par["max_depth"]),
  eta = opt_results$Best_Par["eta"],
  gamma = opt_results$Best_Par["gamma"],
  colsample_bytree = opt_results$Best_Par["colsample_bytree"],
  min_child_weight = as.integer(opt_results$Best_Par["min_child_weight"]),
  subsample = opt_results$Best_Par["subsample"],
  alpha = opt_results$Best_Par["alpha"],
  lambda = opt_results$Best_Par["lambda"],
  scale_pos_weight = opt_results$Best_Par["scale_pos_weight"]  
)

  bst_model <- xgb.train(
    params = best_params,
    data = xgb.DMatrix(data = trainData, label = trainLabel),
    nrounds = as.integer(opt_results$Best_Par["nrounds"]),
    watchlist = list(train = xgb.DMatrix(data = trainData, label = trainLabel), test = xgb.DMatrix(data = ValidationData, label = ValidationLabel)),
    print_every_n = 1,
  )
  

  #Confusion matrix and performance metrics
  #Training group
  preds_train <- predict(bst_model, xgb.DMatrix(trainData))
  pred_labels_train <- ifelse(preds_train > 0.5, 1, 0)
  conf_matrix_train <- confusionMatrix(factor(pred_labels_train), factor(trainLabel))
  print(conf_matrix_train)
  accuracy_train <- conf_matrix_train$overall["Accuracy"]
  precision_train <- conf_matrix_train$byClass["Pos Pred Value"]
  recall_train <- conf_matrix_train$byClass["Sensitivity"]
  f1_train <- 2 * (precision * recall) / (precision + recall)
  
  conf_matrix_train<-as.matrix(conf_matrix_train)
  conf_matrix_train[,1]<-conf_matrix_train[,1]/93
  conf_matrix_train[,2]<-conf_matrix_train[,2]/75  
  library(pheatmap)
  pdf("ConfusionMatrix_CRC_train.pdf",height=20,width=20)
  ph<-pheatmap::pheatmap(conf_matrix_train, 
                         cluster_cols = FALSE, cluster_rows = FALSE,
                         show_colnames =T, fontsize = 20,show_rownames = T,
                         fontsize_row = 20,display_numbers = TRUE,
                         border=FALSE,color = colorRampPalette(c(rep("white",1),rep("#cb181c",1)))(500))
  ph
  dev.off() 
  
  #Validation group
  preds_validation <- predict(bst_model, xgb.DMatrix(ValidationData))
  pred_labels_validation <- ifelse(preds_validation > 0.5, 1, 0)
  conf_matrix_validation <- confusionMatrix(factor(pred_labels_validation), factor(ValidationLabel))
  print(conf_matrix_validation)
  accuracy_validation <- conf_matrix_validation$overall["Accuracy"]
  precision_validation <- conf_matrix_validation$byClass["Pos Pred Value"]
  recall_validation <- conf_matrix_validation$byClass["Sensitivity"]
  f1_validation <- 2 * (precision * recall) / (precision + recall)
  
  conf_matrix_validation<-as.matrix(conf_matrix_validation)
  conf_matrix_validation[,1]<-conf_matrix_validation[,1]/32
  conf_matrix_validation[,2]<-conf_matrix_validation[,2]/32
  pdf("ConfusionMatrix_CRC_validation.pdf",height=20,width=20)
  ph<-pheatmap::pheatmap(conf_matrix_validation, 
                         cluster_cols = FALSE, cluster_rows = FALSE,
                         show_colnames =T, fontsize = 20,show_rownames = T,
                         fontsize_row = 20,display_numbers = TRUE,
                         border=FALSE,color = colorRampPalette(c(rep("white",1),rep("#cb181c",1)))(500))
  ph
  dev.off() 
  
  #Training curves (top100 iterations)
  train_logloss <- bst_model$evaluation_log$train_auc
  train_logloss<-train_logloss[c(1:100)]
  eval_logloss <- bst_model$evaluation_log$test_auc
  eval_logloss<-eval_logloss[c(1:100)]
  pdf("Loss_curve_top100Epoch.pdf")
  plot(1:100, train_logloss, type = "l", col = "blue", ylim = range(c(train_logloss, eval_logloss)), 
       xlab = "Iteration", ylab = "Logloss", main = "Training and Validation Logloss")
  lines(1:100, eval_logloss, col = "red")
  legend("topright", legend = c("Train Logloss", "Validation Logloss"), col = c("blue", "red"), lty = 1)
  dev.off() 
  
  
  #Calculate CRC scores on colorectal adenoma samples 
  preds_adenoma <- predict(bst_model, xgb.DMatrix(adenoma_samples))
  Detection_number_CRA<-length(which(preds_adenoma>=0.5))