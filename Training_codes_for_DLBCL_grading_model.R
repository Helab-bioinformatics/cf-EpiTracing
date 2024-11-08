  library(xgboost)
  library(caret)
  library(ParBayesianOptimization)
  library(e1071)
  library(pROC)
  library(DMwR)
  library(ggplot2)
  library(lattice)
  library(rBayesianOptimization)
  library(data.table)
  library(dplyr)
  
  DLBCL_grading_training<-read.table("DLBCL_grading_training.txt")
  DLBCL_grading_validation<-read.table("DLBCL_grading_validation.txt")
  train_data<-DLBCL_grading_training[,c(1:408)]
  train_label<-DLBCL_grading_training[,409]
  validation_data<-DLBCL_grading_validation[,c(1:408)]
  validation_label<-DLBCL_grading_validation[,409]
  train_matrix <- xgb.DMatrix(data = as.matrix(train_data), label = train_label)
  validation_matrix <- xgb.DMatrix(data = as.matrix(validation_data), label = validation_label)
  
  bayes_eval <- function(max_depth, eta, gamma, colsample_bytree, min_child_weight, subsample, alpha, lambda, scale_pos_weight, nrounds) {
    params <- list(
      objective = "reg:squarederror", 
      eval_metric = "rmse",              
      max_depth = as.integer(max_depth),
      eta = eta,
      gamma = gamma,
      colsample_bytree = colsample_bytree,
      min_child_weight = as.integer(min_child_weight),
      subsample = subsample,
      alpha = alpha,
      lambda = lambda,
      scale_pos_weight = scale_pos_weight,
    )
    
    cv <- xgb.cv(
      params = params,
      data = train_matrix,
      nfold = 10,  
      nrounds = as.integer(nrounds),
      early_stopping_rounds = 50,  
      verbose = 0
    )
        list(Score = min(cv$evaluation_log$test_rmse_mean), Pred = NULL)
  }
  
  bounds <- list(
    max_depth = c(3L, 10L),
    eta = c(0.01, 0.3),
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
    init_points = 20,
    n_iter = 30  
  )
  
  best_params <- list(
    objective = "reg:squarederror",
    eval_metric = "rmse",
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
    data = train_matrix,
    nrounds = as.integer(opt_results$Best_Par["nrounds"]),
    watchlist = list(train = train_matrix,test=validation_matrix),
    print_every_n = 1,
  )
  
  preds_train <- predict(bst_model, train_matrix)
  train_rmse <- sqrt(mean((train_label - preds_train)^2))  
  train_mae <- mean(abs(preds_train - train_label))
  train_correlation<-cor(train_label, preds_train, method = "spearman")
  
  preds_validation<-predict(bst_model, validation_matrix)
  validation_rmse <- sqrt(mean((validation_label - preds_validation)^2))  
  validation_mae <- mean(abs(preds_validation - validation_label))
  validation_correlation<-cor(validation_label, preds_validation, method = "spearman")

  #Top 100 iteration training curves
  train_logloss <- bst_model$evaluation_log$train_rmse
  train_logloss<-train_logloss[c(1:100)]
  eval_logloss <- bst_model$evaluation_log$test_rmse
  eval_logloss<-eval_logloss[c(1:100)]
  pdf("Loss_curve_top10.pdf")
  plot(1:100, train_logloss, type = "l", col = "blue", ylim = range(c(train_logloss, eval_logloss)), 
       xlab = "Iteration", ylab = "Logloss", main = "Training and Validation Logloss")
  lines(1:100, eval_logloss, col = "red")
  legend("topright", legend = c("Train Logloss", "Validation Logloss"), col = c("blue", "red"), lty = 1)
  dev.off()
  
  ###Classification threshold from training group samples
  Train_matrix<-matrix(ncol=2,nrow=nrow(train_matrix))  
  Train_matrix[,1]<-train_label
  Train_matrix[,2]<-predict(bst_model, train_matrix)
  rownames(Train_matrix)<-rownames(train_data)
  Train_matrix_0<-Train_matrix[which(Train_matrix[,1]==0),]
  Train_matrix_1<-Train_matrix[which(Train_matrix[,1]==1),]
  Train_matrix_2<-Train_matrix[which(Train_matrix[,1]==2),]
  Train_matrix_3<-Train_matrix[which(Train_matrix[,1]==3),]
  Train_matrix_4<-Train_matrix[which(Train_matrix[,1]==4),]
  
  Train_matrix_01<-rbind(Train_matrix_0,Train_matrix_1)
  Train_matrix_12<-rbind(Train_matrix_1,Train_matrix_2)
  Train_matrix_23<-rbind(Train_matrix_2,Train_matrix_3)
  Train_matrix_34<-rbind(Train_matrix_3,Train_matrix_4)
  
  threshold1 <- roc(Train_matrix_01[,1], Train_matrix_01[,2])
  threshold2 <- roc(Train_matrix_12[,1], Train_matrix_12[,2])
  threshold3 <- roc(Train_matrix_23[,1], Train_matrix_23[,2])
  threshold4 <- roc(Train_matrix_34[,1], Train_matrix_34[,2])
  
  best_threshold1 <- coords(threshold1, "best", ret = "threshold", transpose = FALSE)
  best_threshold2 <- coords(threshold2, "best", ret = "threshold", transpose = FALSE)
  best_threshold3 <- coords(threshold3, "best", ret = "threshold", transpose = FALSE)
  best_threshold4 <- coords(threshold4, "best", ret = "threshold", transpose = FALSE)
  
  Validation_matrix<-matrix(ncol=3,nrow=nrow(validation_data))  
  Validation_matrix[,1]<-validation_label
  Validation_matrix[,2]<-predict(bst_model, validation_matrix)
  rownames(Validation_matrix)<-rownames(validation_data)
  Validation_matrix[, 2] <- as.numeric(Validation_matrix[, 2])
  indices1 <- which(Validation_matrix[, 2] < as.numeric(best_threshold1))
  indices2 <- which(Validation_matrix[, 2] > as.numeric(best_threshold1) & Validation_matrix[, 2] < as.numeric(best_threshold2))
  indices3 <- which(Validation_matrix[, 2] > as.numeric(best_threshold2) & Validation_matrix[, 2] < as.numeric(best_threshold3))
  indices4 <- which(Validation_matrix[, 2] > as.numeric(best_threshold3) & Validation_matrix[, 2] < as.numeric(best_threshold4))
  indices5 <- which(Validation_matrix[, 2] > as.numeric(best_threshold4))
  
  Validation_matrix[indices1,3]<- "0"
  Validation_matrix[indices2,3]<- "1"
  Validation_matrix[indices3,3]<- "2"
  Validation_matrix[indices4,3]<- "3"
  Validation_matrix[indices5,3]<- "4"
  
  #confusion matrix for validation
  True_label<-as.character(Validation_matrix[,1])
  Prediction<-as.character(Validation_matrix[,3])
  True_label[which(True_label=="0")]<-"A"
  True_label[which(True_label=="1"|True_label=="2")]<-"B"
  True_label[which(True_label=="3"|True_label=="4")]<-"C"
  
  Prediction[which(Prediction=="0")]<-"A"
  Prediction[which(Prediction=="1"|Prediction=="2")]<-"B"
  Prediction[which(Prediction=="3"|Prediction=="4")]<-"C"
  
  conf_matrix <- confusionMatrix(factor(c(True_label)),factor(c(Prediction)))
  accuracy<-conf_matrix[["overall"]][["Accuracy"]]
  precision <- conf_matrix$byClass["Pos Pred Value"]
  recall <- conf_matrix$byClass["Sensitivity"]
  f1 <- 2 * (precision * recall) / (precision + recall)
  
  conf_matrix<-as.matrix(conf_matrix)
  conf_matrix[,1]<-conf_matrix[,1]/32
  conf_matrix[,2]<-conf_matrix[,2]/18
  conf_matrix[,3]<-conf_matrix[,3]/19
  
  pdf("ComfusionMatrix_validation.pdf")
  ph<-pheatmap::pheatmap(conf_matrix, 
                         cluster_cols = FALSE, cluster_rows = FALSE,
                         show_colnames =T, fontsize = 20,show_rownames = T,
                         fontsize_row = 20,display_numbers = TRUE,
                         border=FALSE,color = colorRampPalette(c(rep("white",1),rep("#cb181c",1)))(500))
  ph
  dev.off() 
