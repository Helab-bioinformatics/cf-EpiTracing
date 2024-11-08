  #Training codes for XGBoost tumor burden evaluation model
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
  
  Tumor_burden_training<-read.table("Tumor_burden_training.txt")
  Tumor_burden_validation<-read.table("Tumor_burden_validation.txt")
  trainData<-Tumor_burden_training[,c(1:408)]
  trainLabel<-Tumor_burden_training[,409]
  ValidationData<-Tumor_burden_validation[,c(1:408)]
  ValidationLabel<-Tumor_burden_validation[,409]

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
    )
    
    cv <- xgb.cv(
      params = params,
      data = xgb.DMatrix(data = trainData, label = trainLabel),
      nfold = 10,  
      nrounds = as.integer(nrounds),
      early_stopping_rounds = 50,  
      verbose = 0
    )
    
    list(Score = min(cv$evaluation_log$test_rmse_mean), Pred = NULL)
  }
  
  bounds <- list(
    max_depth = c(3L, 15L),  
    eta = c(0.0001, 0.3),  
    gamma = c(0, 10),  
    colsample_bytree = c(0.4, 1),
    min_child_weight = c(1L, 20L),  
    subsample = c(0.5, 1),
    alpha = c(0, 20),
    lambda = c(0, 20),
    scale_pos_weight = c(1, 10),
    nrounds = c(100L, 1500L) 
  )
  
  opt_results <- BayesianOptimization(
    FUN = bayes_eval,
    bounds = bounds,
    init_points = 10,
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
    lambda = opt_results$Best_Par["lambda"]
  )

  bst_model <- xgb.train(
    params = best_params,
    data = xgb.DMatrix(data = trainData, label = trainLabel),
    nrounds = as.integer(opt_results$Best_Par["nrounds"]),
    watchlist = list(train = xgb.DMatrix(data = trainData, label = trainLabel),validate = xgb.DMatrix(data = ValidationData, label = ValidationLabel)),
    print_every_n = 1,
  )
  
  #Training curves
  train_rmse <- bst_model$evaluation_log$train_rmse[1:100]
  test_rmse <- bst_model$evaluation_log$train_rmse[1:100]
  iterations <- 1:length(train_rmse)
  
  df_loss <- data.frame(
    Iteration = iterations,
    Train_RMSE = train_rmse,
    Test_RMSE = test_rmse
  )

  pdf("20241108_Regression_bias_training_curves_top100.pdf")
  ggplot(df_loss, aes(x = Iteration)) +
    geom_line(aes(y = Train_RMSE, color = "Train RMSE")) +
    geom_line(aes(y = Test_RMSE, color = "Test RMSE")) +
    labs(title = "Loss Curve (Train vs Test RMSE)", x = "Iteration", y = "RMSE") +
    scale_color_manual(values = c("Train RMSE" = "blue", "Test RMSE" = "red")) +
    theme_minimal()
  dev.off()
  
  Validation_tumorBurden<-as.data.frame(ValidationLabel)
  Validation_tumorBurden$Group<-"Validation"
  Validation_tumorBurden$label<-c(rep("Healthy",32),rep("CRC",18))
  preds_validation <- predict(bst_model, xgb.DMatrix(ValidationData))
  Validation_tumorBurden$Prediction<-preds_validation
  colnames(Validation_tumorBurden)[1]<-"True"
  
  Train_tumorBurden<-as.data.frame(trainLabel)
  Train_tumorBurden$Group<-"Training"
  Train_tumorBurden$label<-c(rep("Healthy",93),rep("CRC",36))
  preds_train <- predict(bst_model, xgb.DMatrix(trainData))
  Train_tumorBurden$Prediction<-preds_train
  colnames(Train_tumorBurden)[1]<-"True"
  
  Train_tumorBurden$Prediction[which(Train_tumorBurden$Prediction<1)]<-0
  cor_training <- cor.test(Train_tumorBurden$True,Train_tumorBurden$Prediction,method = "spearman")
  correlation_coefficient_training <- cor_training$estimate
  p_value_training <- cor_training$p.value
  rmse_training <- sqrt(mean((Train_tumorBurden$True - Train_tumorBurden$Prediction) ^ 2))
  mae_training <- mean(abs(Train_tumorBurden$Prediction - Train_tumorBurden$True))
  print(correlation_coefficient_training)
  print(p_value_training)
  print(rmse_training)
  print(mae_training)
  
  Validate_tumorBurden<-Validation_tumorBurden
  Validate_tumorBurden$Prediction[which(Validate_tumorBurden$Prediction<1)]<-0
  cor_validation <- cor.test(Validate_tumorBurden$True,Validate_tumorBurden$Prediction,method = "spearman")
  correlation_coefficient_validation <- cor_validation$estimate
  p_value_validation <- cor_validation$p.value
  rmse_validation <- sqrt(mean((Validate_tumorBurden$True - Validate_tumorBurden$Prediction) ^ 2))
  mae_validation <- mean(abs(Validate_tumorBurden$Prediction - Validate_tumorBurden$True))
  print(correlation_coefficient_validation)
  print(p_value_validation)
  print(rmse_validation)
  print(mae_validation)
  
  #Correlation plot
  library(ggExtra)
  Whole_tumorBurden<-rbind(Train_tumorBurden,Validate_tumorBurden)
  pdf("CorrelationPlot_xgbPrediction_tumorBurden.pdf")
  p<-ggplot(data = Whole_tumorBurden, mapping = aes(x = True, y = Prediction,group=Group)) + geom_point(data = Whole_tumorBurden, aes(size = 3, colour=Whole_tumorBurden$Group))+ stat_smooth(method = 'lm')
  ggMarginal(p,type="density",groupColour = TRUE,groupFill = TRUE)
  dev.off()
  