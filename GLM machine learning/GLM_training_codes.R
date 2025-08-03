###CRC-non CRC GLM classifier
# Set working directory
setwd("/path/to/your/data")  # Please specify your working directory

# Load pre-selected feature data matrix (ensure the file exists)
data_matrix <- read.table("your_data_file.txt", header = TRUE)  # Replace with the actual data file name

# Load necessary libraries
install.packages(c("caret", "pROC", "caTools", "glmnet"))  # Install required packages if not installed
library(caret)
library(pROC)
library(caTools)
library(glmnet)

# Define indices for different sample types (modify as per your data)
healthy_indices <- 1:a  # Replace with actual indices for healthy samples
crc_indices <- (a+1):b   # Replace with actual indices for CRC samples
chd_indices <- (b+1):c   # Replace with actual indices for CHD samples
lymphoma_indices <- (c+1):d  # Replace with actual indices for lymphoma samples

# Split data into training and testing sets (80% training, 20% testing)
train_healthy <- sample(healthy_indices, size = floor(0.8 * length(healthy_indices)))
train_crc <- sample(crc_indices, size = floor(0.8 * length(crc_indices)))
train_chd <- sample(chd_indices, size = floor(0.8 * length(chd_indices)))
train_lymphoma <- sample(lymphoma_indices, size = floor(0.8 * length(lymphoma_indices)))

# Combine training indices and define test indices (remaining data points)
train_indices <- c(train_healthy, train_crc, train_chd, train_lymphoma)
test_indices <- setdiff(1:d, train_indices)  # Adjust '564' to the actual number of samples in your dataset

# Subset the data matrix for training and testing
train_matrix <- data_matrix[, train_indices]
test_matrix <- data_matrix[, test_indices]

# Define the sample types for the training data
# Replace "Healthy", "CRC", etc., based on your data's categories
sampletype <- c(rep("Healthy", length(train_healthy)), 
                rep("CRC", length(train_crc)), 
                rep("Healthy", length(healthy_indices) - length(train_healthy)),
                rep("Healthy", length(chd_indices)), 
                rep("Lymphoma", length(lymphoma_indices)))

# Prepare the training and test data matrices
data_train <- as.data.frame(t(train_matrix))
data_test <- as.data.frame(t(test_matrix))

# Dynamically assign the 'condition' labels (e.g., "Healthy", "CRC", etc.)
data_train$condition <- factor(sampletype[train_indices])
data_test$condition <- factor(sampletype[test_indices])

# Prepare the model input (X) and output (y)
x_train <- as.matrix(data_train[, -which(names(data_train) == "condition")])
y_train <- data_train$condition
x_test <- as.matrix(data_test[, -which(names(data_test) == "condition")])
y_test <- data_test$condition

# Train a Lasso logistic regression model (alpha = 1 for Lasso regularization)
model_lasso <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1)

# Output the best lambda value (regularization strength)
best_lambda <- model_lasso$lambda.min
cat("Best lambda:", best_lambda, "\n")

# Make predictions on test and train datasets
test_pred_probs <- predict(model_lasso, newx = x_test, type = "response", s = "lambda.min")
train_pred_probs <- predict(model_lasso, newx = x_train, type = "response", s = "lambda.min")

# Calculate ROC and AUC for both train and test sets
roc_result_test <- roc(y_test, test_pred_probs)
roc_result_train <- roc(y_train, train_pred_probs)

# Calculate sensitivity and specificity from confusion matrix
pred_labels_test <- ifelse(test_pred_probs > 0.5, 1, 0)
conf_matrix_test <- confusionMatrix(factor(pred_labels_test), factor(y_test))

# Calculate AUC confidence intervals using bootstrap method
auc_ci1 <- ci(roc_result_test, method = "bootstrap")
auc_ci2 <- ci(roc_result_train, method = "bootstrap")

# Print AUC and confidence intervals
cat("Test-95% CI for AUC: [", auc_ci1[1], ", ", auc_ci1[2], "]\n", sep = "")
cat("Train-95% CI for AUC: [", auc_ci2[1], ", ", auc_ci2[2], "]\n", sep = "")

# Plot ROC curve
pdf("ROC_Curve_CRC_binary.pdf", height = 20, width = 40)
plot(roc_result_train, main = "ROC Curve", col = "blue", lwd = 2)
plot(roc_result_test, add = TRUE, col = "red", lwd = 2)
abline(a = 0, b = 1, col = "black", lty = 2)
legend("bottomright", legend = c("Train ROC", "Test ROC"), col = c("blue", "red"), lwd = 2)
dev.off()

