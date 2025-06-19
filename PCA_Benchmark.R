#Creating the benchmark PCA model + evaluation of predictive performance

#Creating the model
library(dfms)
library(pheatmap)
data_rv <- read.csv('C:/Users/justj/Documents/Thesis files/realised_volatility_data.csv', row.names=1)
data_rv <- as.matrix(data_rv)

for (i in colnames(data_rv)) {
  data_rv[, i] <- scale(data_rv[, i])
}

cor_data <- cor(data_rv)

pheatmap(cor_data,
         main = 'Correlation Matrix Heatmap',
         display_numbers = FALSE,
         fontsize_row = 8,
         fontsize_col = 8)


eigen_data <- eigen(cor_data)
eigenvectors <- eigen_data$vectors
print(eigenvectors)
eigenvalues <- eigen_data$values
print(eigenvalues)

variance_explained <- eigenvalues / sum(eigenvalues)
cumulative_variance <- cumsum(variance_explained)
print(data.frame(Eigenvalue = eigenvalues, Variance = variance_explained, Cumulative = cumulative_variance))
plot(eigenvalues, type = 'b', pch = 19, 
     xlab = 'Principal Component', 
     ylab = 'Eigenvalue', 
     main = 'Scree Plot')
pc_scores <- as.matrix(data_rv) %*% eigenvectors

ICr(data_rv)

#Contributions to PC1
loadings <- eigenvectors[, 1]  
print(loadings)
variable_names <- colnames(data_rv)  
loadings <- setNames(loadings, variable_names)
loadings <- sort(abs(loadings), decreasing = TRUE)  
print(loadings)
dominant_variables <- data.frame(
  Variable = names(loadings),
  Loading = loadings
)
print(dominant_variables)

#Barplot of loadings
barplot(
  loadings,
  names.arg = names(loadings),
  las = 2, 
  col = 'thistle',
  main = 'Variable Contributions to PC1',
  xlab = 'Variables',
  ylab = 'Loadings (Absolute Value)'
)

#Barplot of percentage contribution
squared_loadings <- loadings^2 
contribution_percentage <- (squared_loadings / sum(squared_loadings)) * 100 
sorted_contribution <- sort(contribution_percentage, decreasing = TRUE)

contribution_df <- data.frame(
  Variable = names(sorted_contribution),
  Contribution = sorted_contribution
)
print(contribution_df)

barplot(
  sorted_contribution,
  names.arg = names(sorted_contribution),
  las = 2,
  col = 'thistle',
  main = 'Percentage Contribution to PC1',
  xlab = 'Variables',
  ylab = 'Percentage Contribution'
)

#Contributions to PC2
loadings <- eigenvectors[, 2]  
print(loadings)
variable_names <- colnames(data_rv)  
loadings <- setNames(loadings, variable_names)
loadings <- sort(abs(loadings), decreasing = TRUE)  
print(loadings)
dominant_variables <- data.frame(
  Variable = names(loadings),
  Loading = loadings
)
print(dominant_variables)

#Barplot of loadings
barplot(
  loadings,
  names.arg = names(loadings),
  las = 2, 
  col = 'cyan',
  main = 'Variable Contributions to PC2',
  xlab = 'Variables',
  ylab = 'Loadings (Absolute Value)'
)

#Barplot of percentage contribution
squared_loadings <- loadings^2 
contribution_percentage <- (squared_loadings / sum(squared_loadings)) * 100 
sorted_contribution <- sort(contribution_percentage, decreasing = TRUE)

contribution_df <- data.frame(
  Variable = names(sorted_contribution),
  Contribution = sorted_contribution
)
print(contribution_df)

barplot(
  sorted_contribution,
  names.arg = names(sorted_contribution),
  las = 2,
  col = 'cyan',
  main = 'Percentage Contribution to PC2',
  xlab = 'Variables',
  ylab = 'Percentage Contribution'
)

pc_data <- as.data.frame(scale(pc_scores[, 1:2]))

plot(pc_data[,1], type = 'l', col = 'black', lwd = 1,
     main = 'PC1 over Time',
     xlab = 'Quarter',
     ylab = 'PC1')

plot(pc_data[,2], type = 'l', col = 'black', lwd = 1,
     main = 'PC2 over Time',
     xlab = 'Quarter',
     ylab = 'PC2')


#Reconstruction of real data from 2 PCs
top_scores <- pc_scores[, 1:2]            
top_eigenvectors <- eigenvectors[, 1:2] 
X_approx_std <- top_scores %*% t(top_eigenvectors)
data_raw <- read.csv('C:/Users/justj/Documents/Thesis files/realised_volatility_data.csv', row.names=1)
data_raw <- as.matrix(data_raw)
orig_means <- colMeans(data_raw)
orig_sds <- apply(data_raw, 2, sd)
X_approx_unscaled <- sweep(X_approx_std, 2, orig_sds, FUN = '*')
X_approx_unscaled <- sweep(X_approx_unscaled, 2, orig_means, FUN = '+')


#Predictive performance
data_rv <- read.csv('C:/Users/justj/Documents/Thesis files/realised_volatility_data.csv', row.names = 1)
data_rv <- as.matrix(data_rv)
y <- as.matrix(data_rv[1:240, ])       
y_test <- as.matrix(data_rv[241:298, ])  
h <- nrow(y_test)
y_scaled <- scale(y)
means <- attr(y_scaled, 'scaled:center')
sds <- attr(y_scaled, 'scaled:scale')

cor_data <- cor(y_scaled)
eigen_data <- eigen(cor_data)
eigenvectors <- eigen_data$vectors
eigenvalues <- eigen_data$values
pc_scores_train <- y_scaled %*% eigenvectors
y_test_scaled <- scale(y_test, center = means, scale = sds)
pc_scores_test <- y_test_scaled %*% eigenvectors
top_eigenvectors <- eigenvectors[, 1:2]
X_train_approx_std <- pc_scores_train[, 1:2] %*% t(top_eigenvectors)
X_test_approx_std <- pc_scores_test[, 1:2] %*% t(top_eigenvectors)
X_train_approx <- sweep(X_train_approx_std, 2, sds, '*')
X_train_approx <- sweep(X_train_approx, 2, means, '+')
X_test_approx <- sweep(X_test_approx_std, 2, sds, '*')
X_test_approx <- sweep(X_test_approx, 2, means, '+')
reconstruction_rmse <- sqrt(colMeans((y_test - X_test_approx)^2))
reconstruction_mae <- colMeans(abs(y_test - X_test_approx))

print(reconstruction_rmse)
print(reconstruction_mae)

