# Load necessary packages
library(ggplot2)
library(tidyr)

# Define the base directory containing the seed_n folders
base_dir <- "Simulations/Dimension known/Errors"

# Retrieve all seed directories
seed_dirs <- list.files(base_dir, pattern = "^seed_")

# Initialize lists to store errors
errors_Bhat_p4 <- list()
errors_gamma_p4 <- list()
errors_Bhat_p30 <- list()
errors_gamma_p30 <- list()

# Loop through each seed directory to load errors
for (seed_dir in seed_dirs) {
  # Load each error file for the current seed
  error_files <- list.files(file.path(base_dir, seed_dir), pattern = "^Errors_", full.names = TRUE)
  
  for (error_file in error_files) {
    load(error_file)  # Load the `Errors` list
    
    if (length(Errors$Gardes$Bhat) != 0) {
      # Case for p = 4
      errors_Bhat_p4 <- lapply(c("CTI","Gardes","TIREX1","TIREX2"), function(model) rbind(errors_Bhat_p4[[model]], Errors[[model]]$Bhat))
      errors_gamma_p4 <- lapply(names(Errors), function(model) rbind(errors_gamma_p4[[model]], Errors[[model]]$gamma))
    } else {
      # Case for p = 30
      errors_Bhat_p30 <- lapply(c("CTI","TIREX1","TIREX2"), function(model) c(errors_Bhat_p30[[model]], Errors[[model]]$Bhat))
      errors_gamma_p30 <- lapply(names(Errors), function(model) c(errors_gamma_p30[[model]], Errors[[model]]$gamma))
    }
  }
}

# Compute the means for p = 4
means_Bhat_p4 <- lapply(errors_Bhat_p4, function(matrix) apply(matrix, 2, mean, na.rm = TRUE))
means_gamma_p4 <- lapply(errors_gamma_p4, function(matrix) apply(matrix, 2, mean, na.rm = TRUE))

# Create heatmaps for each model with p = 4
for (model in names(means_Bhat_p4)) {
  # Convert the mean error matrix to a data frame for ggplot2
  mean_Bhat_df <- melt(as.data.frame(means_Bhat_p4[[model]]))
  mean_gamma_df <- melt(as.data.frame(means_gamma_p4[[model]]))
  
  # Heatmap for Bhat error
  ggplot(mean_Bhat_df, aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(value, 2)), color = "black") +
    scale_fill_gradient(low = "white", high = "black") +
    labs(title = paste("Mean Error Bhat for Model", model), x = "Alpha Exponents", y = "h Exponents") +
    theme_minimal()
  
  # Heatmap for gamma error
  ggplot(mean_gamma_df, aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(value, 2)), color = "black") +
    scale_fill_gradient(low = "white", high = "black") +
    labs(title = paste("Mean Error Gamma for Model", model), x = "Alpha Exponents", y = "h Exponents") +
    theme_minimal()
}

# Find the optimal hyperparameters for each model (p = 4)
min_params_Bhat_p4 <- lapply(means_Bhat_p4, function(matrix) {
  which(matrix == min(matrix), arr.ind = TRUE)
})
min_params_gamma_p4 <- lapply(means_gamma_p4, function(matrix) {
  which(matrix == min(matrix), arr.ind = TRUE)
})

# Box plots for the best hyperparameter pairs (p = 4)
for (model in names(min_params_Bhat_p4)) {
  alpha_idx <- min_params_Bhat_p4[[model]][1]
  h_idx <- min_params_Bhat_p4[[model]][2]
  
  # Box plot of Bhat error
  ggplot(data.frame(value = errors_Bhat_p4[[model]][, alpha_idx, h_idx]), aes(y = value)) +
    geom_boxplot() +
    labs(title = paste("Boxplot of Bhat Error for Model", model, "at Optimal Parameters"),
         y = "Error Bhat") +
    theme_minimal()
  
  # Box plot of gamma error
  alpha_idx <- min_params_gamma_p4[[model]][1]
  h_idx <- min_params_gamma_p4[[model]][2]
  ggplot(data.frame(value = errors_gamma_p4[[model]][, alpha_idx, h_idx]), aes(y = value)) +
    geom_boxplot() +
    labs(title = paste("Boxplot of Gamma Error for Model", model, "at Optimal Parameters"),
         y = "Error Gamma") +
    theme_minimal()
}

# Box plots for p = 30
for (model in names(errors_Bhat_p30)) {
  # Box plot of Bhat error
  ggplot(data.frame(value = errors_Bhat_p30[[model]]), aes(y = value)) +
    geom_boxplot() +
    labs(title = paste("Boxplot of Bhat Error for Model", model, "when p = 30"),
         y = "Error Bhat") +
    theme_minimal()
  
  # Box plot of gamma error
  ggplot(data.frame(value = errors_gamma_p30[[model]]), aes(y = value)) +
    geom_boxplot() +
    labs(title = paste("Boxplot of Gamma Error for Model", model, "when p = 30"),
         y = "Error Gamma") +
    theme_minimal()
}