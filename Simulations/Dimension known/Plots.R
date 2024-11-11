# Load necessary packages
library(ggplot2)
library(tidyr)

# Define the base directory containing the seed_n folders
base_dir <- "Simulations/Dimension known/Errors"


# Retrieve all directories
models <- list("model1p4","model2p4","model3p4","model4p4",
                   "model1p30","model2p30","model3p30","model4p30")
seed_dirs <- list.files(base_dir, pattern = "^seed_")

# Initialize a list of lists to store the error matrices by method
red_methods_p4 <- c("CTI", "Gardes", "TIREX1", "TIREX2")
errors_Bhat_p4 <- setNames(vector("list", length(red_methods_p4)), red_methods_p4)

red_methods_p30 <- c("CTI", "TIREX1", "TIREX2")
errors_Bhat_p30 <- setNames(vector("list", length(red_methods_p30)), red_methods_p30)

methods_p4 <- c("CTI", "Gardes", "TIREX1", "TIREX2", "B0", "Id")
errors_gamma_p4 <- setNames(vector("list", length(methods_p4)), methods_p4)

methods_p30 <- c("CTI", "TIREX1", "TIREX2", "B0")
errors_gamma_p30 <- setNames(vector("list", length(methods_p30)), methods_p30)



for(model in models){
  
  # Loop through each seed directory to load errors
  for (seed_dir in seed_dirs) {
    # Load each error file for the current seed
    error_file <- paste0("Simulations/Dimension known/Errors/",seed_dir,"/Errors_",model,"-",seed_dir,".RData")
    
    load(error_file)  # Load the `Errors` list
    
    if (length(Errors$Gardes$Bhat) != 0) {
      
      ###################
      # Case for p = 4  #
      ###################
      
      # Loop through each method to store Bhat error matrices
      for (method in methods_p4) {
        if (!is.null(Errors[[method]]$Bhat)) {
          # Append the current Bhat error matrix to the list for the corresponding method
          errors_Bhat_p4[[method]] <- append(errors_Bhat_p4[[method]], list(Errors[[method]]$Bhat))
        }
        if (!is.null(Errors[[method]]$gamma)) {
          # Append the current Bhat error matrix to the list for the corresponding method
          errors_gamma_p4[[method]] <- append(errors_gamma_p4[[method]], list(Errors[[method]]$gamma))
        }
      }
    
      
      # Calculate the mean matrix for each method across all seeds
      means_Bhat_p4 <- lapply(errors_Bhat_p4, function(method_errors) {
        
        if (length(method_errors) > 1) {
          mean_matrix <- Reduce("+", method_errors) / length(method_errors)
        } else {
          mean_matrix <- method_errors[[1]]
        }
        mean_matrix
      }) 
      
      
      means_gamma_p4 <- lapply(errors_gamma_p4, function(method_errors) {
        
        if (length(method_errors) > 1) {
          mean_matrix <- Reduce("+", method_errors) / length(method_errors)
        } else {
          mean_matrix <- method_errors[[1]]
        }
        mean_matrix
      })   
      
      
      # Create heatmaps for each method with p = 4
      
      output_dir <- "Simulations/Dimension known/Plots"
      
      for (method in names(means_gamma_p4)) {
        # Convert the mean error matrix to a data frame for ggplot2
        if (!is.null(means_Bhat_p4[[method]])) {
          mean_Bhat_df_long <- data.frame(
            Var1 = rep(seq_len(nrow(means_Bhat_p4[[method]])), ncol(means_Bhat_p4[[method]])),
            Var2 = rep(seq_len(ncol(means_Bhat_p4[[method]])), each = nrow(means_Bhat_p4[[method]])),
            value = as.vector(means_Bhat_p4[[method]])
          )
          
          # Heatmap for Bhat error
          plot_Bhat <- ggplot(mean_Bhat_df_long, aes(Var1, Var2, fill = value)) +
            geom_tile(color = "white") +
            geom_text(aes(label = round(value, 3)), color = "black") + 
            scale_fill_gradient(low = "white", high = "darkgray") +
            labs(title = paste("Mean Error Bhat for method", method), x = "Alpha Exponents", y = "h Exponents") +
            theme_minimal()
          
          ggsave(
            filename = file.path(output_dir, paste0("Mean_Error_Bhat_",model,"-", method, ".png")),
            plot = plot_Bhat,
            width = 8, height = 6
          )
        }
        
        mean_gamma_df_long <- data.frame(
          Var1 = rep(seq_len(nrow(means_gamma_p4[[method]])), ncol(means_gamma_p4[[method]])),
          Var2 = rep(seq_len(ncol(means_gamma_p4[[method]])), each = nrow(means_gamma_p4[[method]])),
          value = as.vector(means_gamma_p4[[method]])
        )
        
        # Heatmap for gamma error
        plot_gamma <- ggplot(mean_gamma_df_long, aes(Var1, Var2, fill = value)) +
          geom_tile(color = "white") +
          geom_text(aes(label = round(value, 3)), color = "black") +
          scale_fill_gradient(low = "white", high = "darkgray") +
          labs(title = paste("Mean Error Gamma for method", method), x = "Alpha Exponents", y = "h Exponents") +
          theme_minimal()
        
        ggsave(
          filename = file.path(output_dir, paste0("Mean_Error_Gamma_",model,"-", method, ".png")),
          plot = plot_gamma,
          width = 8, height = 6
        )
      }
      
      
      
      
      # Find the optimal hyperparameters for each method (p = 4)
      min_params_gamma_p4 <- lapply(means_gamma_p4, function(matrix) {
        which(matrix == min(matrix), arr.ind = TRUE)
      })
      
      # Create data frames for Bhat and gamma errors with optimal hyperparameters
      all_errors_Bhat <- data.frame()
      all_errors_gamma <- data.frame()
      
      for (method in names(means_Bhat_p4)) {
        # Get the optimal alpha and h indices for Bhat and gamma
        alpha_idx <- min_params_gamma_p4[[method]][1]
        h_idx <- min_params_gamma_p4[[method]][2]
        
        # Extract Bhat errors for optimal parameters and add to the combined data frame
        errors_Bhat_best <- unlist(lapply(errors_Bhat_p4[[method]], function(matrix) matrix[alpha_idx, h_idx]))
        all_errors_Bhat <- rbind(all_errors_Bhat, data.frame(value = errors_Bhat_best, method = method, error_type = "Bhat"))
        
        # Extract gamma errors for optimal parameters and add to the combined data frame
        errors_gamma_best <- unlist(lapply(errors_gamma_p4[[method]], function(matrix) matrix[alpha_idx, h_idx]))
        all_errors_gamma <- rbind(all_errors_gamma, data.frame(value = errors_gamma_best, method = method, error_type = "Gamma"))
      }
      
      # Combine Bhat and gamma errors into a single data frame
      all_errors <- rbind(all_errors_Bhat, all_errors_gamma)
      
      # Plot all boxplots on a single graphic, with separate panels for Bhat and gamma errors
      boxplots <- ggplot(all_errors, aes(x = method, y = value, fill = method)) +
        geom_boxplot() +
        facet_wrap(~error_type, scales = "free_y") +
        labs(title = "Comparison of Bhat and Gamma Errors Across Methods at Optimal Parameters",
             y = "Error Value",
             x = "Method") +
        theme_minimal() +
        theme(legend.position = "none")
      
      ggsave(
        filename = file.path(output_dir, paste0("Boxplots-",model,".png")),
        plot = boxplots,
        width = 8, height = 6
      )  
      
      
      
      
    } else {
      
      ###################
      # Case for p = 30 #
      ###################
      
      # Loop through each method to store Bhat error matrices
      for (method in methods_p30) {
        if (!is.null(Errors[[method]]$Bhat)) {
          # Append the current Bhat error matrix to the list for the corresponding method
          errors_Bhat_p30[[method]] <- Errors[[method]]$Bhat
        }
        if (!is.null(Errors[[method]]$gamma)) {
          # Append the current Bhat error matrix to the list for the corresponding method
          errors_gamma_p30[[method]] <- Errors[[method]]$gamma
        }
      }
    
      
  # Create data frames for Bhat and gamma errors with optimal hyperparameters
  all_errors_Bhat <- data.frame()
  all_errors_gamma <- data.frame()
  
  for (method in names(errors_Bhat_p30)) {
    
    all_errors_Bhat <- rbind(all_errors_Bhat, data.frame(value = errors_Bhat_p30[[method]], method = method, error_type = "Bhat"))
    
    all_errors_gamma <- rbind(all_errors_gamma, data.frame(value = errors_gamma_p30[[method]], method = method, error_type = "Gamma"))
  }
  
  # Combine Bhat and gamma errors into a single data frame
  all_errors <- rbind(all_errors_Bhat, all_errors_gamma)
  
  # Plot all boxplots on a single graphic, with separate panels for Bhat and gamma errors
  boxplots <- ggplot(all_errors, aes(x = method, y = value, fill = method)) +
    geom_boxplot() +
    facet_wrap(~error_type, scales = "free_y") +
    labs(title = "Comparison of Bhat and Gamma Errors Across Methods",
         y = "Error Value",
         x = "Method") +
    theme_minimal() +
    theme(legend.position = "none")
  
  ggsave(
    filename = file.path(output_dir, paste0("Boxplots-",model,".png")),
    plot = boxplots,
    width = 8, height = 6
  )
    }
  }
}
