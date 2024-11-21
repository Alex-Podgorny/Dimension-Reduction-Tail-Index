# Load necessary packages
library(ggplot2)


# Define the directory containing the results
output_dir <- "Simulations/Dimension known/Best_params/Results"

# Define the base directory containing the seed_n folders
base_dir <- "Simulations/Dimension known/Best_params/Errors/"


# Retrieve all directories
models <- list("model1p4","model2p4","model3p4","model4p4",
               "model1p30","model2p30","model3p30","model4p30")


seed_dirs <- list.files(base_dir, pattern = "^seed_")



# Initialize a list of lists to store the error matrices by method
methods <- c("CTI", "Gardes", "TIREX1", "TIREX2", "B0", "Id")

for(model in models){
  
  errors_Bhat <- list()
  errors_gamma <- list()
  
  # Loop through each seed directory to load errors
  for (seed_dir in seed_dirs) {
    # Load each error file for the current seed
    error_file <- paste0(base_dir,seed_dir,"/Errors_",model,"-",seed_dir,".RData")
    
    load(error_file)  # Load the `Errors` list

    # Loop through each method to store Bhat error matrices
    for (method in methods) {
        if (!is.null(Errors[[method]]$Bhat)) {
          # Append the current Bhat error matrix to the list for the corresponding method
          errors_Bhat[[method]] <- append(errors_Bhat[[method]], list(Errors[[method]]$Bhat))
        }
        if (!is.null(Errors[[method]]$gamma)) {
          # Append the current Bhat error matrix to the list for the corresponding method
          errors_gamma[[method]] <- append(errors_gamma[[method]], list(Errors[[method]]$gamma))
        }
      }
  }
      
  # Calculate the mean matrix for each method across all seeds
  means_Bhat <- lapply(errors_Bhat, function(method_errors) {
        
        if (length(method_errors) > 1) {
          mean_matrix <- Reduce("+", method_errors) / length(method_errors)
        } else {
          mean_matrix <- method_errors[[1]]
        }
        mean_matrix
      }) 
      
      
   means_gamma <- lapply(errors_gamma, function(method_errors) {
        
        if (length(method_errors) > 1) {
          mean_matrix <- Reduce("+", method_errors) / length(method_errors)
        } else {
          mean_matrix <- method_errors[[1]]
        }
        mean_matrix
      })  
      
      
    # Create heatmaps for each method 
      
    for (method in methods) {
        # Convert the mean error matrix to a data frame for ggplot2
        if (!is.null(means_Bhat[[method]])) {
          mean_Bhat_df_long <- data.frame(
            Var1 = rep(seq_len(nrow(means_Bhat[[method]])), ncol(means_Bhat[[method]])),
            Var2 = rep(seq_len(ncol(means_Bhat[[method]])), each = nrow(means_Bhat[[method]])),
            value = as.vector(means_Bhat[[method]])
          )
          
          # Heatmap for Bhat error
          plot_Bhat <- ggplot(mean_Bhat_df_long, aes(Var1, Var2, fill = value)) +
            geom_tile(color = "white") +
            geom_text(aes(label = round(value, 3)), color = "black") + 
            scale_fill_gradient2(  low = "yellow",
                                   mid = "white",
                                   high = "purple",
                                   midpoint = median(mean_Bhat_df_long$value),
                                   transform = "log10") +
            scale_x_continuous(breaks = 1:6, labels = c(0.2,0.25,0.3,0.35,0.4,0.45)) +
            scale_y_continuous(breaks = 1:6, labels = c(0.01,0.05,0.1,0.15,0.2,0.3)) +
            labs(title = paste("Mean Error Bhat for method", method), x = "alpha exponents", y = "h exponents") +
            theme_minimal()
          
          ggsave(
            filename = file.path(output_dir, paste0("Mean_Error_Bhat_",model,"-", method, ".png")),
            plot = plot_Bhat,
            width = 8, height = 6
          )
        }
      
      if (!is.null(means_gamma[[method]])){
        mean_gamma_df_long <- data.frame(
          Var1 = rep(seq_len(nrow(means_gamma[[method]])), ncol(means_gamma[[method]])),
          Var2 = rep(seq_len(ncol(means_gamma[[method]])), each = nrow(means_gamma[[method]])),
          value = as.vector(means_gamma[[method]])
        )
        
        # Heatmap for gamma error
        plot_gamma <- ggplot(mean_gamma_df_long, aes(Var1, Var2, fill = value)) +
          geom_tile(color = "white") +
          geom_text(aes(label = round(value, 3)), color = "black") +
          scale_fill_gradient2(  low = "yellow",
                                   mid = "white",
                                   high = "purple",
                                   midpoint = median(mean_gamma_df_long$value),
                                   transform = "log10") +
          scale_x_continuous(breaks = 1:6, labels =  c(0.2,0.25,0.3,0.35,0.4,0.45)) +
          scale_y_continuous(breaks = 1:6, labels = c(0.01,0.05,0.1,0.15,0.2,0.3)) +
          labs(title = paste("Mean Error Gamma for method", method), x = "alpha exponents", y = "h exponents") +
          theme_minimal()
        
        
        ggsave(
          filename = file.path(output_dir, paste0("Mean_Error_Gamma_",model,"-", method, ".png")),
          plot = plot_gamma,
          width = 8, height = 6
        )
      }
    }
      
      
      
      
      # Find the optimal hyperparameters for each method (p = 4)
      min_params_gamma  <- lapply(means_gamma , function(matrix) {
        which(matrix == min(matrix), arr.ind = TRUE)[1,]
      })
      
      # Create data frames for Bhat and gamma errors with optimal hyperparameters
      all_errors_Bhat <- data.frame()
      all_errors_gamma <- data.frame()
      
      for (method in methods) {
        if (!is.null(means_Bhat[[method]])){
          # Get the optimal alpha and h indices for Bhat and gamma
          alpha_idx <- min_params_gamma[[method]][1]
          h_idx <- min_params_gamma[[method]][2]
          
          # Extract Bhat errors for optimal parameters and add to the combined data frame
          errors_Bhat_best <- unlist(lapply(errors_Bhat[[method]], function(matrix) matrix[alpha_idx, h_idx]))
          all_errors_Bhat <- rbind(all_errors_Bhat, data.frame(value = errors_Bhat_best, method = method, error_type = "Bhat"))
        }
        if (!is.null(means_gamma[[method]])){
          # Get the optimal alpha and h indices for Bhat and gamma
          alpha_idx <- min_params_gamma[[method]][1]
          h_idx <- min_params_gamma[[method]][2]
          
          # Extract gamma errors for optimal parameters and add to the combined data frame
          errors_gamma_best <- unlist(lapply(errors_gamma[[method]], function(matrix) matrix[alpha_idx, h_idx]))
          all_errors_gamma <- rbind(all_errors_gamma, data.frame(value = errors_gamma_best, method = method, error_type = "Gamma"))
        }
      }
      # Combine Bhat and gamma errors into a single data frame
      all_errors <- rbind(all_errors_Bhat, all_errors_gamma)
      
      save(all_errors,file=file.path(output_dir, paste0("Errors-",model,".RData")))
      
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
        filename = file.path(output_dir, paste0("Boxplots-",model,".pdf")),
        plot = boxplots,
        width = 8, height = 6
      )  
  } 



