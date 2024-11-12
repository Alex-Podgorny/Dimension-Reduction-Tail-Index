# Load necessary packages
library(ggplot2)
library(tidyr)

# Define the directory containing the results
output_dir <- "Simulations/Dimension known/Fixed_params/Results"

# Define the base directory containing the seed_n folders
base_dir <- "Simulations/Dimension known/Fixed_params/Errors/"


# Retrieve all directories
models <- list("model1p4","model2p4","model3p4","model4p4",
               "model1p30","model2p30","model3p30","model4p30")

seed_dirs <- list.files(base_dir, pattern = "^seed_")


# Initialize a list of lists to store the error matrices by method
methods <- c("CTI", "Gardes", "TIREX1", "TIREX2", "B0", "Id")
errors_Bhat <- list()
errors_gamma <- list()

for(model in models){
  
  # Loop through each seed directory to load errors
  for (seed_dir in seed_dirs) {
    # Load each error file for the current seed
    error_file <- paste0(base_dir,seed_dir,"/Errors_",model,"-",seed_dir,".RData")
    
    load(error_file)  # Load the `Errors` list
    
    # Loop through each method to store Bhat error matrices
    for (method in methods) {
      if (!is.null(Errors[[method]]$Bhat)) {
        # Append the current Bhat error matrix to the list for the corresponding method
        errors_Bhat[[method]] <- c(errors_Bhat[[method]],Errors[[method]]$Bhat)
      }
      if (!is.null(Errors[[method]]$gamma)) {
        # Append the current Bhat error matrix to the list for the corresponding method
        errors_gamma[[method]] <- c(errors_gamma[[method]],Errors[[method]]$gamma) 
      }
    }
  }
  
  # Create data frames for Bhat and gamma
  all_errors_Bhat <- data.frame()
  all_errors_gamma <- data.frame()
  
  for (method in methods) {
    
    if(!is.null(Errors[[method]]$Bhat)){
      all_errors_Bhat <- rbind(all_errors_Bhat, data.frame(value = errors_Bhat[[method]], method = method, error_type = "Bhat"))
    }
    if(!is.null(Errors[[method]]$gamma)){
      all_errors_gamma <- rbind(all_errors_gamma, data.frame(value = errors_gamma[[method]], method = method, error_type = "Gamma"))
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
    filename = file.path(output_dir, paste0("Boxplots-",model,".png")),
    plot = boxplots,
    width = 8, height = 6
  )  
  
} 

