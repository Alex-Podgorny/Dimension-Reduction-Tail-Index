
# Define the directory containing the results
output_dir <- "Simulations/Dimension known/Fixed_params/Results/"

# Define the base directory containing the seed_n folders
base_dir <- "Simulations/Dimension known/Fixed_params/Errors/"


# Retrieve all directories
models <- list("model1p4","model2p4","model3p4","model4p4",
               "model1p30","model2p30","model3p30","model4p30")

seed_dirs <- list.files(base_dir, pattern = "^seed_")



# Initialize a list of lists to store the error matrices by method
errors_Bhat <- c()
errors_gamma <- c()

for(model in models){
  
  # Loop through each seed directory to load errors
  for (seed_dir in seed_dirs) {
    # Load each error file for the current seed
    error_file <- paste0(base_dir,seed_dir,"/Errors_",model,"-",seed_dir,".RData")
    
    load(error_file)  # Load the `Errors` list
    
    errors_Bhat[[model]] <- c(errors_Bhat[[model]],Errors$Bhat)
    errors_gamma[[model]] <- c(errors_gamma[[model]],Errors$gamma) 
  }
  
} 

save(errors_Bhat,file=paste0(output_dir, "Errors_Bhat.RData"))
save(errors_gamma,file=paste0(output_dir, paste0("Errors_gamma.RData")))


