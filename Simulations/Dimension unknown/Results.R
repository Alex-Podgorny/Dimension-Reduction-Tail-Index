
# Define the directory containing the results
output_dir <- "Simulations/Dimension unknown"

# Define the base directory containing the seed_n folders
base_dir <- "Simulations/Dimension unknown/Errors/"


# Retrieve all directories
models <- list("model1p4","model2p4","model3p4","model4p4",
               "model1p30","model2p30","model3p30","model4p30")


seed_dirs <- list.files(base_dir, pattern = "^seed_")


# Initialize matrix results
Matrix_Choices <- matrix(NA,nrow = 8, ncol=3)
Matrix_Errors <- matrix(NA,nrow=8,ncol = 4)
i <- 0

for(model in models){
  
  i <- i + 1
  
  # Initialize vectors to store the error matrices by method
  errors_1 <- c()
  errors_2 <- c()
  errors_3 <- c()
  errors_qhat <- c()
  qhat <- c()
  
  # Loop through each seed directory to load errors
  for (seed_dir in seed_dirs) {
    # Load each error file for the current seed
    error_file <- paste0(base_dir,seed_dir,"/Errors_",model,"-",seed_dir,".RData")
    
    load(error_file)  # Load the `Errors` list
    
    errors_1 <- c(errors_1,Errors$Error_q[1])
    errors_2 <- c(errors_2,Errors$Error_q[2])
    errors_3 <- c(errors_3,Errors$Error_q[3])
    errors_qhat <- c(errors_qhat,Errors$Error_q[Errors$q_hat])
    
    qhat <- c(qhat,Errors$q_hat)
    
  }
  
  Matrix_Choices[i,] <- c(mean(qhat==1),mean(qhat==2),mean(qhat==3))
  
  Matrix_Errors[i,] <- c(mean(errors_1),
                         mean(errors_2),
                         mean(errors_3),
                         mean(errors_qhat))
} 

Results <- list(Matrix_Choices = Matrix_Choices, Matrix_Errors = Matrix_Errors)

save(Results, file = file.path(output_dir, "Resuts.RData"))
