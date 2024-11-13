
# Define the directory containing the results
output_dir <- "Simulations/Dimension unknown"

# Define the base directory containing the seed_n folders
base_dir <- "Simulations/Dimension unknown/Errors/"


# Retrieve all directories
models <- list("model1p4","model2p4","model3p4","model4p4",
               "model1p30","model2p30","model3p30","model4p30")


seed_dirs <- list.files(base_dir, pattern = "^seed_")



# Initialize vectors to store the error matrices by method
errors_1 <- c()
errors_2 <- c()
errors_3 <- c()
errors_qhat <- c()
q_hat <- c()

# Initialize matrix results
Choice <- matrix(NA,nrow = 8, ncol=3)
Errors <- matrix(NA,nrow=8,ncol = 8)
i <- 0

for(model in models){
  
  i <- i + 1
  
  # Loop through each seed directory to load errors
  for (seed_dir in seed_dirs) {
    # Load each error file for the current seed
    error_file <- paste0(base_dir,seed_dir,"/Errors_",model,"-",seed_dir,".RData")
    
    load(error_file)  # Load the `Errors` list
    
    errors_1 <- c(errors_1,Errors$Error_q[1])
    errors_2 <- c(errors_2,Errors$Error_q[2])
    errors_3 <- c(errors_3,Errors$Error_q[3])
    errors_qhat <- c(errors_qhat,Errors$Error_q[Errors$q_hat])
    
    q_hat <- c(q_hat,Errors$q_hat)
    
  }
  
  Choice[i,] <- c(mean(q_hat==1),mean(q_hat==2),mean(q_hat==3))
  
  Errors[i,] <- c(mean(errors_1),sd(errors_1),
                   mean(errors_2),sd(errors_2),
                   mean(errors_3),sd(errors_3),
                   mean(errors_qhat),sd(errors_qhat)
                    )
} 

Results <- list(Choice = Choice, Errors = Errors)

save(Results, filename = file.path(output_dir, "Resuts.RData"))
