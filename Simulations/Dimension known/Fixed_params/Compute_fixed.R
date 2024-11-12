# Loading required package for random sequence generation
library(randtoolbox)

# Loading all R scripts from the 'Methods' directory
# These scripts define various estimation methods used in the analysis
methods_func <- list.files(path = "Methods", pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
sapply(methods_func, source)

# Set the random seed for reproducibility
# (the seed could also be passed from the command line)
# seed = as.numeric(commandArgs(trailingOnly = TRUE))
seed <- 2 

# Create a directory for saving error results, organized by seed value
output_dir <- "Simulations/Dimension known/Fixed_params/Errors/seed_"
dir.create(paste(ouput_dir, seed, sep=""))


# Loop over each dataset in the generated data folder for the given seed
for (data_name in list.files(path = paste("Simulations/Generated data/seed_", seed, sep="")) ) {
  
  # Define file name for saving error results associated with the current dataset
  nameFile <- paste("Errors_", data_name, sep="")
  
  # Load the generated dataset (contains X, y, and model characteristics)
  load(paste("Simulations/Generated data/seed_", seed, "/", data_name, sep=""))
  
  # Extract data and characteristics from the loaded dataset
  X <- data$X           # Matrix of covariates
  y <- data$y           # Response variable vector
  p <- ncol(X)          # Dimension of X
  n <- length(y)        # Sample size
  q <- data$carac$q     # Dimension of the CTI subspace
  B_0 <- data$carac$B_0 # True base of the reduction subspace
  xi <- data$carac$xi   # Function for tail index
  
  # Define epsilon and identify points in compact subset X0
  epsilon <- (1 - 0.9^(1/p)) / 2
  X0 <- which(apply(X, 1, function(x) min(min(x), min(abs(x - 1)))) > epsilon)
  
  # Generate grid points for estimating tail index over X0
  Grid_X0 <- randtoolbox::halton(10000, p) * (1 - 2 * epsilon) + epsilon
  
  
  # Set fixed values for alpha and h based on the sample size
  alpha <- n^(-0.2)
  h <- n^(-0.3/q) / 2
  
  # Case for p = 4: Multiple values of hyperparameters h and alpha
  if (p == 4) {
        
    # Calculate estimation errors for each method and store them in respective matrices
        
    # # CTI method
    # Bhat_CTI <- CTI(X, y, X0, q, alpha, h)
    # Error_Bhat_CTI <- norm(Bhat_CTI %*% t(Bhat_CTI) - B_0 %*% t(B_0), "2")
    # Error_gamma_CTI <- mean((local_Hill(X, y, Grid_X0, Bhat_CTI, alpha, h) - apply(t(t(B_0) %*% t(Grid_X0)), 1, xi))^2,na.rm=TRUE))
    
    # # Gardes method
    # Bhat_TDR <- Gardes(X, y, X0, q, alpha, h)
    # Error_Bhat_TDR <- norm(Bhat_TDR %*% t(Bhat_TDR) - B_0 %*% t(B_0), "2")
    # Error_gamma_TDR <- mean((local_Hill(X, y, Grid_X0, Bhat_TDR, alpha, h) - apply(t(t(B_0) %*% t(Grid_X0)), 1, xi))^2,na.rm=TRUE))
    
    # TIREX1 method
    Bhat_T1 <- TIREX1(X, y, X0, q, alpha)
    Error_Bhat_T1 <- norm(Bhat_T1 %*% t(Bhat_T1) - B_0 %*% t(B_0), "2")
    Error_gamma_T1 <- mean((local_Hill(X, y, Grid_X0, Bhat_T1, alpha, h) - apply(t(t(B_0) %*% t(Grid_X0)), 1, xi))^2,na.rm=TRUE)
    
    # TIREX2 method
    Bhat_T2 <- TIREX2(X, y, X0, q, alpha)
    Error_Bhat_T2 <- norm(Bhat_T2 %*% t(Bhat_T2) - B_0 %*% t(B_0), "2")
    Error_gamma_T2 <- mean((local_Hill(X, y, Grid_X0, Bhat_T2, alpha, h) - apply(t(t(B_0) %*% t(Grid_X0)), 1, xi))^2,na.rm=TRUE)
    
    # Baseline errors for B0 and Identity matrix methods
    Error_gamma_B0 <- mean((local_Hill(X, y,Grid_X0, B_0, alpha, h) - apply(t(t(B_0) %*% t(Grid_X0)), 1, xi))^2,na.rm=TRUE)
    Error_gamma_Id <- mean((local_Hill(X, y,Grid_X0, diag(1, p), alpha, n^(-2/(9 * p)) / 2) - apply(t(t(B_0) %*% t(Grid_X0)), 1, xi))^2,na.rm=TRUE)
    
    # Save all error results in a list
    Errors <- list(
      #CTI = list(Bhat = Error_Bhat_CTI, gamma = Error_gamma_CTI),
      #Gardes = list(Bhat = Error_Bhat_TDR, gamma = Error_gamma_TDR),
      TIREX1 = list(Bhat = Error_Bhat_T1, gamma = Error_gamma_T1),
      TIREX2 = list(Bhat = Error_Bhat_T2, gamma = Error_gamma_T2),
      B0 = list(gamma = Error_gamma_B0),
      Id = list(gamma = Error_gamma_Id)
    )
    
    # Save the error results to a file
    save(Errors, file = paste0(output_dir, seed, "/", nameFile))
  }
  
  # Case for p = 30: Single value for hyperparameters alpha and h (Gardes's method is not computable)
  if (p == 30) {
    
    # Calculate estimation errors for each method and store them in respective matrices
    
    # # CTI method
    # Bhat_CTI <- CTI(X, y, X0, q, alpha, h)
    # Error_Bhat_CTI <- norm(Bhat_CTI %*% t(Bhat_CTI) - B_0 %*% t(B_0), "2")
    # Error_gamma_CTI <- mean((local_Hill(X, y, Grid_X0, Bhat_CTI, alpha, h) - apply(t(t(B_0) %*% t(Grid_X0)), 1, xi))^2,na.rm=TRUE))
    
    # TIREX1 method
    Bhat_T1 <- TIREX1(X, y, X0, q, alpha)
    Error_Bhat_T1 <- norm(Bhat_T1 %*% t(Bhat_T1) - B_0 %*% t(B_0), "2")
    Error_gamma_T1 <- mean((local_Hill(X, y, Grid_X0, Bhat_T1, alpha, h) - apply(t(t(B_0) %*% t(Grid_X0)), 1, xi))^2,na.rm=TRUE)
    
    # TIREX2 method
    Bhat_T2 <- TIREX2(X, y, X0, q, alpha)
    Error_Bhat_T2 <- norm(Bhat_T2 %*% t(Bhat_T2) - B_0 %*% t(B_0), "2")
    Error_gamma_T2 <- mean((local_Hill(X, y, Grid_X0, Bhat_T2, alpha, h) - apply(t(t(B_0) %*% t(Grid_X0)), 1, xi))^2,na.rm=TRUE)
    
    # Baseline errors for B0 and Identity matrix methods
    Error_gamma_B0 <- mean((local_Hill(X, y,Grid_X0, B_0, alpha, h) - apply(t(t(B_0) %*% t(Grid_X0)), 1, xi))^2,na.rm=TRUE)
    
    # Save all error results in a list
    Errors <- list(
      #CTI = list(Bhat = Error_Bhat_CTI, gamma = Error_gamma_CTI),
      TIREX1 = list(Bhat = Error_Bhat_T1, gamma = Error_gamma_T1),
      TIREX2 = list(Bhat = Error_Bhat_T2, gamma = Error_gamma_T2),
      B0 = list(gamma = Error_gamma_B0)
    )
    
    # Save the error results to a file
    save(Errors, file = paste0(output_dir, seed, "/", nameFile))
  
  }
  
}
