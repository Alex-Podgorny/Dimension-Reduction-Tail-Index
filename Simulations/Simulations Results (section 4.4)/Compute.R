# Loading required package
library(randtoolbox)
library(pracma)

# Loading all R scripts from the 'Methods' directory
# These scripts define various estimation methods used in the analysis
methods_func <- list.files(path = "Methods", pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
sapply(methods_func, source)

# Set the random seed for reproducibility
# (the seed could also be passed from the command line)
#seed = as.numeric(commandArgs(trailingOnly = TRUE))
seed <- 1 
set.seed(seed)

# Create a directory for saving error results, organized by seed value
output <- "Simulations/Simulations Results (section 4.4)"

# Define hyperparameter grids for alpha and h exponents
alpha_exponants <- c(0.25,0.3,0.35,0.4,0.45,0.5)
h_exponants <- c(0,0.05,0.1,0.15,0.2,0.3,0.4)

# Loop over each dataset in the generated data folder for the given seed
for (data_name in list.files(path = paste("Simulations/Simulated Data (section 4.3)/Generated data/seed_", seed, sep="")) ) {
  
  # Define file name for saving error results associated with the current dataset
  nameFile <- paste("Errors_", data_name, sep="")
  
  # Load the generated dataset (contains X, y, and model characteristics)
  load(paste("Simulations/Simulated Data (section 4.3)/Generated data/seed_", seed, "/", data_name, sep=""))
  
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

  
  # Initialize matrices to store error results for each method
  Matrix_errors_Bhat_CTI <- matrix(NA, nrow = length(alpha_exposants), ncol = length(h_exposants))
  Matrix_errors_Bhat_TDR <- matrix(NA, nrow = length(alpha_exposants), ncol = length(h_exposants))
  Matrix_errors_Bhat_T1 <- matrix(NA, nrow = length(alpha_exposants), ncol = length(h_exposants))
  Matrix_errors_Bhat_T2 <- matrix(NA, nrow = length(alpha_exposants), ncol = length(h_exposants))
  
  Matrix_errors_gamma_CTI <- matrix(NA, nrow = length(alpha_exposants), ncol = length(h_exposants))
  Matrix_errors_gamma_TDR <- matrix(NA, nrow = length(alpha_exposants), ncol = length(h_exposants))
  Matrix_errors_gamma_T1 <- matrix(NA, nrow = length(alpha_exposants), ncol = length(h_exposants))
  Matrix_errors_gamma_T2 <- matrix(NA, nrow = length(alpha_exposants), ncol = length(h_exposants))
  Matrix_errors_gamma_B0 <- matrix(NA, nrow = length(alpha_exposants), ncol = length(h_exposants))
  Matrix_errors_gamma_Id <- matrix(NA, nrow = length(alpha_exposants), ncol = length(h_exposants))
  
  # Case for p = 4: Multiple values of hyperparameters h and alpha
  if (p == 4) {
    
    # Loop over each combination of alpha and h exponents
    for (i in 1:length(alpha_exposants)) {
      for (j in 1:length(h_exposants)) {
        
        # Set values for hyperparameters alpha and h based on sample size n
        a <- alpha_exposants[i]
        b <- h_exposants[j]
        cat("a =", a, "b =", b)
              
        alpha <- n^(-a)
        h <- n^(-b/q) / 2
        n0 <- n/20
        
        # Calculate estimation errors for each method and store them in respective matrices
        
        # CTI method
        Bhat_CTI <- CTI(X, y, intersect(X0,1:n0), q, alpha, h)
        Matrix_errors_Bhat_CTI[i, j] <- norm(Bhat_CTI - B_0, "F")
        gamma_hat <- local_Hill(X, y, Grid_X0, Bhat_CTI, alpha, h)
        Matrix_errors_gamma_CTI[i, j] <- mean((gamma_hat - apply(t(t(B_0) %*% t(Grid_X0)), 1, xi))^2,na.rm=TRUE)
        
        # Gardes method
        Bhat_TDR <- Gardes(X, y, 1:100, q, alpha, h)
        Matrix_errors_Bhat_TDR[i, j] <- norm(Bhat_TDR - B_0, "F")
        gamma_hat <- local_Hill(X, y, Grid_X0, Bhat_TDR, alpha, h)
        Matrix_errors_gamma_TDR[i, j] <- mean((gamma_hat - apply(t(t(B_0) %*% t(Grid_X0)), 1, xi))^2,na.rm=TRUE)
        
        # TIREX1 method
        Bhat_T1 <- TIREX1(X, y, 1:n, q, alpha)
        Matrix_errors_Bhat_T1[i, j] <- norm(Normalize(Bhat_T1,q) - B_0, "F")
        gamma_hat <- local_Hill(X, y, Grid_X0, Normalize(Bhat_T1,q), alpha, h)
        Matrix_errors_gamma_T1[i, j] <- mean((gamma_hat - apply(t(t(B_0) %*% t(Grid_X0)), 1, xi))^2,na.rm=TRUE)
        
        # TIREX2 method
        Bhat_T2 <- TIREX2(X, y, 1:n, q, alpha)
        Matrix_errors_Bhat_T2[i, j] <- norm(Normalize(Bhat_T2,q) - B_0, "F")
        gamma_hat <- local_Hill(X, y, Grid_X0, Normalize(Bhat_T2,q), alpha, h)
        Matrix_errors_gamma_T2[i, j] <- mean((gamma_hat - apply(t(t(B_0) %*% t(Grid_X0)), 1, xi))^2,na.rm=TRUE)
        
        # Baseline errors for B0 
        Matrix_errors_gamma_B0[i,j] <- mean((local_Hill(X, y,Grid_X0, B_0, alpha, h) - apply(t(t(B_0) %*% t(Grid_X0)), 1, xi))^2,na.rm=TRUE)
        Matrix_errors_gamma_Id[i,j] <- mean((local_Hill(X, y,Grid_X0, diag(1, p), alpha, n^(-b/p) / 2) - apply(t(t(B_0) %*% t(Grid_X0)), 1, xi))^2,na.rm=TRUE)
        
      }
    }
    
    # Save all error matrices for each method in a list
    Errors <- list(
      CTI = list(Bhat = Matrix_errors_Bhat_CTI, gamma = Matrix_errors_gamma_CTI),
      Gardes = list(Bhat = Matrix_errors_Bhat_TDR, gamma = Matrix_errors_gamma_TDR),
      TIREX1 = list(Bhat = Matrix_errors_Bhat_T1, gamma = Matrix_errors_gamma_T1),
      TIREX2 = list(Bhat = Matrix_errors_Bhat_T2, gamma = Matrix_errors_gamma_T2),
      B0 = list(gamma = Matrix_errors_gamma_B0),
      Id = list(gamma = Matrix_errors_gamma_Id)
    )
    
    # Save the error results to a file
    save(Errors, file = paste0(output, seed, "/", nameFile))
    
  }
  
  # Case for p = 30: Single value for hyperparameters alpha and h (Gardes's method is not computable)
  if (p == 30) {
    
    # Loop over each combination of alpha and h exponents
    for (i in 1:length(alpha_exposants)) {
      for (j in 1:length(h_exposants)) {
        
        # Set values for hyperparameters alpha and h based on sample size n
        a <- alpha_exposants[i]
        b <- h_exposants[j]
        cat("a =", a, "b =", b)
        
        alpha <- n^(-a)
        h <- n^(-b/q) / 2
        n0 <- n/20
        
        # Calculate estimation errors for each method and store them in respective matrices
        
        # CTI method
        Bhat_CTI <- CTI(X, y, intersect(X0,1:n0), q, alpha, h)
        Matrix_errors_Bhat_CTI[i, j] <- norm(Bhat_CTI - B_0, "F")
        gamma_hat <- local_Hill(X, y, Grid_X0, Bhat_CTI, alpha, h)
        Matrix_errors_gamma_CTI[i, j] <- mean((gamma_hat - apply(t(t(B_0) %*% t(Grid_X0)), 1, xi))^2,na.rm=TRUE)
        
        # TIREX1 method
        Bhat_T1 <- TIREX1(X, y, 1:n, q, alpha)
        Matrix_errors_Bhat_T1[i, j] <- norm(Normalize(Bhat_T1,q) - B_0, "F")
        gamma_hat <- local_Hill(X, y, Grid_X0, Normalize(Bhat_T1,q), alpha, h)
        Matrix_errors_gamma_T1[i, j] <- mean((gamma_hat - apply(t(t(B_0) %*% t(Grid_X0)), 1, xi))^2,na.rm=TRUE)
        
        # TIREX2 method
        Bhat_T2 <- TIREX2(X, y, 1:n, q, alpha)
        Matrix_errors_Bhat_T2[i, j] <- norm(Normalize(Bhat_T2,q) - B_0, "F")
        gamma_hat <- local_Hill(X, y, Grid_X0, Normalize(Bhat_T2,q), alpha, h)
        Matrix_errors_gamma_T2[i, j] <- mean((gamma_hat - apply(t(t(B_0) %*% t(Grid_X0)), 1, xi))^2,na.rm=TRUE)
        
        # Baseline errors for B0 
        Matrix_errors_gamma_B0[i,j] <- mean((local_Hill(X, y,Grid_X0, B_0, alpha, h) - apply(t(t(B_0) %*% t(Grid_X0)), 1, xi))^2,na.rm=TRUE)
        Matrix_errors_gamma_Id[i,j] <- mean((local_Hill(X, y,Grid_X0, diag(1, p), alpha, n^(-b/p) / 2) - apply(t(t(B_0) %*% t(Grid_X0)), 1, xi))^2,na.rm=TRUE)
      }
    }
    
    # Save all error matrices for each method in a list
    Errors <- list(
      CTI = list(Bhat = Matrix_errors_Bhat_CTI, gamma = Matrix_errors_gamma_CTI),
      TIREX1 = list(Bhat = Matrix_errors_Bhat_T1, gamma = Matrix_errors_gamma_T1),
      TIREX2 = list(Bhat = Matrix_errors_Bhat_T2, gamma = Matrix_errors_gamma_T2),
      B0 = list(gamma = Matrix_errors_gamma_B0),
      Id = list(gamma = Matrix_errors_gamma_Id)
    )
    
    # Save the error results to a file
    save(Errors, file = paste0(output, seed, "/", nameFile))
    
  }
  
}

