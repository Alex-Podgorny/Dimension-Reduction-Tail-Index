# Loading required package for random sequence generation
library(randtoolbox)

# Loading all R scripts from the 'Methods' directory
# These scripts define various estimation methods used in the analysis
methods_func <- list.files(path = "Methods", pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
sapply(methods_func, source)

# Set the random seed for reproducibility
# (the seed could also be passed from the command line)
# seed = as.numeric(commandArgs(trailingOnly = TRUE))
seed <- 1
set.seed(seed)

# Create a directory for saving error results, organized by seed value
output_dir <- "Simulations/Dimension known/Fixed_params/Errors/seed_"
dir.create(paste(output_dir, seed, sep=""))


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
  alpha <- n^(-0.3)
  b <- 0.2
  h <- n^(-b/q) / 2
  n0 <- ceiling(n*h^q*alpha*log(n))
  
        
  # CTI method
  Bhat_CTI <- CTI(X, y, X0[1:n0], q, alpha, h)
  Error_Bhat_CTI <- norm(Bhat_CTI %*% t(Bhat_CTI) - B_0 %*% t(B_0), "2")
  gamma_hat <- local_Hill(X, y, Grid_X0, Bhat_CTI, alpha, h)
  Error_gamma_CTI <- mean((gamma_hat - apply(t(t(B_0) %*% t(Grid_X0)), 1, xi))^2,na.rm=TRUE)
    
  Errors <- list(Bhat = Error_Bhat_CTI, gamma = Error_gamma_CTI)

  # Save the error results to a file
  save(Errors, file = paste0(output_dir, seed, "/", nameFile))
  
}
