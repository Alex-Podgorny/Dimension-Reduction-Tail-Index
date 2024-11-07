#' Estimated Psi function (objective function to minimize)
#'
#' @param X Matrix of covariates of dimension p.
#' @param y Response variable vector.
#' @param X0 Indicator vector for points X in compact subset X0.
#' @param B (p*q)-Matrix for dimension reduction from p to q.
#' @param interm_lvl Proportion of the largest y values used for estimation.
#' @param bandwith Distance to select the observations closest to the point of interest z.
#' @return Mean of estimated tail-index 
#' @export

Psi_function <- function(X, y, X0, B, interm_lvl, bandwith) {
  
  # Check that X0 is a vector of valid indices
  if (!is.vector(X0) || !all(X0 == as.integer(X0))) stop("Error: 'X0' must be a vector of integer indices.")
  if (any(X0 < 1) || any(X0 > nrow(X))) stop("Error: 'X0' contains indices that are out of bounds for 'X'.")
  
  # Calculate the mean of the conditional tail-index on X0
  mean(local_Hill(X, y, X[X0, ], B, interm_lvl, bandwith), na.rm = TRUE)
}



#' Estimation of the base of the central tail-index (CTI) subspace
#'
#' @param X Matrix of covariates of dimension p.
#' @param y Response variable vector.
#' @param X0 Indicator vector for points X in compact subset X0.
#' @param q Dimension of the CTI subspace.
#' @param interm_lvl Proportion of the largest y values used for estimation.
#' @param bandwith Distance to select the observations closest to the point of interest z.
#' @param num_init_evals Number of initial random evaluations for starting points.
#' @param num_starts Number of starting points to use for local optimization.
#' @return Estimated matrix `Bhat` representing the base of the CTI subspace.
#' @export

CTI <- function(X, y, X0, q, interm_lvl, bandwith, num_init_evals = 1000, num_starts = 3) {
  
  if (!requireNamespace("randtoolbox", quietly = TRUE)) {
    stop("Package 'randtoolbox' is required but not installed.")
  }
  
  # Function to compute the orthonormal matrix from a flat vector using QR decomposition
  orthonormalize <- function(A) {
    qr.Q(qr(matrix(A, nrow = p, ncol = q)))
  }
  
  # Define the objective function with Gram-Schmidt orthonormalization
  objective_fn <- function(B) {
    B <- orthonormalize(B)
    Psi_function(X, y, X0, B, interm_lvl, bandwith)
  }
  
  p <- ncol(X)
  
  # Generate Halton sequence for initial evaluations and reshape as matrices
  init_evals <- randtoolbox::halton(num_init_evals, p * q)
  
  # Reshape each row of init_evals to a p*q matrix, apply orthonormalization, and flatten back to a vector
  Inits <- t(apply(init_evals, 1, function(row) c(orthonormalize(row))))
  
  # Compute objective function values for each initialized matrix
  Value_inits <- apply(Inits, 1, objective_fn)
  
  # Initialize lists to store optimized matrices and their objective values
  Bhat_s <- vector("list", num_starts)
  Value_s <- numeric(num_starts)
  
  # Select the indices of the starting points with the lowest objective values
  best_starts <- order(Value_inits)[1:num_starts]
  
  # Optimize locally for each selected starting point
  for (s in seq_len(num_starts)) {
    # Perform local optimization using `optim` starting from each chosen initial point
    result <- optim(
      par = Inits[best_starts[s], ], 
      fn = objective_fn, 
      control = list(trace = 1, maxit = 100, reltol = 1e-3)
    )
    
    # Store the normalized optimized matrix and its objective value
    Bhat_s[[s]] <- orthonormalize(result$par)
    Value_s[s] <- result$value
  }
  
  # Select the optimized matrix with the lowest objective value
  Bhat <- Bhat_s[[which.min(Value_s)]]
  return(Bhat)
}