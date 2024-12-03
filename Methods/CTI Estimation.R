#' Estimated Psi function (objective function to minimize)
#'
#' @param X Matrix of covariates of dimension p.
#' @param y Response variable vector.
#' @param N0 Indicator vector specifying a sample within the compact subset \( X_0 \).
#' @param B (p*q)-Matrix for dimension reduction from p to q.
#' @param interm_lvl Proportion of the largest y values used for estimation.
#' @param bandwidth Distance to select the observations closest to the point of interest z.
#' @return Mean of estimated tail-index 
#' @export

Psi_function <- function(X, y, N0, B, interm_lvl, bandwidth, mink) {
  
  # Check that X0 is a vector of valid indices
  if (!is.vector(N0) || !all(N0 == as.integer(N0))) stop("Error: 'N0' must be a vector of integer indices.")
  if (any(N0 < 1) || any(N0 > nrow(X))) stop("Error: 'N0' contains indices that are out of bounds for 'X'.")
  
  # Calculate the mean of the conditional tail-index on X0
  mean(local_Hill(X, y, X[N0, ], B, interm_lvl, bandwidth,mink), na.rm = TRUE)
}



#' Estimation of the Base of the Central Tail-Index (CTI) Subspace
#'
#' This function estimates the base of the central tail-index (CTI) subspace, 
#' which is a low-dimensional subspace of interest in extreme value analysis. 
#' The CTI subspace is identified using a objective function `Psi` based on 
#' the largest values of the response variable `y`.
#'
#' @param X A matrix of covariates with `p` columns (dimensions).
#' @param y A vector of response variables. The largest values of `y` are used
#'   for estimation based on the `interm_lvl` parameter.
#' @param N0 Indicator vector specifying a sample within the compact subset \( X_0 \).
#' @param q Dimension of the CTI subspace to be estimated.
#' @param interm_lvl Proportion of the largest values of `y` to consider in 
#'   the estimation. Must be a value between 0 and 1.
#' @param bandwidth Bandwidth parameter that determines the distance used to
#'   select observations closest to the point of interest.
#' @param control A list of control parameters for the optimization process:
#'   \describe{
#'     \item{`num_init_evals`}{Number of initial random evaluations to generate starting points.}
#'     \item{`num_starts`}{Number of starting points for the optimization.}
#'     \item{`tol`}{Tolerance level for the convergence of the optimization algorithm.}
#'     \item{`mink`}{Minimum number of points should be used to estimate the tail-index.}
#'   }
#' @return The estimated matrix `Bhat`, which represents the base of the CTI subspace.
#' @export
CTI <- function(X, y, N0, q, interm_lvl, bandwidth, control = list()) {
  
  # Default control parameters
  default_control <- list(
    num_init_evals = 100,  
    num_starts = 1,        
    tol = 1e-2,            
    mink = 1               
  )
  
  control <- modifyList(default_control, control)
  
  
  # Step 1: Define the objective function to minimize
  # The objective function computes the criterion Psi_function for a given basis `B`.
  # It ensures that `B` is orthonormalized and evaluates its performance on the data.
  objective_fn <- function(B) {
    B <- orthonormalize(B, q)  # Ensure the basis matrix `B` is orthonormal
    Psi_function(X, y, N0, B, interm_lvl, bandwidth, mink = control$mink)  # Compute the criterion
  }
  
  # Step 2: Estimate the base of the CTI subspace using the Minimization function
  Bhat <- Minimization(
    objective_fn,                   # Objective function to minimize
    c(ncol(X), q),                  # Dimensions of the matrix to optimize (p x q)
    control$num_init_evals,         # Number of initial evaluations
    control$num_starts,             # Number of optimization starting points
    control$tol                     # Tolerance for convergence
  )
  
  # Step 3: Return the estimated basis matrix
  return(Bhat)
}
