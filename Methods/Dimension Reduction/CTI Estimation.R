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
  
  # Define the objective function
  objective_fn <- function(B) {
    B <- qr.Q(qr(matrix(B, nrow = p, ncol = q)))
    Psi_function(X, y, X0, B, interm_lvl, bandwith)
  }
  
  # Minimization of the objective function
  Bhat <- Minimization(objective_fn,c(ncol(X),q),num_init_evals,num_starts)
  
  return(Bhat)
}