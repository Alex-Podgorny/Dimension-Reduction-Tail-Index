#' Estimation of the base of the reduction subspace with TIREX1 method (Algbalou et al. 2024)
#'
#' @param X Matrix of covariates of dimension p.
#' @param y Response variable vector.
#' @param X0 Indicator vector for points X in compact subset X0.
#' @param q Dimension of the reduction subspace.
#' @param interm_lvl Proportion of the largest y values used for estimation.
#' @return Estimated matrix `Bhat` representing the base of the reduction subspace.
#' @export

TIREX1 <- function(X, y, X0, q, interm_lvl) {
  
  # Check if X0 is a vector of valid indices
  if (!is.vector(X0) || !all(X0 == as.integer(X0)) || any(X0 < 1) || any(X0 > nrow(X))) {
    stop("Error: 'X0' must be a vector of valid integer indices for rows in 'X'.")
  }
  # Check if 'q' is a positive integer less than or equal to ncol(X)
  if (!is.numeric(q) || q <= 0 || q > ncol(X)) {
    stop("Error: 'q' must be a positive integer less than or equal to the number of columns in 'X'.")
  }
  # Check if interm_lvl is a proportion between 0 and 1
  if (!is.numeric(interm_lvl) || interm_lvl <= 0 || interm_lvl > 1) {
    stop("Error: 'interm_lvl' must be a numeric value between 0 and 1.")
  }
  
  # Number of observations to use based on interm_lvl
  k <- floor(interm_lvl * nrow(X[X0,]))
  
  # Center and scale the selected subset of X
  Z <- scale(X[X0,])
  
  # Sort response variable 'y' and reorder Z accordingly
  res_ord <- sort(y[X0], index.return = TRUE)
  Z_ord <- Z[res_ord$ix,]
  
  # Initialize M_hat_i for the cumulative covariance matrix
  M_hat_i <- 0
  n <- nrow(Z_ord)
  
  # Calculation of formula (6.1) of Algbalou et al. 2024
  for (j in 1:k) {
    S_j <- apply(matrix(Z_ord[(n - j + 1):n,], nrow = j), 2, sum)
    M_hat_i <- M_hat_i + S_j %*% t(S_j)
  }
  
  # Compute eigen decomposition of M_hat_i
  res_Meigen <- eigen(M_hat_i)
  
  # Extract the top 'q' eigenvectors to form Bhat
  Bhat <- matrix(res_Meigen$vectors[, 1:q], ncol = q)
  
  return(Bhat)
}

#' Estimation of the base of the reduction subspace with TIREX2 method (Algbalou et al. 2024)
#'
#' @param X Matrix of covariates of dimension p.
#' @param y Response variable vector.
#' @param X0 Indicator vector for points X in compact subset X0.
#' @param q Dimension of the reduction subspace.
#' @param interm_lvl Proportion of the largest y values used for estimation.
#' @return Estimated matrix `Bhat` representing the base of the reduction subspace.
#' @export

TIREX2 <- function(X, y, X0, q, interm_lvl) {
  
  # Check if X0 is a vector of valid indices
  if (!is.vector(X0) || !all(X0 == as.integer(X0)) || any(X0 < 1) || any(X0 > nrow(X))) {
    stop("Error: 'X0' must be a vector of valid integer indices for rows in 'X'.")
  }
  # Check if 'q' is a positive integer less than or equal to ncol(X)
  if (!is.numeric(q) || q <= 0 || q > ncol(X)) {
    stop("Error: 'q' must be a positive integer less than or equal to the number of columns in 'X'.")
  }
  # Check if interm_lvl is a proportion between 0 and 1
  if (!is.numeric(interm_lvl) || interm_lvl <= 0 || interm_lvl > 1) {
    stop("Error: 'interm_lvl' must be a numeric value between 0 and 1.")
  }
  
  # Number of observations to use based on interm_lvl
  k <- floor(interm_lvl * nrow(X[X0,]))
  
  # Center and scale the selected subset of X
  Z <- scale(X[X0,])
  
  # Sort response variable 'y' and reorder Z accordingly
  res_ord <- sort(y[X0], index.return = TRUE)
  Z_ord <- Z[res_ord$ix,]
  
  # Initialize M_hat_i for the cumulative covariance-like matrix
  M_hat_i <- 0
  n <- nrow(Z_ord)
  
  # Calculation of formula (6.2) of Algbalou et al. 2024
  for (j in 1:k) {
    T_j <- 0
    for (i in 1:j) {
      T_j <- T_j + t(matrix(Z_ord[(n - i + 1),], nrow = 1)) %*% matrix(Z_ord[(n - i + 1),], nrow = 1)
    }
    M_hat_i <- M_hat_i + t(T_j - j * diag(1, ncol(X))) %*% (T_j - j * diag(1, ncol(X)))
  }
  
  # Compute eigen decomposition of M_hat_i
  res_Meigen <- eigen(M_hat_i)
  
  # Extract the top 'q' eigenvectors to form Bhat
  Bhat <- matrix(res_Meigen$vectors[, 1:q], ncol = q)
  
  return(Bhat)
}