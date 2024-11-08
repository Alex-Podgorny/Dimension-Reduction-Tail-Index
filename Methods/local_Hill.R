#' Local_Hill estimator of the conditional tail-index with dimension reduction.
#'
#' @param X Matrix of covariates of dimension p.
#' @param y Response variable vector.
#' @param Z Matrix of points of interest.
#' @param B (p*q)-Matrix for dimension reduction from p to q.
#' @param interm_lvl Proportion of the largest y values used for estimation.
#' @param bandwith Distance to select the observations closest to the point of interest z.
#' @return Vector of estimated tail-index. 
#' @export


local_Hill = function(X, y, Z, B, interm_lvl, bandwith) {
  
  # Check if inputs have correct dimensions and values
  if (ncol(X) != ncol(Z)) stop("Error: X and Z must have the same number of columns (dimension p).")
  if (nrow(B) != ncol(X)) stop("Error: B must have the same number of rows as columns in X (dimension p).")
  if (interm_lvl <= 0 || interm_lvl > 1) stop("Error: 'interm_lvl' must be between 0 and 1.")
  if (bandwith <= 0) stop("Error: 'bandwith' must be positive.")
  
  
  # Product storage for B' * X
  B.X <- lapply(1:ncol(B), function(d) c(B[, d] %*% t(X)))
  
  # Initialize vector for estimated tail indexes
  Estimated_Tail_index <- rep(NA, nrow(Z))
  
  for (j in 1:nrow(Z)) {
    # Point of interest z
    z = t(Z[j, ])
    
    # Selection of the observations closest to the point of interest z
    Indicator = which(abs(B.X[[1]] - as.numeric(B[, 1] %*% t(z))) < bandwith)
    if (ncol(B) != 1) {
      for (d in 2:ncol(B)) {
        Indicator = intersect(Indicator, which(abs(B.X[[d]] - as.numeric(B[, d] %*% t(z))) < bandwith))
      }
    }
    
    # Select y values for the current point of interest
    y_z = y[Indicator]
    
    # Check if y_z is empty
    if (length(y_z) == 0) {
      warning(paste("Warning: No observations selected for point of interest z =", j, ". Setting tail index to NaN."))
      Estimated_Tail_index[j] <- NaN
      next
    }
    
    # Sort y_z and compute the number of points M_z and k_z
    ysort_z = sort(y_z)
    M_z = length(ysort_z)
    k_z = floor(M_z * interm_lvl)
    
    # Check if k_z is 0
    if (k_z == 0) {
      warning(paste("Warning: k_z = 0 for point of interest z =", j, ". Setting tail index to NaN."))
      Estimated_Tail_index[j] <- NaN
      next
    }
    
    # local-Hill estimator for z
    Estimated_Tail_index[j] <- mean(log(ysort_z[(M_z - k_z + 1):M_z]) - log(ysort_z[M_z - k_z]))
  }
  
  return(Estimated_Tail_index)
}