#' Local_Hill estimator of the conditional tail-index with dimension reduction (cf definition 5).
#' The interm_lvl and bandwidth inputs correspond to alpha and h respectively in definition 5.
#'
#' @param X Matrix of covariates of dimension p.
#' @param y Response variable vector.
#' @param Z Matrix of points of interest.
#' @param B (p*q)-Matrix for dimension reduction from p to q.
#' @param interm_lvl Proportion of the largest y values used for estimation.
#' @param bandwidth Distance to select the observations closest to the point of interest z.
#' @param mink Minimum number of points should be used to estimate the tail-index.
#' @return Vector of estimated tail-index. 
#' @export


local_Hill = function(X, y, Z, B, interm_lvl = nrow(X)^(-0.3), bandwidth = nrow(X)^(-0.2/ncol(B))/2 , mink = 1) {
  
  # Check if inputs have correct dimensions and values
  if (ncol(X) != ncol(Z)) stop("Error: X and Z must have the same number of columns (dimension p).")
  if (nrow(B) != ncol(X)) stop("Error: B must have the same number of rows as columns in X (dimension p).")
  if (interm_lvl <= 0 || interm_lvl > 1) stop("Error: 'interm_lvl' must be between 0 and 1.")
  if (bandwidth <= 0) stop("Error: 'bandwidth' must be positive.")
  
  
  # Product storage for B' * X
  B.X <- lapply(1:ncol(B), function(d) c(B[, d] %*% t(X)))
  
  # Initialize vector for estimated tail indexes
  Estimated_Tail_index <- rep(NA, nrow(Z))
  
  for (j in 1:nrow(Z)) {
    # Point of interest z
    z = t(Z[j, ])
    
    # Selection of the observations closest to the point of interest z
    Indicator = which(abs(B.X[[1]] - as.numeric(B[, 1] %*% t(z))) < bandwidth)
    if (ncol(B) != 1) {
      for (d in 2:ncol(B)) {
        Indicator = intersect(Indicator, which(abs(B.X[[d]] - as.numeric(B[, d] %*% t(z))) < bandwidth))
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
    k_z = ceiling(M_z * interm_lvl)
    
    # Check if k_z is 0
    if (k_z < mink) {
      warning(paste("Warning: k_z <", mink ,"for point of interest z =", j, ". Setting tail index to NaN."))
      Estimated_Tail_index[j] <- NaN
      next
    }
    
    # local-Hill estimator for z
    Estimated_Tail_index[j] <- mean(log(ysort_z[(M_z - k_z + 1):M_z]) - log(ysort_z[M_z - k_z]))
  }
  
  return(Estimated_Tail_index)
}