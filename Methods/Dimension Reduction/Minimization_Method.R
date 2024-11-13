#' Minimization method for objective function with (p*q)-matrix argument
#'
#' @param objective_fn Objective function to minimize.
#' @param dim vector c(p,q) of the matrix dimension.
#' @param num_init_evals Number of initial random evaluations for starting points.
#' @param num_starts Number of starting points to use for local optimization.
#' @return Estimated matrix `Bhat` representing the base of the CTI subspace.
#' @export


Minimization <- function(objective_fn,dim,num_init_evals, num_starts){

  p = dim[1]
  q = dim[2]

  # Function to compute the orthonormal matrix from a flat vector using QR decomposition
  orthonormalize <- function(A) {
    qr.Q(qr(matrix(A, nrow = p, ncol = q)))
  }

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



# Minimization <- function(objective_fn,dim,c_0,tol,num_init_evals,num_starts){
#   
#   p <- dim[1]
#   q <- dim[2]
#   
#   # Function to compute the orthonormal matrix from a flat vector using QR decomposition
#   orthonormalize <- function(A) {
#       qr.Q(qr(matrix(A, nrow = p, ncol = q)))
#     }
# 
#   # Generate Halton sequence for initial evaluations and reshape as matrices
#   init_evals <- randtoolbox::halton(num_init_evals, p * q)
# 
#   # Reshape each row of init_evals to a p*q matrix, apply orthonormalization, and flatten back to a vector
#   Inits <- t(apply(init_evals, 1, function(row) c(orthonormalize(row))))
# 
#   # Compute objective function values for each initialized matrix
#   Value_inits <- apply(Inits, 1, objective_fn)
# 
# 
#   # Select the indices of the starting points with the lowest objective values
#   best_start <- order(Value_inits)[1]
# 
#   B_new <- matrix(Inits[best_start,],ncol=q)
#   c <- c_0
#   Dif = 1
#   
#   while(Dif > tol){
#     B_old = B_new
# 
#     fitness_c = function(B){
#       B = matrix(B,ncol=q)
#       objective_fn(B) + c*norm(t(B)%*%B - diag(1,q),"2")
#     }
# 
#     Min_c = psoptim(c(B_old),fitness_c,control = list(trace=1,maxit=100))
#     B_new = matrix(Min_c$par,ncol=q)
# 
#     c = 2*c
#     Dif = norm(B_old-B_new,"2")
#   }
#   return(B_new)
# }