#' row-wise Minimization Starting from a Given Matrix
#'
#' This function performs a row-wise minimization of an objective function,
#' starting from an initial matrix `B_start`. The optimization proceeds iteratively,
#' adjusting rows of rows while preserving orthonormality constraints on the matrix.
#'
#' @param B_start Initial matrix for the optimization process.
#' @param objective_fn Objective function to minimize. It should accept a matrix as input.
#' @param dim Vector `c(p, q)` specifying the dimensions of the input matrix, where `p` is the
#'   number of rows and `q` is the number of columns.
#' @param tol Convergence tolerance for the stopping criterion. Iterations stop when
#'   the difference between successive matrices is less than this value.
#' @return An orthonormalized matrix that minimizes the given objective function.
#'
Minimization_from_start <- function(B_start, objective_fn, dim, tol = 1e-2) {
  
  p <- dim[1]  # Number of rows in the matrix
  q <- dim[2]  # Number of columns in the matrix
  
  # Step 1: Orthonormalize the starting matrix
  B_new <- orthonormalize(B_start, q)
  
  # Initialize the stopping criterion
  eps <- 1
  
  # Step 2: Iterative optimization loop
  while (eps > tol) {
    
    # Save the current matrix for comparison later
    B_old <- B_new
    
    # Step 3: Optimize the matrix row by row
    for (i in 1:p) {
      
      cat(i,"th row", sep ="", "\n")  # Print the current row being optimized
      
      # Define the objective function for the current row
      min_fn <- function(row) {
        # Update the current row in the matrix
        B_new[i, ] <- row
        # Evaluate the objective function with the updated row
        objective_fn(B_new)
      }
      
      # Perform the optimization on the current row
      result <- optim(
        par = B_new[i, ],  # Initial row
        fn = min_fn,                                               # Objective function
        method = "SANN",                                           # Simulated Annealing
        control = list(maxit = 20, tmax = 30)                      # Optimization settings
      )
      
      # Update the row with the optimized result
      B_new[i, ] <- result$par
      
      # Re-orthonormalize the matrix after each row update
      B_new <- orthonormalize(B_new, q)
    }
    
    # Step 4: Update the stopping criterion
    eps <- norm(B_old%*%t(B_old) - B_new%*%t(B_new), "2")
  }
  
  # Step 5: Return the optimized matrix
  return(B_new)
}


#' Minimization method for objective function with (p*q)-matrix argument
#'
#' This function performs minimization of an objective function where the 
#' argument is a matrix of dimensions `p*q`. It initializes several starting 
#' points, evaluates the function, and uses a row-wise optimization strategy 
#' to identify the optimal matrix.
#'
#' @param objective_fn Objective function to minimize.
#' @param dim A vector `c(p, q)` specifying the dimensions of the input matrix.
#' @param num_init_evals Number of initial random evaluations to generate starting points.
#' @param num_starts Number of starting points used for optimization.
#' @param tol Tolerance level for convergence in the optimization algorithm.
#' @return Estimated matrix `Bhat`, representing the base of the CTI subspace.
#' @export
Minimization <- function(objective_fn, dim, num_init_evals, num_starts, tol = 1e-2) {
  
  p <- dim[1]  # Number of rows in the matrix (p)
  q <- dim[2]  # Number of columns in the matrix (q)
  
  # Step 1: Generate initial random matrices using a Halton sequence
  init_evals <- 2*randtoolbox::halton(num_init_evals, p * q)-1
  
  # Step 2: Compute the objective function for each initialized matrix
  value_inits <- apply(init_evals, 1, objective_fn)
  
  # Step 3: Initialize storage for results
  B_s <- vector("list", num_starts)   # List to store optimized matrices
  Value_s <- numeric(num_starts)     # Vector to store function values for each optimization
  
  # Step 4: Perform optimization from multiple starting points
  for (s in 1:num_starts) {
    # Select the starting matrix with the best (smallest) objective value so far
    B_start <- init_evals[which.min(value_inits),]
    
    # Optimize the objective function starting from `B_start`
    B_s[[s]] <- Minimization_from_start(
      B_start,       # Initial matrix for optimization
      objective_fn,  # Objective function to minimize
      dim,           # Matrix dimensions
      tol            # Convergence tolerance
    )
    
    # Compute the objective function value for the optimized matrix
    Value_s[s] <- objective_fn(B_s[[s]])
  }
  
  # Step 5: Select the best optimized matrix
  B <- B_s[[which.min(Value_s)]]
  return(B)
}