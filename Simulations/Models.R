
# Four models are defined for the distribution of the random pair (X, Y).
# The covariate vector X follows a uniform distribution on the space [0,1]^p with p in {4,30}.
# The variable Y is defined as Y = U^(-γ(x)) * ℓ(U^-1, x), where U is a standard uniform variable.

# Two slowly varying functions ℓ are considered:
# ℓ^(1)(u^-1, x) = [1 + exp(B_1^T * x - u^-1)]^(-1) and 
# ℓ^(2)(u^-1, x) = [exp(-u/2) * B_1^T * x]^(-1),
# where B_1 = (0,0,5,5,0,...,0)^T in R^p.

ell_1 = function(u,t){(1 + exp(t - u^(-1)))^(-1)}

ell_2 = function(u,t){exp(-u/2)*t}

# For the first three models, the CTI subspace has dimension q = 1, with basis B_0 = (2,1,0,...,0)^T / sqrt(5).
# Two functions γ(x) = ξ_{B_0}(B_0^T * x) are used, defined as:
# ξ_{B_0}^{(1)}(z) = 0.1 + 0.9 * (sqrt(5) * z / 3)^4 and
# ξ_{B_0}^{(2)}(z) = 0.1 + 0.9 * |cos(2z)|, for all z in R.
# For the fourth model, the CTI subspace has dimension q = 2, with basis
# B_0 = (e_1, e_2) in R^{p x 2}, where e_1 = (1,0,...,0)^T and e_2 = (0,1,0,...,0)^T.
# The function γ(x) is specified as :
# ξ_{B_0}^{(3)}(z) = 0.1 + 1.8 * [(e_1^T * z - 0.5)^2 + (e_2^T * z - 0.5)^2].


xi_1 = function(t){
  return(0.1+ 0.9*(t/(3/sqrt(5)))^4)
}

xi_2 = function(t){
  return(0.1 + 0.9*abs(cos(2*t)))
}

xi_3 = function(t){
  0.1 + 1.8*((t[1]-0.5)^2 + (t[2]-0.5)^2)
}



# MODELS

models = list(
  model1p4 = list(
    name = "Model1",
    p=4,
    q=1,
    n=2000,
    B_0 = matrix(c(2,1,0,0),nrow=4,ncol=1)/sqrt(5),
    B_1 = matrix(c(0,0,5,5),nrow=4,ncol=1),
    xi = xi_1,
    ell = ell_1
  ),

  model2p4 = list(
    name = "Model2",
    p=4,
    q=1,
    n=2000,
    B_0 = matrix(c(2,1,0,0),nrow=4,ncol=1)/sqrt(5),
    B_1 = matrix(c(0,0,5,5),nrow=4,ncol=1),
    xi = xi_1,
    ell = ell_2
  ),

  model3p4 = list(
    name = "Model3",
    p=4,
    q=1,
    n=2000,
    B_0 = matrix(c(2,1,0,0),nrow=4,ncol=1)/sqrt(5),
    B_1 = matrix(c(0,0,5,5),nrow=4,ncol=1),
    xi = xi_2,
    ell = ell_1
  ),

  model4p4 = list(
    name = "Model4",
    p=4,
    q=2,
    n=2000,
    B_0 = diag(1,4,2),
    B_1 = matrix(c(0,0,5,5),nrow=4,ncol=1),
    xi = xi_3,
    ell = ell_1
  ),

  model1p30 = list(
    name = "Model1",
    p=30,
    q=1,
    n=4000,
    B_0 = matrix(c(2,1,0,0,rep(0,26)),nrow=30,ncol=1)/sqrt(5),
    B_1 = matrix(c(0,0,5,5,rep(0,26)),nrow=30,ncol=1),
    xi = xi_1,
    ell = ell_1
  ),

  model2p30 = list(
      name = "Model2",
      p=30,
      q=1,
      n=4000,
      B_0 = matrix(c(2,1,0,0,rep(0,26)),nrow=30,ncol=1)/sqrt(5),
      B_1 = matrix(c(0,0,5,5,rep(0,26)),nrow=30,ncol=1),
      xi = xi_1,
      ell = ell_2
    ),

  model3p30  = list(
    name = "Model3",
    p=30,
    q=1,
    n=4000,
    B_0 = matrix(c(2,1,0,0,rep(0,26)),nrow=30,ncol=1)/sqrt(5),
    B_1 = matrix(c(0,0,5,5,rep(0,26)),nrow=30,ncol=1),
    xi = xi_2,
    ell = ell_1
  ),

  model4p30 = list(
    name = "Model4",
    p=30,
    q=2,
    n=4000,
    B_0 = diag(1,30,2),
    B_1 = matrix(c(0,0,5,5,rep(0,26)),nrow=30,ncol=1),
    xi = xi_3,
    ell = ell_1
  )
)

  


  
  
