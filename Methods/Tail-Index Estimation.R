#' Local_Hill estimator of the conditional tail-index with dimension reduction
#'
#' @param X Matrix of covariates of dimension p.
#' @param y Response variable vector.
#' @param Z Matrix of points of interest.
#' @param B (p*q)-Matrix for dimension reduction from p to q.
#' @param interm_lvl Proportion of the largest y values used for estimation.
#' @param bandwith Distance to select the observations closest to the point of interest z.
#' @return Vector of estimated tail-index. 
#' @export


local_Hill = function(X,y,Z,B,interm_lvl,bandwith){
  
  #Product storage $B^{\top} \cdot X$.
  B.X = list()
  for(d in 1:ncol(B)){
    B.X[[d]] = c(B[,d]%*%t(X))
  }
  
  # Tail index estimation for the points of interest Z.
  Estimated_Tail_index = c()
  
  for(j in 1:nrow(Z)){
    #Point of interest z
    z = t(Z[j,])
    
    #Selection of the observations closest to the point of interest z.
    Indic = which(abs(B.X[[1]] - as.numeric(B[,1]%*%t(z)))<bandwith)
    if(ncol(B)!=1){
      for(d in 2:ncol(B)){
        Indicd = which(abs(B.X[[d]] - as.numeric(B[,d]%*%t(z)))<bandwith)
        Indic = intersect(Indic,Indicd)
      }
    }

    Y_z = c(Y)[Indic]
    
    # local-Hill estimator for z
    Ysort_z = sort(Y_x)
    M_z = length(Ysort_z)
    k_z = M*interm_lvl
    
    Estimated_Tail_index[j] <-mean(log(Ysort_x[M_z-1:k_z+1])-log(Ysort_x[M_z-k_z]))
  }
  
  return(Estimated_Tail_index)
}