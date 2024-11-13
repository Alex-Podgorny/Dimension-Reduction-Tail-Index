#' Estimation of the base of the central TDR subspace (Gardes 2018)
#'
#' @param X Matrix of covariates of dimension p.
#' @param y Response variable vector.
#' @param X0 Indicator vector for points X in compact subset X0.
#' @param q Dimension of the CTI subspace.
#' @param interm_lvl Proportion of the largest y values used for estimation.
#' @param bandwidth Distance to select the observations closest to the point of interest z.
#' @param N0 sub_sample size
#' @return Estimated matrix `Bhat` representing the base of the central TDR subspace.
#' @export

Gardes = function(X,y,N0,q,interm_lvl,bandwidth,num_init_evals = 100, num_starts = 3){
  
  
  vu<-matrix(runif(p*(p-q),-1,1),ncol=(p-q)) # For build orthogonal matrix#
  W<-X[N0,] # Random subsample in X0#
  
  
  indvar<-matrix(NA,nrow=(p-q),ncol=2^(p-q))
  for (i in 1:(p-q)) {
    indvar[i,]<-rep(c(rep((2*i-1),2^(p-i-q)),rep((2*i),2^(p-i-q))),2^(i-1))
  }
  
  #Intermediate function
  foptimtmp<-function(B,u) {
    
    A<-matrix(NA,ncol=p,nrow=p)
    A[,1:q]<- qr.Q(qr(matrix(B, nrow = p, ncol = q)))
    for (i in (q+1):p) {
      A[,i]<-A[,1]+(vu[,(i-q)])^3
    }
    Aortho<- qr.Q(qr(A))
    
    Indnul<-(1:nrow(X))*apply(X%*%B==matrix(u,byrow=T,ncol=q,nrow=nrow(X)),1,prod)
    Indnul<-Indnul[Indnul!=0]
    Xnew<-X[-Indnul,]
    ynew<-y[-Indnul]
    
    
    indid = list()
    indi = c(abs((Xnew%*%B)[,1] - u[1])<bandwidth)
    if(q!=1){
      for(d in 2:q){
        indid[[d]] = c(abs((Xnew%*%B)[,d] - u[d])<bandwidth)
        indi = indi*indid[[d]]
      }
    }
    
    matfK<-indi
    
    y_u = c(y)[as.logical(indi)]
    ysort_u = sort(y_u)
    M = length(ysort_u)
    
    yn<- ysort_u[M-M*interm_lvl]
    
    resa<-rep(NA,2^(p-q))
    #indi<-abs(Xnew%*%B-u)<bandwidth
    indi2<-c(1:nrow(Xnew))*indi
    if (sum(indi2!=0)<=1) {resa<-rep(Inf,2^(p-q))} else {
      indi2<-indi2[indi2!=0]
      matpivot2<-apply(Xnew[indi2,],2,median)
      
      matind<-matrix(NA,nrow=nrow(Xnew),ncol=2*(p-1))
      
      for (j in seq(1,(2*(p-1)-1),by=2)) {
        matXA<-Xnew%*%Aortho[,(ceiling(j/2)+1)]
        matPX<-matpivot2%*%Aortho[,(ceiling(j/2)+1)]
        matind[,j]<-(matXA<c(matPX))
        matind[,(j+1)]<-(matXA>=c(matPX))
      }
      
      for (j in 1:2^(p-q)) {
        indi3<-indi
        for (g in 1:(p-q)) {
          ig<-indvar[g,j]
          indi3<-indi3*matind[,ig]
        }
        resa[j]<-sum((ynew>yn)*indi3*matfK+1e-6)/(interm_lvl*sum(indi3*matfK)+1e-6)-1
      }
    }
    
    resa
  }	
  
  # Define the objective function
  objective_fn <- function(B) {
    vecu<-W%*%B
    foptimtmp_B = function(u){foptimtmp(B,u)}
    E<-apply(vecu,1,foptimtmp_B)
    sum((apply(E,1,mean))^2)
  }
  
  # Minimization of the objective function (function from Minimization_Method.R)
  Bhat <- Minimization(objective_fn,c(ncol(X),q),num_init_evals,num_starts)
  
  return(Bhat)
  
}
