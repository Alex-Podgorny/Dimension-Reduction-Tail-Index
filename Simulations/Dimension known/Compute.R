# Loading package
library(randtoolbox)



# Loading methods functions
methods_func <- list.files(path = "Methods", pattern = "\\.R$", full.names = TRUE,recursive = TRUE)

sapply(methods_func, source)

#set seed
#seed = as.numeric(commandArgs(trailingOnly = TRUE))
seed = 1 

dir.create(paste("Simulations/Dimension known/Errors/seed_",seed,sep=""))

for (data_name in list.files(path = paste("Simulations/Generated data/seed_",seed,sep="")) ){
  
  nameFile <- paste("Errors_",data_name,sep="")
  
  load(paste("Simulations/Generated data/seed_",seed,"/",data_name,sep=""))
  
  X <- data$X
  y <- data$y
  p <- ncol(X)
  n <- length(y)
  q <- data$carac$q
  B_0 <- data$carac$B_0
  xi <- data$carac$xi
  
  epsilon <- (1 - 0.9^(1/p))/2
  X0 <- which(apply(X,1,function(x) min(min(x),min(abs(x-1)))) > epsilon)
  Grid_X0 <- randtoolbox::halton(10000, p) * (1-2*epsilon) + epsilon
  
  alpha_exposants <- c(0.2,0.3,0.4,0.5)
  h_exposants <- c(0.2,0.3,0.4,0.5)
  
  Matrix_errors_Bhat_CTI <- matrix(NA,nrow = length(alpha_exposants) , ncol = length(h_exposants))
  Matrix_errors_Bhat_TDR <- matrix(NA,nrow = length(alpha_exposants) , ncol = length(h_exposants))
  Matrix_errors_Bhat_T1 <- matrix(NA,nrow = length(alpha_exposants) , ncol = length(h_exposants))
  Matrix_errors_Bhat_T2 <- matrix(NA,nrow = length(alpha_exposants) , ncol = length(h_exposants))
  
  Matrix_errors_gamma_CTI <- matrix(NA,nrow = length(alpha_exposants) , ncol = length(h_exposants))
  Matrix_errors_gamma_TDR <- matrix(NA,nrow = length(alpha_exposants) , ncol = length(h_exposants))
  Matrix_errors_gamma_T1 <- matrix(NA,nrow = length(alpha_exposants) , ncol = length(h_exposants))
  Matrix_errors_gamma_T2 <- matrix(NA,nrow = length(alpha_exposants) , ncol = length(h_exposants))
  
  
  for(i in 1:length(alpha_exposants)){
    for(j in 1:length(h_exposants)){
      alpha <- n^(-alpha_exposants[i])
      b <- h_exposants[j]
      h <- n^(-b/q)/2
      
      
      Bhat_CTI <- CTI(X, y, X0, q, alpha, h)
      Matrix_errors_Bhat_CTI[i,j] <- norm(Bhat_CTI%*%t(Bhat_CTI) - B_0%*%t(B_0),"2")
      Matrix_errors_gamma_CTI[i,j] <- mean((local_Hill(X,y,Grid_X0,Bhat_CTI,alpha,h) - apply(t(t(B_0)%*%t(Grid_X0)),1,xi))^2)
      
      Bhat_TDR <- Gardes(X, y, X0, q, alpha, h)
      Matrix_errors_Bhat_TDR[i,j] <- norm(Bhat_TDR%*%t(Bhat_TDR) - B_0%*%t(B_0),"2")
      Matrix_errors_gamma_TDR[i,j] <- mean((local_Hill(X,y,Grid_X0,Bhat_TDR,alpha,h) - apply(t(t(B_0)%*%t(Grid_X0)),1,xi))^2)
      
      Bhat_T1 <- TIREX1(X, y, X0, q, alpha)
      Matrix_errors_Bhat_T1[i,j] <- norm(Bhat_T1%*%t(Bhat_T1) - B_0%*%t(B_0),"2")
      Matrix_errors_gamma_T1[i,j] <- mean((local_Hill(X,y,Grid_X0,Bhat_T1,alpha,h) - apply(t(t(B_0)%*%t(Grid_X0)),1,xi))^2)
      
      Bhat_T2 <- TIREX2(X, y, X0, q, alpha)
      Matrix_errors_Bhat_T21[i,j] <- norm(Bhat_T2%*%t(Bhat_T2) - B_0%*%t(B_0),"2")
      Matrix_errors_gamma_T2[i,j] <- mean((local_Hill(X,y,Grid_X0,Bhat_T2,alpha,h) - apply(t(t(B_0)%*%t(Grid_X0)),1,xi))^2)
      
    }
  }
    
  Errors <- list(
    CTI = list(Bhat = Matrix_errors_Bhat_CTI,gamma = Matrix_errors_gamma_CTI),
    Gardes = list(Bhat = Matrix_errors_Bhat_TDR,gamma = Matrix_errors_gamma_TDR),
    TIREX1 = list(Bhat = Matrix_errors_Bhat_T1,gamma = Matrix_errors_gamma_T1),
    TIREX2 = list(Bhat = Matrix_errors_Bhat_T2,gamma = Matrix_errors_gamma_T2)
  )
  
  save(Errors,file=paste("Simulations/Dimension known/Errors/seed_",seed,"/",nameFile,sep=""))
  
}