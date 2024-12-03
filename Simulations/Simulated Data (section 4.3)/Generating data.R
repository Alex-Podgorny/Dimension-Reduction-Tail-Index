
for (seed in 1:100){
  dir.create(paste("Simulations/Generated data/seed_",seed,sep=""))
  
  set.seed(seed)
  
  n_max <- 4000 # max sample size
  
  p_max = 30 # max dimension of X
  
  U_max = runif(n_max)
  X_max = matrix(runif(n_max*p_max),ncol=p_max)
  
  
  source("Simulations/Models.R")
  
  for(model in names(models)){
    p = models[[model]][["p"]] ; q = models[[model]][["q"]] ; n = models[[model]][["n"]]
    B_0 = models[[model]][["B_0"]] ; B_1 = models[[model]][["B_1"]]
    xi = models[[model]][["xi"]] ; ell = models[[model]][["ell"]]
    nameFile = paste(
      model,
      "-seed_",seed,
      ".RData",
      sep=""
    )
    
    X <- X_max[1:n,1:p]
    U <- U_max[1:n]
    y <- U^(-apply(t(t(B_0)%*%t(X)),1,xi))*apply(cbind(U,t(t(B_1)%*%t(X))),1,function(col) ell(col[1],col[2]))
    
    data <- list(X=X,y=y,carac = list(xi=xi,q=q,B_0=B_0))
    
    save(data,file=paste("Simulations/Generated data/seed_",seed,"/",nameFile,sep=""))  
  }
}




