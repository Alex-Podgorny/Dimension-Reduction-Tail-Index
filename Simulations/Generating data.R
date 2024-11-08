
#seed = as.numeric(commandArgs(trailingOnly = TRUE))
seed = 1

set.seed(seed)

n <- 1000 # sample size

pmax = 30

U = runif(n)
X_max = matrix(runif(n*pmax),ncol=pmax)


source("Simulations/Models.R")

for(model in names(models)){
  name = models[[model]][["name"]]
  p = models[[model]][["p"]] ; q = models[[model]][["q"]]
  B_0 = models[[model]][["B_0"]] ; B_1 = models[[model]][["B_1"]]
  xi = models[[model]][["xi"]] ; ell = models[[model]][["ell"]]
  nameFile = paste(
    "Data-",name,
    "-p_",p,
    "-q_",q,
    "-n_",n,
    "-seed_",seed,
    ".RData",
    sep=""
  )
  
  X <- X_max[,1:p]
  y <- U^(-apply(t(t(B_0)%*%t(X)),1,xi))*apply(cbind(U,t(t(B_1)%*%t(X))),1,function(col) ell(col[1],col[-1]))
  
  Sim <- list(X=X,y=y)
  
  save(Sim,file=paste("Simulations/Generated data/",nameFile,sep=""))  
}


