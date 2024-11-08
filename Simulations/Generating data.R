
#seed = as.numeric(commandArgs(trailingOnly = TRUE))
seed = 1

dir.create(paste("Simulations/Generated data/seed_",seed,sep=""))

set.seed(seed)

n_max <- 10000 # max sample size

p_max = 30 # max dimension of X

U_max = runif(n_max)
X_max = matrix(runif(n_max*p_max),ncol=p_max)


source("Simulations/Models.R")

for(model in names(models)){
  name = models[[model]][["name"]]
  p = models[[model]][["p"]] ; q = models[[model]][["q"]] ; n = models[[model]][["n"]]
  B_0 = models[[model]][["B_0"]] ; B_1 = models[[model]][["B_1"]]
  xi = models[[model]][["xi"]] ; ell = models[[model]][["ell"]]
  nameFile = paste(
    name,
    "-p_",p,
    "-q_",q,
    "-n_",n,
    "-seed_",seed,
    ".RData",
    sep=""
  )
  
  X <- X_max[1:n,1:p]
  U <- U_max[1:n]
  y <- U^(-apply(t(t(B_0)%*%t(X)),1,xi))*apply(cbind(U,t(t(B_1)%*%t(X))),1,function(col) ell(col[1],col[-1]))
  
  data <- list(X=X,y=y,carac = list(name=name,xi=xi,q=q,B_0=B_0))
  
  save(data,file=paste("Simulations/Generated data/seed_",seed,"/",nameFile,sep=""))  
}


